#!/usr/bin/env python

import collections
import os
import re
import loki_source


class Source_dbsnp(loki_source.Source):
	
	
	##################################################
	# private class data
	
	
	_chmList = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','PAR','MT')
	
	
	##################################################
	# private class data
	
	
	def _identifyLatestSNPContig(self, filenames):
		reFile = re.compile('^b([0-9]+)_SNPContigLocusId(.*)\.bcp\.gz$', re.IGNORECASE)
		bestbuild = 0
		bestfile = None
		for filename in filenames:
			match = reFile.match(filename)
			if match:
				filebuild = int(match.group(1))
				if filebuild > bestbuild:
					bestbuild = filebuild
					bestfile = filename
		#foreach file in path
		return bestfile
	#_identifyLatestSNPContig()
	
	
	##################################################
	# source interface
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a2 (2012-09-04)'
	#getVersionString()
	
	
	@classmethod
	def getOptions(cls):
		return {
			'loci'   : '[all|validated]  --  store all or only validated SNP loci (default: all)',
			'merges' : '[yes|no]  --  process and store RS# merge history (default: yes)',
			'roles'  : '[yes|no]  --  process and store SNP roles (default: yes)',
		}
	#getOptions()
	
	
	def validateOptions(self, options):
		for o,v in options.iteritems():
			v = v.strip().lower()
			if o == 'loci':
				if 'all'.startswith(v):
					v = 'all'
				elif 'validated'.startswith(v):
					v = 'validated'
				else:
					return "loci must be 'all' or 'validated'"
			elif o in ('roles','merges'):
				if 'yes'.startswith(v):
					v = 'yes'
				elif 'no'.startswith(v):
					v = 'no'
				else:
					return "%s must be 'yes' or 'no'" % o
			else:
				return "unknown option '%s'" % o
			options[o] = v
		return True
	#validateOptions()
	
	
	def download(self, options):
		# define a callback to identify the latest SNPContigLocusId file
		def remFilesCallback(ftp):
			remFiles = dict()
			for chm in self._chmList:
				remFiles['chr_%s.txt.gz' % chm] = '/snp/organisms/human_9606/chr_rpts/chr_%s.txt.gz' % chm
			
			if options.get('merges',True):
				remFiles['RsMergeArch.bcp.gz'] = '/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz'
			
			if options.get('roles',True):
				remFiles['SnpFunctionCode.bcp.gz'] = '/snp/organisms/human_9606/database/shared_data/SnpFunctionCode.bcp.gz'
				path = '/snp/organisms/human_9606/database/organism_data'
				ftp.cwd(path)
				bestfile = self._identifyLatestSNPContig(ftp.nlst())
				if bestfile:
					remFiles[bestfile] = '%s/%s' % (path,bestfile)
			
			return remFiles
		#remFilesCallback
		
		# download the latest source files
		self.downloadFilesFromFTP('ftp.ncbi.nih.gov', remFilesCallback)
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# process merge report (no header!)
		if options.get('merges','yes') == 'yes':
			self.log("processing SNP merge records ...")
			mergeFile = self.zfile('RsMergeArch.bcp.gz') #TODO:context manager,iterator
			numMerge = 0
			listMerge = list()
			for line in mergeFile:
				words = line.split("\t")
				if len(words) <= 6:
					continue
				rsOld = long(words[0])
				#rsNew = long(words[1])
				rsCur = long(words[6])
				
				numMerge += 1
				listMerge.append( (rsOld,rsCur) )
				
				# write to the database after each million, to keep memory usage down
				if len(listMerge) >= 1000000:
					self.log(" %1.f million so far\n" % (numMerge/1000000.0)) #TODO: time estimate
					self.log("writing SNP merge records to the database ...")
					self.addSNPMerges(listMerge)
					listMerge = list()
					self.log(" OK\n")
					self.log("processing SNP merge records ...")
			#foreach line in mergeFile
			self.log(" OK: %d merged RS#s\n" % numMerge)
			
			# write any remaining records
			if listMerge:
				self.log("writing SNP merge records to the database ...")
				self.addSNPMerges(listMerge)
				self.log(" OK\n")
			listMerge = None
		#if merges
		
		# process SNP role function codes
		if options.get('roles','yes') == 'yes':
			""" /* from dbSNP_main_table.sql.gz */
CREATE TABLE [SnpFunctionCode]
(
[code] [tinyint] NOT NULL ,
[abbrev] [varchar](20) NOT NULL ,
[descrip] [varchar](255) NOT NULL ,
[create_time] [smalldatetime] NOT NULL ,
[top_level_class] [char](5) NOT NULL ,
[is_coding] [tinyint] NOT NULL ,
[is_exon] [bit] NULL ,
[var_prop_effect_code] [int] NULL ,
[var_prop_gene_loc_code] [int] NULL ,
[SO_id] [varchar](32) NULL
)
"""
			self.log("processing SNP role codes ...")
			roleID = dict()
			codeFile = self.zfile('SnpFunctionCode.bcp.gz')
			for line in codeFile:
				words = line.split('\t')
				code = int(words[0])
				name = words[1]
				desc = words[2]
				coding = int(words[5]) if (len(words) > 5 and words[5] != '') else None
				exon = int(words[6]) if (len(words) > 6 and words[6] != '') else None
				
				roleID[code] = self.addRole(name, desc, coding, exon)
			#foreach line in codeFile
			self.log(" OK: %d codes\n" % len(roleID))
			
			# process SNP roles
			""" /* from human_9606_table.sql.gz */
CREATE TABLE [b137_SNPContigLocusId]
(
[snp_id] [int] NULL ,
[contig_acc] [varchar](32) NOT NULL ,
[contig_ver] [tinyint] NULL ,
[asn_from] [int] NULL ,
[asn_to] [int] NULL ,
[locus_id] [int] NULL ,
[locus_symbol] [varchar](64) NULL ,
[mrna_acc] [varchar](32) NOT NULL ,
[mrna_ver] [smallint] NOT NULL ,
[protein_acc] [varchar](32) NULL ,
[protein_ver] [smallint] NULL ,
[fxn_class] [int] NULL ,
[reading_frame] [int] NULL ,
[allele] [varchar](255) NULL ,
[residue] [varchar](1000) NULL ,
[aa_position] [int] NULL ,
[build_id] [varchar](4) NOT NULL ,
[ctg_id] [int] NULL ,
[mrna_start] [int] NULL ,
[mrna_stop] [int] NULL ,
[codon] [varchar](1000) NULL ,
[protRes] [char](3) NULL ,
[contig_gi] [int] NULL ,
[mrna_gi] [int] NULL ,
[mrna_orien] [tinyint] NULL ,
[cp_mrna_ver] [int] NULL ,
[cp_mrna_gi] [int] NULL ,
[verComp] [int] NULL
)
"""
			self.log("processing SNP roles ...")
			listRole = list()
			numRole = 0
			setOrphan = set()
			numOrphan = 0
			funcFile = self.zfile(self._identifyLatestSNPContig(os.listdir('.')))
			for line in funcFile:
				words = line.split("\t")
				rs = long(words[0])
				entrez = int(words[5])
				#genesymbol = words[6]
				code = int(words[11])
				
				if code not in roleID:
					setOrphan.add(code)
					numOrphan += 1
				else:
					listRole.append( (rs,entrez,roleID[code]) )
					numRole += 1
				
				# write to the database after each 2.5 million, to keep memory usage down
				if len(listRole) >= 2500000:
					self.log(" %1.1f million so far\n" % (numRole/1000000.0)) #TODO: time estimate
					self.log("writing SNP roles to the database ...")
					self.addSNPEntrezRoles(listRole)
					listRole = list()
					self.log(" OK\n")
					self.log("processing SNP roles ...")
			#foreach line in funcFile
			self.log(" OK: %d roles\n" % (numRole,))
			self.logPush()
			if setOrphan:
				self.log("WARNING: %d roles (%d codes) unrecognized\n" % (numOrphan,len(setOrphan)))
			setOrphan = None
			self.logPop()
			
			# write any remaining records
			if listRole:
				self.log("writing SNP roles to the database ...")
				self.addSNPEntrezRoles(listRole)
				self.log(" OK\n")
			listRole = None
		#if roles
		
		# process chromosome report files
		grcBuild = None
		snpLociValid = (options.get('loci') == 'validated')
		for fileChm in self._chmList:
			self.log("processing chromosome %s SNPs ..." % fileChm)
			chmFile = self.zfile('chr_%s.txt.gz' % fileChm)
			
			# verify file headers
			header1 = chmFile.next().rstrip()
			chmFile.next()
			chmFile.next()
			header2 = chmFile.next().rstrip()
			header3 = chmFile.next().rstrip()
			chmFile.next()
			chmFile.next()
			if header1 != "dbSNP Chromosome Report":
				raise Exception("ERROR: unrecognized file header '%s'" % header1)
			if header2 != "rs#\tmap\tsnp\tchr\tctg\ttotal\tchr\tctg\tctg\tctg\tctg\tchr\tlocal\tavg\ts.e.\tmax\tvali-\tgeno-\tlink\torig\tupd":
				raise Exception("ERROR: unrecognized file subheader '%s'" % header2)
			if header3 != "\twgt\ttype\thits\thits\thits\t\tacc\tver\tID\tpos\tpos\tloci\thet\thet\tprob\tdated\ttypes\touts\tbuild\tbuild":
				raise Exception("ERROR: unrecognized file subheader '%s'" % header3)
			
			# process lines
			reBuild = re.compile('GRCh([0-9]+)')
			setSNP = set()
			setChrPos = collections.defaultdict(set)
			setBadBuild = set()
			setBadVers = set()
			setBadValid = set()
			setBadChr = set()
			for line in chmFile:
				words = line.split("\t")
				rs = words[0].strip()
				chm = words[6].strip()
				pos = words[11].strip()
				valid = 1 if (int(words[16]) > 0) else 0
				build = reBuild.search(words[21])
				
				if rs != '' and chm != '' and pos != '':
					rs = long(rs)
					pos = long(pos)
					if not build:
						setBadBuild.add(rs)
					elif grcBuild and grcBuild != build.group(1):
						setBadVers.add(rs)
					elif snpLociValid and not valid:
						setBadValid.add(rs)
					elif (fileChm != 'PAR') and (chm != fileChm):
						setBadChr.add(rs)
					elif (fileChm == 'PAR') and (chm != 'X') and (chm != 'Y'):
						setBadChr.add(rs)
					else:
						if not grcBuild:
							grcBuild = build.group(1)
						setSNP.add(rs)
						setChrPos[chm].add( (rs,pos,valid) )
				#if rs/chm/pos provided
			#foreach line in chmFile
			
			# print results
			setBadChr.difference_update(setSNP)
			setBadValid.difference_update(setSNP, setBadChr)
			setBadVers.difference_update(setSNP, setBadChr, setBadValid)
			setBadBuild.difference_update(setSNP, setBadChr, setBadValid, setBadVers)
			self.log(" OK: %d SNPs, %d loci\n" % (len(setSNP),sum(len(setChrPos[chm]) for chm in setChrPos)))
			self.logPush()
			if setBadBuild:
				self.log("WARNING: %d SNPs not mapped to any GRCh build\n" % (len(setBadBuild)))
			if setBadVers:
				self.log("WARNING: %d SNPs mapped to GRCh build version other than %s\n" % (len(setBadVers),grcBuild))
			if setBadValid:
				self.log("WARNING: %d SNPs not validated\n" % (len(setBadValid)))
			if setBadChr:
				self.log("WARNING: %d SNPs on mismatching chromosome\n" % (len(setBadChr)))
			self.logPop()
			
			# store data
			self.log("writing chromosome %s SNPs to the database ..." % fileChm)
			for chm in setChrPos:
				self.addChromosomeSNPLoci(self._loki.chr_num[chm], setChrPos[chm])
			setSNP = setChrPos = setSNP = setBadBuild = setBadVers = setBadValid = setBadChr = None
			self.log(" OK\n")
		#foreach chromosome
		
		# store source metadata
		self.setSourceBuilds(grcBuild, None)
	#update()
	
#Source_dbsnp
