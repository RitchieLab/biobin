#!/usr/bin/env python

import os
import re
from loki import loki_source


class Source_gwas(loki_source.Source):
	
	
	##################################################
	# source interface
	
	
	@classmethod
	def getVersionString(cls):
		return '2.5 (2016-09-19)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
	#	self.downloadFilesFromHTTP('www.genome.gov', {
	#		'gwascatalog.txt': '/admin/gwascatalog.txt',
	#	})
		self.downloadFilesFromHTTP('www.ebi.ac.uk', {
			'gwas_catalog_v1.0-associations.tsv' : '/gwas/api/search/downloads/full'
		}, alwaysDownload=True)
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# process gwas cataog
		# the catalog uses dbSNP positions from b132, which should already be 1-based
		self.log("processing GWAS catalog annotations ...")
		reRS = re.compile('rs([0-9]+)', re.I)
		reChrPos = re.compile('(?:^|[^_])chr([0-9XYMT]+)[:_]([0-9]+)', re.I)
		reSNP = re.compile('(?:^|[^_])(?:chr([0-9XYMT]+)[:_]([0-9]+)|rs([0-9]+))', re.I)
		listNone = [None]
		numInc = numInvalid = 0
		setGwas = set()
		if os.path.exists('gwas_catalog_v1.0-associations.tsv'):
			with open('gwas_catalog_v1.0-associations.tsv','rU') as gwasFile:
				header = gwasFile.next().rstrip()
				cols = list(w.strip() for w in header.decode('latin-1').split("\t"))
				try:
					colPubmedID = cols.index("PUBMEDID")
					colTrait = cols.index("DISEASE/TRAIT")
					colChm = cols.index("CHR_ID")
					colPos = cols.index("CHR_POS")
					colAlleles = cols.index("STRONGEST SNP-RISK ALLELE")
					colSNPs = cols.index("SNPS")
					colRAF = cols.index("RISK ALLELE FREQUENCY")
					colORBeta = cols.index("OR or BETA")
					col95CI = cols.index("95% CI (TEXT)")
				except ValueError as e:
					self.log(" ERROR\n")
					raise Exception("unrecognized file header: %s" % str(e))
				l = 1
				for line in gwasFile:
					l += 1
					line = line.rstrip("\r\n")
					words = list(w.strip() for w in line.decode('latin-1').split("\t"))
					if len(words) <= col95CI:
						# blank line at the end is normal
						if (len(words) > 1) or words[0]:
							numInc += 1
						continue
					elif (' x ' in words[colPos]) or (' x ' in words[colSNPs]):
						# GWAS interaction pairs are not yet supported in LOKI
						numInvalid += 1
						continue
					pubmedID = int(words[colPubmedID]) if words[colPubmedID] else None
					trait = words[colTrait]
					listChm = words[colChm].split(';') if words[colChm] else list()
					listPos = words[colPos].split(';') if words[colPos] else list()
					snps = words[colSNPs] if words[colAlleles].endswith('aplotype') else words[colAlleles]
					listSNPs = reSNP.findall(snps)
					riskAfreq = words[colRAF]
					orBeta = words[colORBeta]
					allele95ci = words[col95CI]
					if (len(listChm) == len(listPos) == 0) and (len(listSNPs) > 0):
						listChm = listPos = list(None for i in xrange(len(listSNPs)))
					if (len(listChm) == len(listPos)) and (len(listChm) > 0) and (len(listSNPs) == 0):
						listSNPs = list((None,None,None) for i in xrange(len(listChm)))
					if len(listChm) == len(listPos) == len(listSNPs):
						for i in xrange(len(listSNPs)):
							rs = int(listSNPs[i][2]) if listSNPs[i][2] else None
							chm = self._loki.chr_num.get(listChm[i]) or self._loki.chr_num.get(listSNPs[i][0])
							pos = long(listPos[i]) if listPos[i] else (long(listSNPs[i][1]) if listSNPs[i][1] else None)
							setGwas.add( (rs,chm,pos,trait,snps,orBeta,allele95ci,riskAfreq,pubmedID) )
					elif len(listChm) == len(listPos):
						for i in xrange(len(listChm)):
							rs = None
							chm = self._loki.chr_num.get(listChm[i])
							pos = long(listPos[i]) if listPos[i] else None
							setGwas.add( (rs,chm,pos,trait,snps,orBeta,allele95ci,riskAfreq,pubmedID) )
						for i in xrange(len(listSNPs)):
							rs = int(listSNPs[i][2]) if listSNPs[i][2] else None
							chm = self._loki.chr_num.get(listSNPs[i][0])
							pos = long(listSNPs[i][1]) if listSNPs[i][1] else None
							setGwas.add( (rs,chm,pos,trait,snps,orBeta,allele95ci,riskAfreq,pubmedID) )
					else:
						numInvalid += 1
				#foreach line
			#with gwasFile
		else:
			with open('gwascatalog.txt','rU') as gwasFile:
				header = gwasFile.next().rstrip()
				if header.startswith("Date Added to Catalog\tPUBMEDID\tFirst Author\tDate\tJournal\tLink\tStudy\tDisease/Trait\tInitial Sample Size\tReplication Sample Size\tRegion\tChr_id\tChr_pos\tReported Gene(s)\tMapped_gene\tUpstream_gene_id\tDownstream_gene_id\tSnp_gene_ids\tUpstream_gene_distance\tDownstream_gene_distance\tStrongest SNP-Risk Allele\tSNPs\tMerged\tSnp_id_current\tContext\tIntergenic\tRisk Allele Frequency\tp-Value\tPvalue_mlog\tp-Value (text)\tOR or beta\t95% CI (text)\t"): # "Platform [SNPs passing QC]\tCNV"
					pass
				elif header.startswith("Date Added to Catalog\tPUBMEDID\tFirst Author\tDate\tJournal\tLink\tStudy\tDisease/Trait\tInitial Sample Description\tReplication Sample Description\tRegion\tChr_id\tChr_pos\tReported Gene(s)\tMapped_gene\tUpstream_gene_id\tDownstream_gene_id\tSnp_gene_ids\tUpstream_gene_distance\tDownstream_gene_distance\tStrongest SNP-Risk Allele\tSNPs\tMerged\tSnp_id_current\tContext\tIntergenic\tRisk Allele Frequency\tp-Value\tPvalue_mlog\tp-Value (text)\tOR or beta\t95% CI (text)\t"): # "Platform [SNPs passing QC]\tCNV"
					pass
				else:
					self.log(" ERROR\n")
					raise Exception("unrecognized file header")
				for line in gwasFile:
					line = line.rstrip("\r\n")
					words = list(w.strip() for w in line.decode('latin-1').split("\t"))
					if len(words) <= 31:
						# blank line at the end is normal
						if (len(words) > 1) or words[0]:
							numInc += 1
						continue
					chm = self._loki.chr_num[words[11]] if (words[11] in self._loki.chr_num) else None
					pos = long(words[12]) if words[12] else None
					trait = words[7]
					snps = words[21] if words[20].endswith('aplotype') else words[20]
					rses = list(int(rs[2:]) for rs in reRS.findall(snps)) or listNone
					orBeta = words[30]
					allele95ci = words[31]
					riskAfreq = words[26]
					pubmedID = int(words[1]) if words[1] else None
					for rs in rses:
						setGwas.add( (rs,chm,pos,trait,snps,orBeta,allele95ci,riskAfreq,pubmedID) )
				#foreach line
			#with gwasFile
		#if path
		self.log(" OK: %d entries (%d incomplete, %d invalid)\n" % (len(setGwas),numInc,numInvalid))
		if setGwas:
			self.log("writing GWAS catalog annotations to the database ...")
			self.addGWASAnnotations(setGwas)
			self.log(" OK\n")
	#update()
	
#Source_gwas
