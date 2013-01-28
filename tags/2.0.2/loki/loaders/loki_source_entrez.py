#!/usr/bin/env python

import collections
import re
from loki import loki_source


class Source_entrez(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a2 (2012-09-12)'
	#getVersionString()
	
	
	@classmethod
	def getOptions(cls):
		return {
			'locus-tags'    : "[yes|no]  --  include a gene's 'Locus Tag' as an alias (default: no)",
			'favor-primary' : "[yes|no]  --  reduce symbol ambiguity by favoring primary symbols (default: yes)",
			'favor-hist'    : "[yes|no]  --  reduce symbol ambiguity by favoring primary symbols (default: yes)",
		}
	#getOptions()
	
	
	def validateOptions(self, options):
		for o,v in options.iteritems():
			v = v.strip().lower()
			if o in ('locus-tags','favor-primary','favor-hist'):
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
		# download the latest source files
		self.downloadFilesFromFTP('ftp.ncbi.nih.gov', {
			'Homo_sapiens.gene_info.gz':       '/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
			'gene2refseq.gz':                  '/gene/DATA/gene2refseq.gz',
			'gene_history.gz':                 '/gene/DATA/gene_history.gz',
			'gene2ensembl.gz':                 '/gene/DATA/gene2ensembl.gz',
			'gene2unigene':                    '/gene/DATA/gene2unigene',
			'gene_refseq_uniprotkb_collab.gz': '/gene/DATA/gene_refseq_uniprotkb_collab.gz',
		})
		self.downloadFilesFromFTP('ftp.uniprot.org', {
			'HUMAN_9606_idmapping_selected.tab.gz': '/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz',
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		ldprofileID = self.addLDProfiles([
			('', 'no LD adjustment', None, None),
		])
		namespaceID = self.addNamespaces([
			('symbol',      0),
			('entrez_gid',  0),
			('refseq_gid',  0),
			('refseq_pid',  1),
			('ensembl_gid', 0),
			('ensembl_pid', 1),
			('hgnc_id',     0),
			('mim_id',      0),
			('hprd_id',     0),
			('vega_id',     0),
			('rgd_id',      0),
			('mirbase_id',  0),
			('unigene_gid', 0),
			('uniprot_gid', 0),
			('uniprot_pid', 1),
		])
		typeID = self.addTypes([
			('gene',),
		])
		
		nsNames = { ns:set() for ns in namespaceID }
		nsNameNames = { ns:set() for ns in namespaceID }
		numNames = numNameNames = numNameRefs = 0
		
		# process genes (no header!)
		self.log("processing genes ...")
		entrezGene = dict()
		entrezChm = dict()
		primaryEntrez = dict()
		xrefNS = {
			'Ensembl_G': 'ensembl_gid',
			'Ensembl_T': 'ensembl_gid',
			'Ensembl_P': 'ensembl_pid',
			'HGNC':      'hgnc_id',
			'MIM':       'mim_id',
			'HPRD':      'hprd_id',
			'Vega':      'vega_id',
			'RGD':       'rgd_id',
			'miRBase':   'mirbase_id',
		}
		geneFile = self.zfile('Homo_sapiens.gene_info.gz') #TODO:context manager,iterator
		for line in geneFile:
			# quickly filter out all non-9606 (human) taxonomies before taking the time to split()
			if line.startswith("9606\t"):
				words = line.rstrip().split("\t")
				entrezID = int(words[1])
				symbol = words[2]
				aliases = words[4].split("|") if words[4] != "-" else list()
				if options.get('locus-tags','no') == 'yes' and words[3] != "-":
					aliases.append(words[3])
				xrefs = words[5].split("|") if words[5] != "-" else list()
				chm = words[6]
				desc = words[8]
				
				entrezGene[entrezID] = (symbol,desc)
				entrezChm[entrezID] = chm
				if symbol not in primaryEntrez:
					primaryEntrez[symbol] = entrezID
				elif primaryEntrez[symbol] != entrezID:
					primaryEntrez[symbol] = False
				
				# entrezID as a name for itself looks funny here, but later on
				# we'll be translating the target entrezID to biopolymer_id and
				# adding more historical entrezID aliases
				nsNames['entrez_gid'].add( (entrezID,entrezID) )
				nsNames['symbol'].add( (entrezID,symbol) )
				for alias in aliases:
					nsNames['symbol'].add( (entrezID,alias) )
				for xref in xrefs:
					xrefDB,xrefID = xref.split(":",1)
					# turn ENSG/ENSP/ENST into Ensembl_X
					if xrefDB == "Ensembl" and xrefID.startswith("ENS") and len(xrefID) > 3:
						xrefDB = "Ensembl_%c" % xrefID[3]
					if xrefDB in xrefNS:
						nsNames[xrefNS[xrefDB]].add( (entrezID,xrefID) )
			#if taxonomy is 9606 (human)
		#foreach line in geneFile
		
		# delete any symbol alias which is also the primary name of exactly one other gene
		if options.get('favor-primary','yes') == 'yes':
			dupe = set()
			for alias in nsNames['symbol']:
				entrezID = alias[0]
				symbol = alias[1]
				if (symbol in primaryEntrez) and (primaryEntrez[symbol] != False) and (primaryEntrez[symbol] != entrezID):
					dupe.add(alias)
			nsNames['symbol'] -= dupe
			dupe = None
		#if favor-primary
		
		# print stats
		numGenes = len(entrezGene)
		numNames0 = numNames
		numNames = sum(len(nsNames[ns]) for ns in nsNames)
		self.log(" OK: %d genes, %d identifiers\n" % (numGenes,numNames-numNames0))
		
		# store genes
		self.log("writing genes to the database ...")
		listEntrez = entrezGene.keys()
		listBID = self.addTypedBiopolymers(typeID['gene'], (entrezGene[entrezID] for entrezID in listEntrez))
		entrezBID = dict(zip(listEntrez,listBID))
		numGenes = len(entrezBID)
		self.log(" OK: %d genes\n" % (numGenes))
		entrezGene = None
		
		# translate target entrezID to biopolymer_id in nsNames
		for ns in nsNames:
			names = set( (entrezBID[name[0]],name[1]) for name in nsNames[ns] if name[0] in entrezBID )
			nsNames[ns] = names
		numNames = sum(len(nsNames[ns]) for ns in nsNames)
		
		# process gene regions
		self.log("processing gene regions ...")
		reBuild = re.compile('GRCh([0-9]+)')
		buildGenes = collections.defaultdict(set)
		buildRegions = collections.defaultdict(set)
		setOrphan = set()
		setBadNC = set()
		setBadBuild = set()
		setBadChr = set()
		refseqBIDs = collections.defaultdict(set)
		regionFile = self.zfile('gene2refseq.gz') #TODO:context manager,iterator
		header = regionFile.next().rstrip()
		if header != "#Format: tax_id GeneID status RNA_nucleotide_accession.version RNA_nucleotide_gi protein_accession.version protein_gi genomic_nucleotide_accession.version genomic_nucleotide_gi start_position_on_the_genomic_accession end_position_on_the_genomic_accession orientation assembly (tab is used as a separator, pound sign - start of a comment)":
			self.log(" ERROR: unrecognized file header\n")
			self.log("%s\n" % header)
		else:
			for line in regionFile:
				# skip non-9606 (human) taxonomies before taking the time to split()
				if not line.startswith("9606\t"):
					continue
				
				# grab relevant columns
				words = line.split("\t")
				entrezID = int(words[1])
				rnaAcc = words[3].rsplit('.',1)[0] if words[3] != "-" else None
				proAcc = words[5].rsplit('.',1)[0] if words[5] != "-" else None
				genAcc = words[7].rsplit('.',1)[0] if words[7] != "-" else None
				posMin = long(words[9]) if words[9] != "-" else None
				posMax = long(words[10]) if words[10] != "-" else None
				build = reBuild.search(words[12].rstrip() if (len(words) > 12 and words[12] != "-") else '')
				
				# skip unrecognized IDs
				if entrezID not in entrezBID:
					setOrphan.add(entrezID)
					continue
				
				# store rna and protein sequence RefSeq IDs
				# (don't store genAcc, there's only one per chromosome)
				if rnaAcc:
					nsNames['refseq_gid'].add( (entrezBID[entrezID],rnaAcc) )
				if proAcc:
					nsNames['refseq_pid'].add( (entrezBID[entrezID],proAcc) )
					refseqBIDs[proAcc].add(entrezBID[entrezID])
				
				# skip non-whole-chromosome regions
				# (refseq accession types: http://www.ncbi.nlm.nih.gov/RefSeq/key.html)
				if not (genAcc and genAcc.startswith('NC_')):
					setBadNC.add(entrezID)
					continue
				elif not build:
					setBadBuild.add(entrezID)
					continue
				
				# skip chromosome mismatches
				if genAcc in ('NC_001807','NC_012920'): #TODO: avoid hardcoding this mapping
					chm = self._loki.chr_num.get('MT')
				else:
					chm = self._loki.chr_num.get(genAcc[3:].lstrip('0'))
				if not chm:
					setBadChr.add(entrezID)
					continue
				elif (entrezID in entrezChm) and (self._loki.chr_name[chm] not in entrezChm[entrezID].split('|')):
					# TODO: we're ignoring any gene region with an ambiguous chromosome
					# (gene_info says one thing, gene2refseq says another); is that right?
					#print "%s %s -> %s" % (entrezID,entrezChm[entrezID],self._loki.chr_name[chm])
					#100293744 X -> Y
					#100302657 3 -> 15
					#100418703 Y -> X
					#100507426 Y -> X
					setBadChr.add(entrezID)
					continue
				
				# store the region by build version number, so we can pick the majority build later
				buildGenes[build.group(1)].add(entrezID)
				buildRegions[build.group(1)].add( (entrezBID[entrezID],chm,posMin,posMax) )
			#foreach line in regionFile
			
			# identify majority build version
			grcBuild = max(buildRegions, key=lambda build: len(buildRegions[build]))
			setBadVers = set()
			for build,genes in buildGenes.iteritems():
				if build != grcBuild:
					setBadVers.update(genes)
			
			# print stats
			setBadVers.difference_update(buildGenes[build])
			setBadChr.difference_update(buildGenes[build], setBadVers)
			setBadBuild.difference_update(buildGenes[build], setBadVers, setBadChr)
			setBadNC.difference_update(buildGenes[build], setBadVers, setBadChr, setBadNC)
			numRegions = len(buildRegions[grcBuild])
			numGenes = len(buildGenes[grcBuild])
			numNames0 = numNames
			numNames = sum(len(nsNames[ns]) for ns in nsNames)
			self.log(" OK: %d regions (%d genes), %d identifiers\n" % (numRegions,numGenes,numNames-numNames0))
			self.logPush()
			if setOrphan:
				self.log("WARNING: %d regions for undefnied EntrezIDs\n" % (len(setOrphan)))
			if setBadNC:
				self.log("WARNING: %d genes not mapped to whole chromosome\n" % (len(setBadNC)))
			if setBadBuild:
				self.log("WARNING: %d genes not mapped to any GRCh build\n" % (len(setBadBuild)))
			if setBadVers:
				self.log("WARNING: %d genes mapped to GRCh build version other than %s\n" % (len(setBadVers),grcBuild))
			if setBadChr:
				self.log("WARNING: %d genes on mismatching chromosome\n" % (len(setBadChr)))
			self.logPop()
			entrezChm = setOrphan = setBadNC = setBadBuild = setBadChr = setBadVers = buildGenes = None
			
			# store gene regions
			self.log("writing gene regions to the database ...")
			numRegions = len(buildRegions[grcBuild])
			self.addBiopolymerLDProfileRegions(ldprofileID[''], buildRegions[grcBuild])
			self.log(" OK: %d regions\n" % (numRegions))
			buildRegions = None
		#if gene regions header ok
		
		# process historical gene names
		self.log("processing historical gene names ...")
		entrezUpdate = {}
		historyEntrez = {}
		histFile = self.zfile('gene_history.gz') #TODO:context manager,iterator
		header = histFile.next().rstrip()
		if header != "#Format: tax_id GeneID Discontinued_GeneID Discontinued_Symbol Discontinue_Date (tab is used as a separator, pound sign - start of a comment)":
			self.log(" ERROR: unrecognized file header\n")
			self.log("%s\n" % header)
		else:
			for line in histFile:
				# quickly filter out all non-9606 (human) taxonomies before taking the time to split()
				if line.startswith("9606\t"):
					words = line.split("\t")
					entrezID = int(words[1]) if words[1] != "-" else None
					oldEntrez = int(words[2]) if words[2] != "-" else None
					oldName = words[3] if words[3] != "-" else None
					
					if entrezID and entrezID in entrezBID:
						if oldEntrez and oldEntrez != entrezID:
							entrezUpdate[oldEntrez] = entrezID
							nsNames['entrez_gid'].add( (entrezBID[entrezID],oldEntrez) )
						if oldName and (oldName not in primaryEntrez or primaryEntrez[oldName] == False):
							if oldName not in historyEntrez:
								historyEntrez[oldName] = entrezID
							elif historyEntrez[oldName] != entrezID:
								historyEntrez[oldName] = False
							nsNames['symbol'].add( (entrezBID[entrezID],oldName) )
				#if taxonomy is 9606 (human)
			#foreach line in histFile
			
			# delete any symbol alias which is also the historical name of exactly one other gene
			if options.get('favor-hist','yes') == 'yes':
				dupe = set()
				for alias in nsNames['symbol']:
					entrezID = alias[0]
					symbol = alias[1]
					if (symbol in historyEntrez) and (historyEntrez[symbol] != False) and (historyEntrez[symbol] != entrezID):
						dupe.add(alias)
				nsNames['symbol'] -= dupe
				dupe = None
			#if favor-hist
			
			# print stats
			numNames0 = numNames
			numNames = sum(len(nsNames[ns]) for ns in nsNames)
			self.log(" OK: %d identifiers\n" % (numNames-numNames0))
		#if historical name header ok
		
		# process ensembl gene names
		self.log("processing ensembl gene names ...")
		ensFile = self.zfile('gene2ensembl.gz') #TODO:context manager,iterator
		header = ensFile.next().rstrip()
		if header != "#Format: tax_id GeneID Ensembl_gene_identifier RNA_nucleotide_accession.version Ensembl_rna_identifier protein_accession.version Ensembl_protein_identifier (tab is used as a separator, pound sign - start of a comment)":
			self.log(" ERROR: unrecognized file header\n")
			self.log("%s\n" % header)
		else:
			for line in ensFile:
				# quickly filter out all non-9606 (human) taxonomies before taking the time to split()
				if line.startswith("9606\t"):
					words = line.split("\t")
					entrezID = int(words[1])
					ensemblG = words[2] if words[2] != "-" else None
					ensemblT = words[4] if words[4] != "-" else None
					ensemblP = words[6] if words[6] != "-" else None
					
					if ensemblG or ensemblT or ensemblP:
						while entrezID and (entrezID in entrezUpdate):
							entrezID = entrezUpdate[entrezID]
						
						if entrezID and (entrezID in entrezBID):
							if ensemblG:
								nsNames['ensembl_gid'].add( (entrezBID[entrezID],ensemblG) )
							if ensemblT:
								nsNames['ensembl_gid'].add( (entrezBID[entrezID],ensemblT) )
							if ensemblP:
								nsNames['ensembl_pid'].add( (entrezBID[entrezID],ensemblP) )
				#if taxonomy is 9606 (human)
			#foreach line in ensFile
			
			# print stats
			numNames0 = numNames
			numNames = sum(len(nsNames[ns]) for ns in nsNames)
			self.log(" OK: %d identifiers\n" % (numNames-numNames0))
		#if ensembl name header ok
		
		# process unigene gene names
		self.log("processing unigene gene names ...")
		with open('gene2unigene','rU') as ugFile:
			header = ugFile.next().rstrip()
			if header != "#Format: GeneID UniGene_cluster (tab is used as a separator, pound sign - start of a comment)":
				self.log(" ERROR: unrecognized file header\n")
				self.log("%s\n" % header)
			else:
				for line in ugFile:
					words = line.rstrip().split("\t")
					entrezID = int(words[0]) if words[0] != "-" else None
					unigeneID = words[1] if words[1] != "-" else None
					
					while entrezID and (entrezID in entrezUpdate):
						entrezID = entrezUpdate[entrezID]
					
					# there will be lots of extraneous mappings for genes of other species
					if entrezID and (entrezID in entrezBID) and unigeneID:
						nsNames['unigene_gid'].add( (entrezBID[entrezID],unigeneID) )
				#foreach line in ugFile
				
				# print stats
				numNames0 = numNames
				numNames = sum(len(nsNames[ns]) for ns in nsNames)
				self.log(" OK: %d identifiers\n" % (numNames-numNames0))
			#if unigene name header ok
		#with ugFile
		
		if 0:
			# process uniprot gene names from entrez
			self.log("processing uniprot gene names ...")
			upFile = self.zfile('gene_refseq_uniprotkb_collab.gz') #TODO:context manager,iterator
			header = upFile.next().rstrip()
			if header != "#Format: NCBI_protein_accession UniProtKB_protein_accession (tab is used as a separator, pound sign - start of a comment)":
				self.log(" ERROR: unrecognized file header\n")
				self.log("%s\n" % header)
			else:
				for line in upFile:
					words = line.split("\t")
					proteinAcc = words[0].rsplit('.',1)[0] if words[0] != "-" else None
					uniprotAcc = words[1] if words[1] != "-" else None
					
					# there will be tons of identifiers missing from refseqBIDs because they're non-human
					if proteinAcc and (proteinAcc in refseqBIDs) and uniprotAcc:
						for biopolymerID in refseqBIDs[proteinAcc]:
							nsNames['uniprot_pid'].add( (biopolymerID,uniprotAcc) )
				#foreach line in upFile
					
				# print stats
				numNames0 = numNames
				numNames = sum(len(nsNames[ns]) for ns in nsNames)
				self.log(" OK: %d identifiers\n" % (numNames-numNames0))
			#if header ok
		
		# process uniprot gene names from uniprot (no header!)
		self.log("processing uniprot gene names ...")
		upFile = self.zfile('HUMAN_9606_idmapping_selected.tab.gz') #TODO:context manager,iterator
		""" /* ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README */
1. UniProtKB-AC
2. UniProtKB-ID
3. GeneID (EntrezGene)
4. RefSeq
5. GI
6. PDB
7. GO
8. IPI
9. UniRef100
10. UniRef90
11. UniRef50
12. UniParc
13. PIR
14. NCBI-taxon
15. MIM
16. UniGene
17. PubMed
18. EMBL
19. EMBL-CDS
20. Ensembl
21. Ensembl_TRS
22. Ensembl_PRO
23. Additional PubMed
"""
		for line in upFile:
			words = line.split("\t")
			uniprotAcc = words[0]
			uniprotID = words[1]
			found = False
			for word2 in words[2].split(';'):
				entrezID = int(word2.strip()) if word2 else None
				if entrezID and (entrezID in entrezBID):
					nsNameNames['uniprot_pid'].add( (namespaceID['entrez_gid'],entrezID,uniprotAcc) )
					nsNameNames['uniprot_gid'].add( (namespaceID['entrez_gid'],entrezID,uniprotID) )
					found = True
			#foreach entrezID mapping
			if not found:
				for word3 in words[3].split(';'):
					refseqID = word3.strip().split('.',1)[0] if word3 else None
					if refseqID:
						nsNameNames['uniprot_pid'].add( (namespaceID['refseq_pid'],refseqID,uniprotAcc) )
						nsNameNames['uniprot_pid'].add( (namespaceID['refseq_gid'],refseqID,uniprotAcc) )
						nsNameNames['uniprot_gid'].add( (namespaceID['refseq_pid'],refseqID,uniprotID) )
						nsNameNames['uniprot_gid'].add( (namespaceID['refseq_gid'],refseqID,uniprotID) )
				#foreach refseq mapping
				for word14 in words[14].split(';'):
					mimID = word14.strip() if word14 else None
					if mimID:
						nsNameNames['uniprot_pid'].add( (namespaceID['mim_id'],mimID,uniprotAcc) )
						nsNameNames['uniprot_gid'].add( (namespaceID['mim_id'],mimID,uniprotID) )
				#foreach mim mapping
				for word15 in words[15].split(';'):
					unigeneID = word15.strip() if word15 else None
					if unigeneID:
						nsNameNames['uniprot_pid'].add( (namespaceID['unigene_gid'],unigeneID,uniprotAcc) )
						nsNameNames['uniprot_gid'].add( (namespaceID['unigene_gid'],unigeneID,uniprotID) )
				#foreach mim mapping
				for word19 in words[19].split(';'):
					ensemblGID = word19.strip() if word19 else None
					if ensemblGID:
						nsNameNames['uniprot_pid'].add( (namespaceID['ensembl_gid'],ensemblGID,uniprotAcc) )
						nsNameNames['uniprot_gid'].add( (namespaceID['ensembl_gid'],ensemblGID,uniprotID) )
				#foreach ensG mapping
				for word20 in words[20].split(';'):
					ensemblTID = word20.strip() if word20 else None
					if ensemblTID:
						nsNameNames['uniprot_pid'].add( (namespaceID['ensembl_gid'],ensemblTID,uniprotAcc) )
						nsNameNames['uniprot_gid'].add( (namespaceID['ensembl_gid'],ensemblTID,uniprotID) )
				#foreach ensT mapping
				for word21 in words[21].split(';'):
					ensemblPID = word21.strip() if word21 else None
					if ensemblPID:
						nsNameNames['uniprot_pid'].add( (namespaceID['ensembl_pid'],ensemblPID,uniprotAcc) )
						nsNameNames['uniprot_gid'].add( (namespaceID['ensembl_pid'],ensemblPID,uniprotID) )
				#foreach ensP mapping
			#if no entrezID match
		#foreach line in upFile
		
		# print stats
		numNames0 = numNames
		numNames = sum(len(nsNames[ns]) for ns in nsNames)
		numNameNames0 = numNameNames
		numNameNames = sum(len(set(n[2] for n in nsNameNames[ns])) for ns in nsNameNames)
		numNameRefs0 = numNameRefs
		numNameRefs = sum(len(nsNameNames[ns]) for ns in nsNameNames)
		self.log(" OK: %d identifiers (%d references)\n" % (numNames-numNames0+numNameNames-numNameNames0,numNameRefs-numNameRefs0))
		
		# store gene names
		self.log("writing gene identifiers to the database ...")
		numNames = 0
		for ns in nsNames:
			if nsNames[ns]:
				numNames += len(nsNames[ns])
				self.addBiopolymerNamespacedNames(namespaceID[ns], nsNames[ns])
		self.log(" OK: %d identifiers\n" % (numNames,))
		nsNames = None
		
		# store gene names
		self.log("writing gene identifier references to the database ...")
		numNameNames = 0
		for ns in nsNameNames:
			if nsNameNames[ns]:
				numNameNames += len(nsNameNames[ns])
				self.addBiopolymerTypedNameNamespacedNames(typeID['gene'], namespaceID[ns], nsNameNames[ns])
		self.log(" OK: %d references\n" % (numNameNames,))
		nsNameNames = None
		
		# store source metadata
		self.setSourceBuilds(grcBuild, None)
	#update()
	
#Source_entrez
