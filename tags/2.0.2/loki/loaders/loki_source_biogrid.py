#!/usr/bin/env python

import zipfile
from loki import loki_source


class Source_biogrid(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-08-30)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromHTTP('thebiogrid.org', {
			'BIOGRID-ORGANISM-LATEST.tab2.zip': '/downloads/archives/Latest%20Release/BIOGRID-ORGANISM-LATEST.tab2.zip',
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('biogrid_id', 0),
			('symbol',     0),
			('entrez_gid', 0),
		])
		typeID = self.addTypes([
			('interaction',),
			('gene',),
		])
		
		# process associations
		self.log("verifying archive file ...")
		pairLabels = dict()
		empty = tuple()
		with zipfile.ZipFile('BIOGRID-ORGANISM-LATEST.tab2.zip','r') as assocZip:
			err = assocZip.testzip()
			if err:
				self.log(" ERROR\n")
				self.log("CRC failed for %s\n" % err)
				return False
			self.log(" OK\n")
			self.log("processing gene interactions ...")
			for info in assocZip.infolist():
				if info.filename.find('Homo_sapiens') >= 0:
					assocFile = assocZip.open(info,'r')
					header = assocFile.next().rstrip()
					observedHeaders = {
						"#BioGRID Interaction ID	Entrez Gene Interactor A	Entrez Gene Interactor B	BioGRID ID Interactor A	BioGRID ID Interactor B	Systematic Name Interactor A	Systematic Name Interactor B	Official Symbol Interactor A	Official Symbol Interactor B	Synonymns Interactor A	Synonyms Interactor B	Experimental System	Experimental System Type	Author	Pubmed ID	Organism Interactor A	Organism Interactor B	Throughput	Score	Modification	Phenotypes	Qualifications	Tags	Source Database",
						"#BioGRID Interaction ID	Entrez Gene Interactor A	Entrez Gene Interactor B	BioGRID ID Interactor A	BioGRID ID Interactor B	Systematic Name Interactor A	Systematic Name Interactor B	Official Symbol Interactor A	Official Symbol Interactor B	Synonyms Interactor A	Synonyms Interactor B	Experimental System	Experimental System Type	Author	Pubmed ID	Organism Interactor A	Organism Interactor B	Throughput	Score	Modification	Phenotypes	Qualifications	Tags	Source Database",
					}
					if header not in observedHeaders:
						self.log(" ERROR\n")
						self.log("unrecognized file header in '%s': %s\n" % (info.filename,header))
						return False
					for line in assocFile:
						words = line.split("\t")
						bgID = int(words[0])
						entrezID1 = int(words[1])
						entrezID2 = int(words[2])
						syst1 = words[5] if words[5] != "-" else None
						syst2 = words[6] if words[6] != "-" else None
						gene1 = words[7]
						gene2 = words[8]
						aliases1 = words[9].split('|') if words[9] != "-" else empty
						aliases2 = words[10].split('|') if words[10] != "-" else empty
						tax1 = words[15]
						tax2 = words[16]
						
						if tax1 == '9606' and tax2 == '9606':
							member1 = (entrezID1, gene1, syst1) + tuple(aliases1)
							member2 = (entrezID2, gene2, syst2) + tuple(aliases2)
							if member1 != member2:
								pair = (member1,member2)
								if pair not in pairLabels:
									pairLabels[pair] = set()
								pairLabels[pair].add(bgID)
						#if interaction is ok
					#foreach line in assocFile
					assocFile.close()
				#if Homo_sapiens file
			#foreach file in assocZip
		#with assocZip
		numAssoc = len(pairLabels)
		numGene = len(set(pair[0] for pair in pairLabels) | set(pair[1] for pair in pairLabels))
		numName = sum(len(pairLabels[pair]) for pair in pairLabels)
		self.log(" OK: %d interactions (%d genes), %d pair identifiers\n" % (numAssoc,numGene,numName))
		
		# store interaction groups
		self.log("writing interaction pairs to the database ...")
		listPair = pairLabels.keys()
		listGID = self.addTypedGroups(typeID['interaction'], (("biogrid:%s" % min(pairLabels[pair]),None) for pair in listPair))
		pairGID = dict(zip(listPair,listGID))
		self.log(" OK\n")
		
		# store interaction labels
		listLabels = []
		for pair in listPair:
			listLabels.extend( (pairGID[pair],label) for label in pairLabels[pair] )
		self.log("writing interaction names to the database ...")
		self.addGroupNamespacedNames(namespaceID['biogrid_id'], listLabels)
		self.log(" OK\n")
		
		# store gene interactions
		self.log("writing gene interactions to the database ...")
		nsAssoc = {
			'symbol':     set(),
			'entrez_gid': set(),
		}
		numAssoc = 0
		for pair in pairLabels:
			numAssoc += 1
			nsAssoc['entrez_gid'].add( (pairGID[pair],numAssoc,pair[0][0]) )
			for n in xrange(1,len(pair[0])):
				nsAssoc['symbol'].add( (pairGID[pair],numAssoc,pair[0][n]) )
			
			numAssoc += 1
			nsAssoc['entrez_gid'].add( (pairGID[pair],numAssoc,pair[1][0]) )
			for n in xrange(1,len(pair[1])):
				nsAssoc['symbol'].add( (pairGID[pair],numAssoc,pair[1][n]) )
		for ns in nsAssoc:
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID[ns], nsAssoc[ns])
		self.log(" OK\n")
		
		# identify pseudo-pathways
		if 0: #TODO
			self.log("identifying implied networks ...")
			geneAssoc = dict()
			for pair in listPair:
				if pair[0] not in geneAssoc:
					geneAssoc[pair[0]] = set()
				geneAssoc[pair[0]].add(pair[1])
				if pair[1] not in geneAssoc:
					geneAssoc[pair[1]] = set()
				geneAssoc[pair[1]].add(pair[0])
			listPath = self.findMaximalCliques(geneAssoc)
			numAssoc = sum(len(path) for path in listPath)
			numGene = len(geneAssoc)
			numGroup = len(listPath)
			self.log(" OK: %d associations (%d genes in %d groups)\n" % (numAssoc,numGene,numGroup))
	#update()
	
#Source_biogrid
