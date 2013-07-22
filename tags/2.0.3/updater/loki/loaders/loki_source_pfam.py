#!/usr/bin/env python

import collections
from loki import loki_source


class Source_pfam(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.1 (2013-05-02)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromFTP('ftp.sanger.ac.uk', {
			'pfamA.txt.gz':                      '/pub/databases/Pfam/current_release/database_files/pfamA.txt.gz',
			'pfamA_reg_full_significant.txt.gz': '/pub/databases/Pfam/current_release/database_files/pfamA_reg_full_significant.txt.gz',
			'pfamseq.txt.gz':                    '/pub/databases/Pfam/current_release/database_files/pfamseq.txt.gz',
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('pfam_id',       0),
			('proteinfamily', 0),
			('uniprot_pid',   1),
		])
		relationshipID = self.addRelationships([
			('',),
		])
		typeID = self.addTypes([
			('proteinfamily',),
			('gene',),
		])
		
		# process protein families
		self.log("processing protein families ...")
		pfamFile = self.zfile('pfamA.txt.gz') #TODO:context manager,iterator
		groupFam = collections.defaultdict(set)
		famAcc = {}
		famID = {}
		famName = {}
		famDesc = {}
		for line in pfamFile:
			words = line.decode('latin-1').split("\t",10)
			pfamNum = int(words[0].strip())
			pfamAcc = words[1].strip()
			pfamID = words[2].strip()
			name = words[4].strip()
			group = words[8].strip()
			desc = words[9].strip()
			
			groupFam[group].add(pfamNum)
			famAcc[pfamNum] = pfamAcc
			famID[pfamNum] = pfamID
			famName[pfamNum] = name
			famDesc[pfamNum] = desc
		numGroup = len(groupFam)
		numFam = len(famName)
		self.log(" OK: %d categories, %d families\n" % (numGroup,numFam))
		
		# store protein families
		self.log("writing protein families to the database ...")
		listGroup = groupFam.keys()
		listGID = self.addTypedGroups(typeID['proteinfamily'], ((group,"") for group in listGroup))
		groupGID = dict(zip(listGroup,listGID))
		listFam = famAcc.keys()
		listGID = self.addTypedGroups(typeID['proteinfamily'], ((famName[fam],famDesc[fam]) for fam in listFam))
		famGID = dict(zip(listFam,listGID))
		self.log(" OK\n")
		
		# store protein family names
		self.log("writing protein family names to the database ...")
		self.addGroupNamespacedNames(namespaceID['pfam_id'], ((groupGID[group],group) for group in listGroup))
		self.addGroupNamespacedNames(namespaceID['pfam_id'], ((famGID[fam],famAcc[fam]) for fam in listFam))
		self.addGroupNamespacedNames(namespaceID['proteinfamily'], ((famGID[fam],famID[fam]) for fam in listFam))
		self.addGroupNamespacedNames(namespaceID['proteinfamily'], ((famGID[fam],famName[fam]) for fam in listFam))
		famName = famDesc = None
		self.log(" OK\n")
		
		# store protein family meta-group links
		self.log("writing protein family links to the database ...")
		for group in groupFam:
			self.addGroupRelationships( (famGID[fam],groupGID[group],relationshipID[''],None) for fam in groupFam[group] )
		groupFam = None
		self.log(" OK\n")
		
		# process protein identifiers
		self.log("processing protein identifiers ...")
		seqFile = self.zfile('pfamseq.txt.gz') #TODO:context manager,iterator
		proNames = dict()
		for line in seqFile:
			words = line.split("\t",10)
			proteinNum = int(words[0])
			uniprotID = words[1]
			uniprotAcc = words[2]
			species = words[9]
			
			if species == 'Homo sapiens (Human)':
				proNames[proteinNum] = (uniprotID,uniprotAcc)
		#foreach protein
		self.log(" OK: %d proteins\n" % (len(proNames),))
		
		# process associations
		self.log("processing protein associations ...")
		assocFile = self.zfile('pfamA_reg_full_significant.txt.gz') #TODO:context manager,iterator
		setAssoc = set()
		numAssoc = numID = 0
		for line in assocFile:
			words = line.split("\t",15)
			pfamNum = int(words[1])
			proteinNum = int(words[2])
			inFull = int(words[14])
			
			if (pfamNum in famGID) and (proteinNum in proNames) and inFull:
				numAssoc += 1
				numID += len(proNames[proteinNum])
				for name in proNames[proteinNum]:
					setAssoc.add( (famGID[pfamNum],numAssoc,name) )
			#if association is ok
		#foreach association
		self.log(" OK: %d associations (%d identifiers)\n" % (numAssoc,numID))
		
		# store gene associations
		self.log("writing gene associations to the database ...")
		self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID['uniprot_pid'], setAssoc)
		self.log(" OK\n")
	#update()
	
#Source_pfam
