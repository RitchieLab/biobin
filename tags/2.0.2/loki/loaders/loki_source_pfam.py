#!/usr/bin/env python

from loki import loki_source


class Source_pfam(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a2 (2012-09-13)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromFTP('ftp.sanger.ac.uk', {
			'pfamA.txt.gz':    '/pub/databases/Pfam/current_release/database_files/pfamA.txt.gz',
			'seq_info.txt.gz': '/pub/databases/Pfam/current_release/database_files/seq_info.txt.gz',
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
		groupAcc = {}
		accID = {}
		accName = {}
		accDesc = {}
		for line in pfamFile:
			words = line.decode('latin-1').split("\t")
			pfamAcc = words[1].strip()
			pfamID = words[2].strip()
			name = words[4].strip()
			group = words[8].strip()
			desc = words[9].strip()
			
			if group not in groupAcc:
				groupAcc[group] = set()
			groupAcc[group].add(pfamAcc)
			accID[pfamAcc] = pfamID
			accName[pfamAcc] = name
			accDesc[pfamAcc] = desc
		numGroup = len(groupAcc)
		numAcc = len(accName)
		self.log(" OK: %d categories, %d families\n" % (numGroup,numAcc))
		
		# store protein families
		self.log("writing protein families to the database ...")
		listGroup = groupAcc.keys()
		listGID = self.addTypedGroups(typeID['proteinfamily'], ((group,"") for group in listGroup))
		groupGID = dict(zip(listGroup,listGID))
		listAcc = accID.keys()
		listGID = self.addTypedGroups(typeID['proteinfamily'], ((accName[pfamAcc],accDesc[pfamAcc]) for pfamAcc in listAcc))
		accGID = dict(zip(listAcc,listGID))
		self.log(" OK\n")
		
		# store protein family names
		self.log("writing protein family names to the database ...")
		self.addGroupNamespacedNames(namespaceID['pfam_id'], ((groupGID[group],group) for group in listGroup))
		self.addGroupNamespacedNames(namespaceID['pfam_id'], ((accGID[pfamAcc],pfamAcc) for pfamAcc in listAcc))
		self.addGroupNamespacedNames(namespaceID['proteinfamily'], ((accGID[pfamAcc],accID[pfamAcc]) for pfamAcc in listAcc))
		self.addGroupNamespacedNames(namespaceID['proteinfamily'], ((accGID[pfamAcc],accName[pfamAcc]) for pfamAcc in listAcc))
		self.log(" OK\n")
		
		# store protein family meta-group links
		self.log("writing protein family links to the database ...")
		for group in groupAcc:
			self.addGroupRelationships( (accGID[pfamAcc],groupGID[group],relationshipID[''],None) for pfamAcc in groupAcc[group] )
		self.log(" OK\n")
		
		# process associations
		self.log("processing gene associations ...")
		assocFile = self.zfile('seq_info.txt.gz') #TODO:context manager,iterator
		setAssoc = set()
		numAssoc = numID = 0
		for line in assocFile:
			words = line.split("\t")
			if len(words) < 6:
				continue
			pfamAcc = words[0]
			uniprotAcc = words[5]
			uniprotID = words[6]
			species = words[8]
			
			if pfamAcc in accGID and species == 'Homo sapiens (Human)':
				numAssoc += 1
				numID += 2
				setAssoc.add( (accGID[pfamAcc],numAssoc,uniprotAcc) )
				setAssoc.add( (accGID[pfamAcc],numAssoc,uniprotID) )
			#if association is ok
		#foreach association
		self.log(" OK: %d associations (%d identifiers)\n" % (numAssoc,numID))
		
		# store gene associations
		self.log("writing gene associations to the database ...")
		self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID['uniprot_pid'], setAssoc)
		self.log(" OK\n")
	#update()
	
#Source_pfam
