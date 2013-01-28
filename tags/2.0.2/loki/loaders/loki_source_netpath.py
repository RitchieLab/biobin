#!/usr/bin/env python

import zipfile
from loki import loki_source


class Source_netpath(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-08-28)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromHTTP('www.netpath.org', {
		#	'NetPath_GeneReg_TSV.zip':  '/data/batch/NetPath_GeneReg_TSV.zip', #Last-Modified: Fri, 31 Oct 2008 17:00:16 GMT
			'NetPath_GeneReg_TSV1.zip': '/data/batch/NetPath_GeneReg_TSV1.zip', #Last-Modified: Sat, 03 Sep 2011 10:07:03 GMT
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('netpath_id', 0),
			('pathway',    0),
			('symbol',     0),
			('entrez_gid', 0),
		])
		typeID = self.addTypes([
			('pathway',),
			('gene',),
		])
		
		# process pathways and associations
		self.log("verifying archive ...")
		pathName = {}
		nsAssoc = {
			'symbol'     : set(),
			'entrez_gid' : set(),
		}
		numAssoc = 0
		with zipfile.ZipFile('NetPath_GeneReg_TSV1.zip','r') as pathZip:
			err = pathZip.testzip()
			if err:
				self.log(" ERROR\n")
				self.log("CRC failed for %s\n" % err)
				return False
			self.log(" OK\n")
			self.log("processing pathways ...")
			for info in pathZip.infolist():
				# there should be only one, but just in case..
				if info.filename == 'NetPath_Gene_regulation_all.txt':
					pathFile = pathZip.open(info,'rU')
					header = pathFile.next().rstrip()
					if header != "Gene regulation id	Pathway name	Pathway ID	Gene name	Entrez gene ID	Regulation	Experiment	PubMed ID":
						self.log(" ERROR\n")
						self.log("unrecognized file header in '%s': %s\n" % (info.filename,header))
						return False
					for line in pathFile:
						words = line.decode('latin-1').split("\t")
						pathway = words[1]
						pathID = words[2]
						gene = words[3].strip()
						entrezID = words[4]
						
						pathName[pathID] = pathway
						numAssoc += 1
						nsAssoc['entrez_gid'].add( (pathID,numAssoc,entrezID) )
						nsAssoc['symbol'].add( (pathID,numAssoc,gene) )
					#foreach line in pathFile
					pathFile.close()
				#if file is the one we want
			#foreach file in pathZip
		#with pathZip
		numPathways = len(pathName)
		numID = sum(len(nsAssoc[ns]) for ns in nsAssoc)
		self.log(" OK: %d pathways, %d associations (%d identifiers)\n" % (numPathways,numAssoc,numID))
		
		# store pathways
		self.log("writing pathways to the database ...")
		listPath = pathName.keys()
		listGID = self.addTypedGroups(typeID['pathway'], ((pathName[pathID],None) for pathID in listPath))
		pathGID = dict(zip(listPath,listGID))
		self.log(" OK\n")
		
		# store pathway names
		self.log("writing pathway names to the database ...")
		self.addGroupNamespacedNames(namespaceID['netpath_id'], ((pathGID[pathID],pathID) for pathID in listPath))
		self.addGroupNamespacedNames(namespaceID['pathway'], ((pathGID[pathID],pathName[pathID]) for pathID in listPath))
		self.log(" OK\n")
		
		# store gene associations
		self.log("writing gene associations to the database ...")
		for ns in nsAssoc:
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID[ns], ((pathGID[assoc[0]],assoc[1],assoc[2]) for assoc in nsAssoc[ns]))
		self.log(" OK\n")
	#update()
	
#Source_netpath
