#!/usr/bin/env python

#import collections
import itertools
import zipfile
import loki_source


class Source_reactome(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-08-28)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromHTTP('www.reactome.org', {
			'ReactomePathways.gmt.zip'                    : '/download/current/ReactomePathways.gmt.zip',
			'uniprot_2_pathways.txt'                      : '/download/current/uniprot_2_pathways.txt',
			'uniprot_2_pathways.stid.txt'                 : '/download/current/uniprot_2_pathways.stid.txt',
			'curated_and_inferred_uniprot_2_pathways.txt' : '/download/current/curated_and_inferred_uniprot_2_pathways.txt',
			'homo_sapiens.interactions.txt.gz'            : '/download/current/homo_sapiens.interactions.txt.gz',
		#	'gene_association.reactome'                   : '/download/current/gene_association.reactome',
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('symbol',       0),
			('entrez_gid',   0),
			('ensembl_gid',  0),
			('ensembl_pid',  1),
			('uniprot_pid',  1),
			('pathway',      0),
			('reactome_id',  0),
		])
		relationshipID = dict()
		typeID = self.addTypes([
			('gene',),
			('pathway',),
		])
		
		# initialize storage
		numAssoc = 0
		nsAssoc = {
			'symbol'      : { 'path':set(), 'react':set() },
			'entrez_gid'  : { 'path':set(), 'react':set() },
			'ensembl_gid' : { 'path':set(), 'react':set() },
			'ensembl_pid' : { 'path':set(), 'react':set() },
			'uniprot_pid' : { 'path':set(), 'react':set() },
		}
		pathReact = dict()
		reactPath = dict()
		#reactLinks = collections.defaultdict(set)
		
		# process pathways
		# <description>\t"Reactome Pathway"\t<symbol1>\t<symbol2>...
		self.log("verifying pathway archive ...")
		setPath = set()
		setID = set()
		with zipfile.ZipFile('ReactomePathways.gmt.zip','r') as pathZip:
			err = pathZip.testzip()
			if err:
				self.log(" ERROR\n")
				self.log("CRC failed for %s\n" % err)
				return False
			self.log(" OK\n")
			self.log("processing pathways ...")
			for info in pathZip.infolist():
				# there should be only one, but just in case..
				if info.filename == 'ReactomePathways.gmt':
					pathFile = pathZip.open(info,'rU')
					for line in pathFile:
						words = line.decode('latin-1').rstrip().split("\t")
						if line.startswith('#') or (len(words) < 3) or (words[1] != "Reactome Pathway"):
							continue
						path = words[0]
						
						if path not in pathReact:
							pathReact[path] = None
						for n in xrange(2, len(words)):
							numAssoc += 1
							setPath.add(path)
							setID.add(words[n])
							nsAssoc['symbol']['path'].add( (path,numAssoc,words[n]) )
						#foreach gene symbol
					#foreach line in pathFile
					pathFile.close()
				#if file is the one we want
			#foreach file in pathZip
			self.log(" OK: %d associations (%d pathways, %d identifiers)\n" % (numAssoc, len(setPath), len(setID)))
		#with pathZip
		
		# process uniprot mappings
		# <uniprot>\t<reactID>\t<description>\t<url>
		self.log("processing stable protein associations ...")
		newAssoc = 0
		setPath = set()
		setID = set()
		with open('uniprot_2_pathways.stid.txt', 'rU') as assocFile:
			for line in assocFile:
				words = line.decode('latin-1').split("\t")
				if line.startswith('#') or (len(words) < 3):
					continue
				uniprotPID = words[0]
				reactID = words[1]
				path = words[2]
				
				pathReact[path] = reactID
				reactPath[reactID] = path
				numAssoc += 1
				newAssoc += 1
				setPath.add(path)
				setID.add(uniprotPID)
				nsAssoc['uniprot_pid']['path'].add( (path,numAssoc,uniprotPID) )
			#foreach line in assocFile
		#with assocFile
		self.log(" OK: %d associations (%d pathways, %d identifiers)\n" % (newAssoc, len(setPath), len(setID)))
		
		# process (unstable?) uniprot mappings
		# <uniprot>\t"UniProt:"<uniprot>\t[<count>]: <description>; <description>; ...\t<url>
		self.log("processing additional protein associations ...")
		newAssoc = 0
		setPath = set()
		setID = set()
		with open('uniprot_2_pathways.txt', 'rU') as assocFile:
			for line in assocFile:
				words = line.decode('latin-1').split("\t")
				if line.startswith('#') or (len(words) < 3):
					continue
				uniprotPID = words[0]
				pathList = words[2]
				if pathList.startswith('['): # "[42 processes]: desc; desc; ..."
					pathList = pathList.split(']:',1)[1]
				
				for path in pathList.split(';'):
					path = path.strip()
					if path not in pathReact:
						pathReact[path] = None
					numAssoc += 1
					newAssoc += 1
					setPath.add(path)
					setID.add(uniprotPID)
					nsAssoc['uniprot_pid']['path'].add( (path,numAssoc,uniprotPID) )
				#foreach pathway description
			#foreach line in assocFile
		#with assocFile
		self.log(" OK: %d associations (%d pathways, %d identifiers)\n" % (newAssoc, len(setPath), len(setID)))
		
		# process curated/inferred uniprot mappings
		if 0:
			# <uniprot>\t"UniProt:"<uniprot>\t[<count>]: description; description; ...\t<url>\t<species>
			self.log("processing curated+inferred protein associations ...")
			newAssoc = 0
			setPath = set()
			setID = set()
			with open('curated_and_inferred_uniprot_2_pathways.txt', 'rU') as assocFile:
				for line in assocFile:
					words = line.decode('latin-1').split("\t")
					if line.startswith('#') or (len(words) < 5) or words[4].rstrip() != "Homo sapiens":
						continue
					uniprotPID = words[0]
					pathList = words[2]
					if pathList.startswith('['): # "[42 processes]: desc; desc; ..."
						pathList = pathList.split(']:',1)[1]
					
					for path in pathList.split(';'):
						path = path.strip()
						if path not in pathReact:
							pathReact[path] = None
						numAssoc += 1
						newAssoc += 1
						setPath.add(path)
						setID.add(uniprotPID)
						nsAssoc['uniprot_pid']['path'].add( (path,numAssoc,uniprotPID) )
				#foreach line in assocFile
			#with assocFile
			self.log(" OK: %d associations (%d pathways, %d identifiers)\n" % (newAssoc, len(setPath), len(setID)))
		
		# process interaction associations
		# <uniprot>\t<ensembl>\t<entrez>\t<uniprot>\t<ensembl>\t<entrez>\t<type>\t<reactID>["<->"<reactID>]
		self.log("processing protein interactions ...")
		newAssoc = 0
		setPath = set()
		setID = set()
		iaFile = self.zfile('homo_sapiens.interactions.txt.gz') #TODO:context manager,iterator
		for line in iaFile:
			words = line.decode('latin-1').split("\t")
			if line.startswith('#') or len(words) < 8:
				continue
			uniprotP1 = words[0][8:]  if words[0].startswith('UniProt:')     else None
			ensemblG1 = words[1][8:]  if words[1].startswith('ENSEMBL:ENSG') else None
			ensemblP1 = words[1][8:]  if words[1].startswith('ENSEMBL:ENSP') else None
			entrezG1  = words[2][12:] if words[2].startswith('Entrez Gene:') else None
			uniprotP2 = words[3][8:]  if words[3].startswith('UniProt:')     else None
			ensemblG2 = words[4][8:]  if words[4].startswith('ENSEMBL:ENSG') else None
			ensemblP2 = words[4][8:]  if words[4].startswith('ENSEMBL:ENSP') else None
			entrezG2  = words[5][12:] if words[5].startswith('Entrez Gene:') else None
			relationship = words[6]
			reactIDs = words[7].split('<->')
			reactID1 = reactIDs[0].split('.',1)[0]
			if reactID1 not in reactPath:
				reactPath[reactID1] = None
			setPath.add(reactID1)
			
			numAssoc += 1
			newAssoc += 1
			if uniprotP1:
				nsAssoc['uniprot_pid']['react'].add( (reactID1,numAssoc,uniprotP1) )
				setID.add(('uniprot_pid',uniprotP1))
			if ensemblG1:
				nsAssoc['ensembl_gid']['react'].add( (reactID1,numAssoc,ensemblG1) )
				setID.add(('ensembl_gid',ensemblG1))
			if ensemblP1:
				nsAssoc['ensembl_pid']['react'].add( (reactID1,numAssoc,ensemblP1) )
				setID.add(('ensembl_pid',ensemblP1))
			if entrezG1:
				nsAssoc['entrez_gid']['react'].add( (reactID1,numAssoc,entrezG1) )
				setID.add(('entrez_gid',entrezG1))
			
			numAssoc += 1
			newAssoc += 1
			if uniprotP2:
				nsAssoc['uniprot_pid']['react'].add( (reactID1,numAssoc,uniprotP2) )
				setID.add(('uniprot_pid',uniprotP2))
			if ensemblG2:
				nsAssoc['ensembl_gid']['react'].add( (reactID1,numAssoc,ensemblG2) )
				setID.add(('ensembl_gid',ensemblG2))
			if ensemblP2:
				nsAssoc['ensembl_pid']['react'].add( (reactID1,numAssoc,ensemblP2) )
				setID.add(('ensembl_pid',ensemblP2))
			if entrezG2:
				nsAssoc['entrez_gid']['react'].add( (reactID1,numAssoc,entrezG2) )
				setID.add(('entrez_gid',entrezG2))
			
			if len(reactIDs) > 1 and reactIDs[0] == reactIDs[1]:
				reactID2 = reactIDs[1].split('.',1)[0]
				if reactID2 not in reactPath:
					reactPath[reactID2] = None
				setPath.add(reactID2)
				
				numAssoc += 1
				newAssoc += 1
				if uniprotP1:
					nsAssoc['uniprot_pid']['react'].add( (reactID2,numAssoc,uniprotP1) )
					setID.add(('uniprot_pid',uniprotP1))
				if ensemblG1:
					nsAssoc['ensembl_gid']['react'].add( (reactID2,numAssoc,ensemblG1) )
					setID.add(('ensembl_gid',ensemblG1))
				if ensemblP1:
					nsAssoc['ensembl_pid']['react'].add( (reactID2,numAssoc,ensemblP1) )
					setID.add(('ensembl_pid',ensemblP1))
				if entrezG1:
					nsAssoc['entrez_gid']['react'].add( (reactID2,numAssoc,entrezG1) )
					setID.add(('entrez_gid',entrezG1))
				
				numAssoc += 1
				newAssoc += 1
				if uniprotP2:
					nsAssoc['uniprot_pid']['react'].add( (reactID2,numAssoc,uniprotP2) )
					setID.add(('uniprot_pid',uniprotP2))
				if ensemblG2:
					nsAssoc['ensembl_gid']['react'].add( (reactID2,numAssoc,ensemblG2) )
					setID.add(('ensembl_gid',ensemblG2))
				if ensemblP2:
					nsAssoc['ensembl_pid']['react'].add( (reactID2,numAssoc,ensemblP2) )
					setID.add(('ensembl_pid',ensemblP2))
				if entrezG2:
					nsAssoc['entrez_gid']['react'].add( (reactID2,numAssoc,entrezG2) )
					setID.add(('entrez_gid',entrezG2))
				
				#TODO: need the links?
				#if relationship not in relationshipID:
				#	relationshipID[relationship] = self.addRelationship(relationship)
				#reactLinks[reactID1].add( (reactID2,relationshipID[relationship]) )
				#reactLinks[reactID2].add( (reactID1,relationshipID[relationship]) )
			#if 2 reactions
		#foreach line in iaFile
		self.log(" OK: %d associations (%d pathways, %d identifiers)\n" % (newAssoc, len(setPath), len(setID)))
		
		# merge pathway-description and reactome-ID paths
		self.log("cross-referencing pathway mappings ...")
		listGroup = list()
		pathGroup = dict()
		reactGroup = dict()
		for path,react in itertools.chain(pathReact.iteritems(), ((path,react) for react,path in reactPath.iteritems())):
			if path and (path in pathGroup):
				if react:
					reactGroup[react] = pathGroup[path]
			elif react and (react in reactGroup):
				if path:
					pathGroup[path] = reactGroup[react]
			else:
				if path:
					pathGroup[path] = len(listGroup)
				if react:
					reactGroup[react] = len(listGroup)
				listGroup.append( (react,path) )
		#foreach pathway
		self.log(" OK: %d total pathways\n" % (len(listGroup),))
		
		# store pathways
		self.log("writing pathways to the database ...")
		listGID = self.addTypedGroups(typeID['pathway'], (((group[1] or group[0]),(group[0] if group[1] else None)) for group in listGroup))
		pathGID = { path:listGID[pathGroup[path]] for path in pathGroup }
		reactGID = { react:listGID[reactGroup[react]] for react in reactGroup }
		self.log(" OK\n")
		
		# store pathway names
		self.log("writing pathway names to the database ...")
		self.addGroupNamespacedNames(namespaceID['pathway'], ((gid,path) for path,gid in pathGID.iteritems()))
		self.addGroupNamespacedNames(namespaceID['reactome_id'], ((gid,react) for react,gid in reactGID.iteritems()))
		self.log(" OK\n")
		
		# store gene associations
		self.log("writing gene associations to the database ...")
		for ns in nsAssoc:
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID[ns], ((pathGID[assoc[0]],assoc[1],assoc[2]) for assoc in nsAssoc[ns]['path']))
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID[ns], ((reactGID[assoc[0]],assoc[1],assoc[2]) for assoc in nsAssoc[ns]['react']))
		self.log(" OK\n")
	#update()
	
#Source_reactome
