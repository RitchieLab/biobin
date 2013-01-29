#!/usr/bin/env python

import loki_source


class Source_kegg(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-08-15)'
	#getVersionString()
	
	
	@classmethod
	def getOptions(cls):
		return {
			'api': '[soap|rest]  --  use the old SOAP or the new REST access API (default: soap)'
		}
	#getOptions()
	
	
	def validateOptions(self, options):
		for o,v in options.iteritems():
			if o == 'api':
				v = v.strip().lower()
				if 'soap'.startswith(v):
					v = 'soap'
				elif 'rest'.startswith(v):
					v = 'rest'
				else:
					return "api must be 'soap' or 'rest'"
				options[o] = v
			else:
				return "unknown option '%s'" % o
		return True
	#validateOptions()
	
	
	def download(self, options):
		if (options.get('api') == 'rest'):
			# download the latest source files
			self.downloadFilesFromHTTP('rest.kegg.jp', {
				'list-pathway-hsa':  '/list/pathway/hsa',
				'link-hsa-pathway':  '/link/hsa/pathway',
			})
		#if api==rest
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('kegg_id',    0),
			('pathway',    0),
			('entrez_gid', 0),
		])
		typeID = self.addTypes([
			('pathway',),
			('gene',),
		])
		
		if (options.get('api') == 'rest'):
			# process pathways
			self.log("processing pathways ...")
			pathName = {}
			with open('list-pathway-hsa','rU') as pathFile:
				for line in pathFile:
					words = line.split("\t")
					pathID = words[0]
					name = words[1].rstrip()
					if name.endswith(" - Homo sapiens (human)"):
						name = name[:-23]
					
					pathName[pathID] = name
				#foreach line in pathFile
			#with pathFile
			self.log(" OK: %d pathways\n" % (len(pathName),))
			
			# store pathways
			self.log("writing pathways to the database ...")
			listPath = pathName.keys()
			listGID = self.addTypedGroups(typeID['pathway'], ((pathName[pathID],None) for pathID in listPath))
			pathGID = dict(zip(listPath,listGID))
			self.log(" OK\n")
			
			# store pathway names
			self.log("writing pathway names to the database ...")
			self.addGroupNamespacedNames(namespaceID['kegg_id'], ((pathGID[pathID],pathID) for pathID in listPath))
			self.addGroupNamespacedNames(namespaceID['pathway'], ((pathGID[pathID],pathName[pathID]) for pathID in listPath))
			self.log(" OK\n")
			
			# process associations
			self.log("processing gene associations ...")
			entrezAssoc = set()
			numAssoc = 0
			with open('link-hsa-pathway','rU') as assocFile:
				for line in assocFile:
					words = line.split("\t")
					pathID = words[0]
					hsaGene = words[1].rstrip()
					
					if (pathID in pathGID) and (hsaGene.startswith("hsa:")):
						numAssoc += 1
						entrezAssoc.add( (pathGID[pathID],numAssoc,hsaGene[4:]) )
					#if pathway and gene are ok
				#foreach line in assocFile
			#with assocFile
			self.log(" OK: %d associations\n" % (numAssoc,))
			
			# store gene associations
			self.log("writing gene associations to the database ...")
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID['entrez_gid'], entrezAssoc)
			self.log(" OK\n")
		else: #default api==soap
			import suds.client
			
			# connect to SOAP/WSDL service
			self.log("connecting to KEGG data service ...")
			service = suds.client.Client('http://soap.genome.jp/KEGG.wsdl').service
			self.log(" OK\n")
			
			# fetch pathway list
			self.log("fetching pathways ...")
			pathName = {}
			for pathway in service.list_pathways('hsa'):
				pathID = pathway['entry_id'][0]
				name = pathway['definition'][0]
				if name.endswith(" - Homo sapiens (human)"):
					name = name[:-23]
				
				pathName[pathID] = name
			#foreach pathway
			self.log(" OK: %d pathways\n" % (len(pathName),))
			
			# store pathways
			self.log("writing pathways to the database ...")
			listPath = pathName.keys()
			listGID = self.addTypedGroups(typeID['pathway'], ((pathName[pathID],None) for pathID in listPath))
			pathGID = dict(zip(listPath,listGID))
			self.log(" OK\n")
			
			# store pathway names
			self.log("writing pathway names to the database ...")
			self.addGroupNamespacedNames(namespaceID['kegg_id'], ((pathGID[pathID],pathID) for pathID in listPath))
			self.addGroupNamespacedNames(namespaceID['pathway'], ((pathGID[pathID],pathName[pathID]) for pathID in listPath))
			self.log(" OK\n")
			
			# fetch genes for each pathway
			self.log("fetching gene associations ...")
			entrezAssoc = set()
			numAssoc = 0
			for pathID in listPath:
					for hsaGene in service.get_genes_by_pathway(pathID):
							numAssoc += 1
							entrezAssoc.add( (pathGID[pathID],numAssoc,hsaGene[4:]) )
					#foreach association
			#foreach pathway
			self.log(" OK: %d associations\n" % (numAssoc,))
			
			# store gene associations
			self.log("writing gene associations to the database ...")
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID['entrez_gid'], entrezAssoc)
			self.log(" OK\n")
		#if rest/soap API
	#update()
	
#Source_kegg
