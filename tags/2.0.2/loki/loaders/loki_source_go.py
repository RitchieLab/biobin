#!/usr/bin/env python

import re
from loki import loki_source


class Source_go(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a3 (2012-10-03)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromFTP('ftp.geneontology.org', {
			'gene_association.goa_human.gz': '/go/gene-associations/gene_association.goa_human.gz',
			'gene_ontology.1_2.obo':         '/go/ontology/obo_format_1_2/gene_ontology.1_2.obo',
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('go_id',       0),
			('ontology',    0),
			('symbol',      0),
			('uniprot_pid', 1),
		])
		relationshipID = self.addRelationships([
			('is_a',),
		])
		typeID = self.addTypes([
			('ontology',),
			('gene',),
		])
		
		# process ontology terms
		self.log("processing ontology terms ...")
		# file format specification: http://www.geneontology.org/GO.format.obo-1_2.shtml
		# correctly handling all the possible escape sequences and special cases
		# in the OBO spec would be somewhat involved, but the previous version
		# of biofilter used a much simpler approach which seemed to work okay in
		# practice, so we'll stick with that for now
		reTrailingEscape = re.compile('(?:^|[^\\\\])(?:\\\\\\\\)*\\\\$')
		empty = tuple()
		goName = {}
		goDef = {}
		goLinks = {}
		#goNS = {}
		#oboProps = {}
		curStanza = curID = curAnon = curObs = curName = curNS = curDef = curLinks = None
		with open('gene_ontology.1_2.obo','rU') as oboFile:
			while True:
				try:
					line = oboFile.next().rstrip()
					parts = line.split('!',1)[0].split(':',1)
					tag = parts[0].strip()
					val = parts[1].strip() if (len(parts) > 1) else None
				except StopIteration:
					line = False
				
				if line == False or tag.startswith('['):
					if (curStanza == 'Term') and curID and (not curAnon) and (not curObs):
						goName[curID] = curName
						goDef[curID] = curDef
						goLinks[curID] = curLinks or empty
				#		goNS[curID] = curNS or (oboProps['default-namespace'][-1] if ('default-namespace' in oboProps) else None)
					if line == False:
						break
					curStanza = tag[1:tag.index(']')]
					curID = curAnon = curObs = curName = curNS = curDef = curLinks = None
				#elif not curStanza:
				#	# before the first stanza, tag-value pairs are global file properties
				#	if tag not in oboProps:
				#		oboProps[tag] = []
				#	oboProps[tag].append(val)
				elif tag == 'id':
					curID = val
				elif tag == 'alt_id':
					pass #TODO
				elif tag == 'def':
					curDef = val
					if val.startswith('"'):
						curDef = ''
						words = val.split('"')
						for w in xrange(1,len(words)):
							curDef += words[w]
							if not reTrailingEscape.search(words[w]):
								break
				elif tag == 'is_anonymous':
					curAnon = (val.lower().split()[0] == 'true')
				elif tag == 'is_obsolete':
					curObs = (val.lower().split()[0] == 'true')
				elif tag == 'replaced_by':
					pass #TODO
				#elif tag == 'namespace':
				#	curNS = val
				elif tag == 'name':
					curName = val
				elif tag == 'synonym':
					pass #TODO
				elif tag == 'xref':
					pass #TODO
				elif tag == 'is_a':
					curLinks = curLinks or set()
					curLinks.add( (val.split()[0], relationshipID['is_a'], -1) )
				elif tag == 'relationship':
					curLinks = curLinks or set()
					words = val.split()
					if words[0] not in relationshipID:
						relationshipID[words[0]] = self.addRelationship(words[0])
					if words[0] == 'part_of':
						contains = -1
					elif words[0] in ('regulates','positively_regulates','negatively_regulates'):
						contains = 0
					else:
						contains = None
					curLinks.add( (words[1], relationshipID[words[0]], contains) )
			#foreach line
		#with oboFile
		numTerms = len(goName)
		numLinks = sum(len(goLinks[goID]) for goID in goLinks)
		self.log(" OK: %d terms, %d links\n" % (numTerms,numLinks))
		
		# store ontology terms
		self.log("writing ontology terms to the database ...")
		listGoID = goName.keys()
		listGID = self.addTypedGroups(typeID['ontology'], ((goName[goID],goDef[goID]) for goID in listGoID))
		goGID = dict(zip(listGoID,listGID))
		self.log(" OK\n")
		
		# store ontology term names
		self.log("writing ontology term names to the database ...")
		self.addGroupNamespacedNames(namespaceID['go_id'], ((goGID[goID],goID) for goID in listGoID))
		self.addGroupNamespacedNames(namespaceID['ontology'], ((goGID[goID],goName[goID]) for goID in listGoID))
		self.log(" OK\n")
		
		# store ontology term links
		self.log("writing ontology term relationships to the database ...")
		listLinks = []
		for goID in goLinks:
			for link in (goLinks[goID] or empty):
				if link[0] in goGID:
					listLinks.append( (goGID[goID],goGID[link[0]],link[1],link[2]) )
		self.addGroupRelationships(listLinks)
		self.log(" OK\n")
		
		# process gene associations
		self.log("processing gene associations ...")
		assocFile = self.zfile('gene_association.goa_human.gz') #TODO:context manager,iterator
		nsAssoc = {
			'uniprot_pid': set(),
			'symbol':      set()
		}
		numAssoc = numID = 0
		for line in assocFile:
			words = line.split('\t')
			if len(words) < 13:
				continue
			xrefDB = words[0]
			xrefID = words[1]
			gene = words[2]
			#assocType = words[3]
			goID = words[4]
			#reference = words[5]
			evidence = words[6]
			#withID = words[7]
			#goType = words[8]
			#desc = words[9]
			aliases = words[10].split('|')
			#xrefType = words[11]
			taxon = words[12]
			#updated = words[13]
			#assigner = words[14]
			#extensions = words[15].split('|')
			#xrefIDsplice = words[16]
			
			# TODO: why ignore IEA?
			if xrefDB == 'UniProtKB' and goID in goGID and evidence != 'IEA' and taxon == 'taxon:9606':
				numAssoc += 1
				numID += 2
				nsAssoc['uniprot_pid'].add( (goGID[goID],numAssoc,xrefID) )
				nsAssoc['symbol'].add( (goGID[goID],numAssoc,gene) )
				for alias in aliases:
					numID += 1
					# aliases might be either symbols or uniprot identifiers, so try them both ways
					nsAssoc['uniprot_pid'].add( (goGID[goID],numAssoc,alias) )
					nsAssoc['symbol'].add( (goGID[goID],numAssoc,alias) )
			#if association is ok
		#foreach association
		self.log(" OK: %d associations (%d identifiers)\n" % (numAssoc,numID))
		
		# store gene associations
		self.log("writing gene associations to the database ...")
		for ns in nsAssoc:
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID[ns], nsAssoc[ns])
		self.log(" OK\n")
	#update()
	
#Source_go
