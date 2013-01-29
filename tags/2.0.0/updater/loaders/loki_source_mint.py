#!/usr/bin/env python

import datetime
import os
import re
import loki_source


class Source_mint(loki_source.Source):
	
	
	##################################################
	# private class methods
	
	
	def _identifyLatestFilename(self, filenames):
		reFile = re.compile('^([0-9]+)-([0-9]+)-([0-9]+)-mint-human.txt$', re.IGNORECASE)
		bestdate = datetime.date.min
		bestfile = None
		for filename in filenames:
			match = reFile.match(filename)
			if match:
				filedate = datetime.date(int(match.group(1)), int(match.group(2)), int(match.group(3)))
				if filedate > bestdate:
					bestdate = filedate
					bestfile = filename
		#foreach filename
		return bestfile
	#_identifyLatestFilename()
	
	
	##################################################
	# source interface
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-08-30)'
	#getVersionString()
	
	
	def download(self, options):
		# define a callback to identify the latest *-mint-human.txt file
		def remFilesCallback(ftp):
			remFiles = {}
			
			path = '/pub/release/txt/current'
			ftp.cwd(path)
			bestfile = self._identifyLatestFilename(ftp.nlst())
			if bestfile:
				remFiles[bestfile] = '%s/%s' % (path,bestfile)
			
			return remFiles
		#remFilesCallback
		
		# download the latest source files
		self.downloadFilesFromFTP('mint.bio.uniroma2.it', remFilesCallback)
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('mint_id',     0),
			('symbol',      0),
			('entrez_gid',  0),
			('refseq_gid',  0),
			('refseq_pid',  1),
			('uniprot_pid', 1),
		])
		typeID = self.addTypes([
			('interaction',),
			('gene',),
		])
		
		# process interation groups
		self.log("processing interaction groups ...")
		mintDesc = dict()
		nsAssoc = {
			'symbol':      set(),
			'entrez_gid':  set(),
			'refseq_gid':  set(),
			'refseq_pid':  set(),
			'uniprot_pid': set(),
		}
		numAssoc = numID = 0
		with open(self._identifyLatestFilename(os.listdir('.')),'rU') as assocFile:
			header = assocFile.next().rstrip()
			if not header.startswith("ID interactors A (baits)\tID interactors B (preys)\tAlt. ID interactors A (baits)\tAlt. ID interactors B (preys)\tAlias(es) interactors A (baits)\tAlias(es) interactors B (preys)\tInteraction detection method(s)\tPublication 1st author(s)\tPublication Identifier(s)\tTaxid interactors A (baits)\tTaxid interactors B (preys)\tInteraction type(s)\tSource database(s)\tInteraction identifier(s)\t"): #Confidence value(s)\texpansion\tbiological roles A (baits)\tbiological role B\texperimental roles A (baits)\texperimental roles B (preys)\tinteractor types A (baits)\tinteractor types B (preys)\txrefs A (baits)\txrefs B (preys)\txrefs Interaction\tAnnotations A (baits)\tAnnotations B (preys)\tInteraction Annotations\tHost organism taxid\tparameters Interaction\tdataset\tCaution Interaction\tbinding sites A (baits)\tbinding sites B (preys)\tptms A (baits)\tptms B (preys)\tmutations A (baits)\tmutations B (preys)\tnegative\tinference\tcuration depth":
				self.log(" ERROR\n")
				self.log("unrecognized file header: %s\n" % header)
				return False
			xrefNS = {
				'entrezgene/locuslink': ('entrez_gid',),
				'refseq':               ('refseq_gid','refseq_pid'),
				'uniprotkb':            ('uniprot_pid',),
			}
			l = 0
			for line in assocFile:
				l += 1
				words = line.split('\t')
				genes = words[0].split(';')
				genes.extend(words[1].split(';'))
				aliases = words[4].split(';')
				aliases.extend(words[5].split(';'))
				method = words[6]
				taxes = words[9].split(';')
				taxes.extend(words[10].split(';'))
				labels = words[13].split('|')
				
				# identify interaction group label
				mint = None
				for label in labels:
					if label.startswith('mint:'):
						mint = label
						break
				mint = mint or "MINT-unlabeled-%d" % l
				mintDesc[mint] = method
				
				# identify interacting genes/proteins
				for n in xrange(0,len(taxes)):
					if taxes[n] == "taxid:9606(Homo sapiens)":
						numAssoc += 1
						# the "gene" is a helpful database cross-reference with a label indicating its type
						xrefDB,xrefID = genes[n].split(':',1)
						if xrefDB in xrefNS:
							numID += 1
							if xrefDB == 'refseq':
								xrefID = xrefID.rsplit('.',1)[0]
							elif xrefDB == 'uniprotkb':
								xrefID = xrefID.rsplit('-',1)[0]
							for ns in xrefNS[xrefDB]:
								nsAssoc[ns].add( (mint,numAssoc,xrefID) )
						# but the "alias" could be of any type and isn't identified,
						# so we'll store copies under each possible type
						# and find out later which one matches something
						numID += 1
						nsAssoc['symbol'].add( (mint,numAssoc,aliases[n]) )
						nsAssoc['refseq_gid'].add( (mint,numAssoc,aliases[n].rsplit('.',1)[0]) )
						nsAssoc['refseq_pid'].add( (mint,numAssoc,aliases[n].rsplit('.',1)[0]) )
						nsAssoc['uniprot_pid'].add( (mint,numAssoc,aliases[n].rsplit('-',1)[0]) )
					#if human
				#foreach interacting gene/protein
			#foreach line in assocFile
		#with assocFile
		self.log(" OK: %d associations (%d identifiers)\n" % (numAssoc,numID))
		
		# store interaction groups
		self.log("writing interaction groups to the database ...")
		listMint = mintDesc.keys()
		listGID = self.addTypedGroups(typeID['interaction'], ((mint,mintDesc[mint]) for mint in listMint))
		mintGID = dict(zip(listMint,listGID))
		self.log(" OK\n")
		
		# store interaction group names
		self.log("writing interaction group names to the database ...")
		self.addGroupNamespacedNames(namespaceID['mint_id'], ((mintGID[mint],mint) for mint in listMint))
		self.log(" OK\n")
		
		# store gene interactions
		self.log("writing gene interactions to the database ...")
		for ns in nsAssoc:
			self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID[ns], ((mintGID[a[0]],a[1],a[2]) for a in nsAssoc[ns]))
		self.log(" OK\n")
	#update()
	
#Source_mint
