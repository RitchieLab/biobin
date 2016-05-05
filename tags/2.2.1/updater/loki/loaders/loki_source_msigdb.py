#!/usr/bin/env python

from loki import loki_source

import collections
import re


class Source_msigdb(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0 (2015-08-14)'
	#getVersionString()
	
	
	@classmethod
	def getOptions(cls):
		return {
			'email'  : '<address>  --  email address registered with GSEA/MSigDB (required)',
		}
	#getOptions()
	
	
	def validateOptions(self, options):
		for o,v in options.iteritems():
			v = v.strip().lower()
			if o == 'email':
				if '@' in v:
					v = v.strip()
				else:
					return "invalid email address"
			else:
				return "unknown option '%s'" % o
			options[o] = v
		return True
	#validateOptions()
	
	
	def download(self, options):
		# post the user's registered email address to get a session cookie
		headers = self.getHTTPHeaders('www.broadinstitute.org', '/gsea/j_spring_security_check', {
			'j_username': options.get('email','anonymous@example.com'),
			'j_password': 'password',
			'login': 'login'
		})
		cookie = re.findall('(?:^|, )JSESSIONID=([0-9A-Za-z]+)(?:,|;|$)', headers['set-cookie'])[-1]
		
		# use the cookie to actually download the file
		self.downloadFilesFromHTTP('www.broadinstitute.org', {
			'msigdb.v5.0.entrez.gmt' : '/gsea/resources/msigdb/5.0/msigdb.v5.0.entrez.gmt'
		}, {'Cookie': 'JSESSIONID='+cookie})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		namespaceID = self.addNamespaces([
			('pathway',    0),
			('entrez_gid', 0),
		])
		typeID = self.addTypes([
			('pathway',),
			('gene',),
		])
		
		# process groups
		self.log("processing gene sets ...")
		setGroup = set()
		entrezAssoc = list()
		numAssoc = 0
		with open('msigdb.v5.0.entrez.gmt','rU') as gsFile:
			for line in gsFile:
				words = line.split("\t")
				label = words[0]
				setGroup.add(label)
				for i in xrange(2,len(words)):
					numAssoc += 1
					entrezAssoc.append( (label,numAssoc,words[i]) )
			#foreach line in gsFile
		#with gsFile
		listGroup = list(setGroup)
		setGroup = None
		self.log(" OK: %d gene sets, %d associations\n" % (len(listGroup),len(entrezAssoc)))
		
		# store groups
		self.log("writing gene sets to the database ...")
		listGID = self.addTypedGroups(typeID['pathway'], ((label,None) for label in listGroup))
		groupGID = dict(zip(listGroup,listGID))
		self.log(" OK\n")
		
		# store pathway names
		self.log("writing gene set names to the database ...")
		self.addGroupNamespacedNames(namespaceID['pathway'], ((gid,label) for label,gid in groupGID.iteritems()))
		self.log(" OK\n")
		
		# store gene associations
		self.log("writing gene associations to the database ...")
		self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID['entrez_gid'], ((groupGID[label],n,entrezID) for label,n,entrezID in entrezAssoc))
		self.log(" OK\n")
	#update()
	
#Source_msigdb
