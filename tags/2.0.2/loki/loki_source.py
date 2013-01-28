#!/usr/bin/env python

import apsw
import datetime
import ftplib
import httplib
import itertools
import os
import sys
import time
import zlib

import loki_db


class Source(object):
	
	
	##################################################
	# constructor
	
	
	def __init__(self, lokidb):
		assert(isinstance(lokidb, loki_db.Database))
		assert(self.__class__.__name__.startswith('Source_'))
		self._loki = lokidb
		self._db = lokidb._db
		self._sourceID = self.addSource(self.getSourceName())
		assert(self._sourceID > 0)
	#__init__()
	
	
	##################################################
	# source interface
	
	
	@classmethod
	def getVersionString(cls):
		# when checked out from SVN, these $-delimited strings are magically kept updated
		rev = '$Revision$'.split()
		date = '$Date$'.split()
		stat = None
		
		if len(rev) > 2:
			version = 'r%s' % rev[1:2]
		else:
			stat = stat or os.stat(sys.modules[cls.__module__].__file__)
			version = '%s' % (stat.st_size,)
		
		if len(date) > 3:
			version += ' (%s %s)' % date[1:3]
		else:
			stat = stat or os.stat(sys.modules[cls.__module__].__file__)
			version += datetime.datetime.utcfromtimestamp(stat.st_mtime).strftime(' (%Y-%m-%d)' if (len(rev) > 2) else ' (%Y-%m-%d %H:%M:%S)')
		
		return version
	#getVersionString()
	
	
	@classmethod
	def getOptions(cls):
		return None
	#getOptions()
	
	
	def validateOptions(self, options):
		for o in options:
			return "unexpected option '%s'" % o
		return True
	#validateOptions()
	
	
	def download(self, options):
		raise Exception("invalid LOKI Source plugin: download() not implemented")
	#download()
	
	
	def update(self, options):
		raise Exception("invalid LOKI Source plugin: update() not implemented")
	#update()
	
	
	##################################################
	# context manager
	
	
	def __enter__(self):
		return self._loki.__enter__()
	#__enter__()
	
	
	def __exit__(self, excType, excVal, traceback):
		return self._loki.__exit__(excType, excVal, traceback)
	#__exit__()
	
	
	##################################################
	# logging
	
	
	def log(self, message=""):
		return self._loki.log(message)
	#log()
	
	
	def logPush(self, message=None):
		return self._loki.logPush(message)
	#logPush()
	
	
	def logPop(self, message=None):
		return self._loki.logPop(message)
	#logPop()
	
	
	##################################################
	# database update
	
	
	def prepareTableForUpdate(self, table):
		return self._loki.prepareTableForUpdate(table)
	#prepareTableUpdate()
	
	
	def prepareTableForQuery(self, table):
		return self._loki.prepareTableForQuery(table)
	#prepareTableQuery()
	
	
	##################################################
	# metadata management
	
	
	def addLDProfile(self, ldprofile, description=None, metric=None, value=None):
		return self.addLDProfiles([(ldprofile,description,metric,value)])[ldprofile]
	#addLDProfile()
	
	
	def addLDProfiles(self, ldprofiles):
		# ldprofiles=[ (ldprofile,description,metric,value), ... ]
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		ret = {}
		# use ABORT to avoid wasting autoincrements on existing rows,
		# and execute() to avoid bailing out of executemany() due to ABORT
		for l in ldprofiles:
			try:
				dbc.execute("INSERT OR ABORT INTO `db`.`ldprofile` (ldprofile,description,metric,value) VALUES (LOWER(?),?,LOWER(?),?); SELECT LAST_INSERT_ROWID()", l)
			except apsw.ConstraintError:
				dbc.execute("SELECT ldprofile_id FROM `db`.`ldprofile` WHERE ldprofile = LOWER(?)", l[0:1])
			for row in dbc:
				ret[l[0]] = row[0]
		return ret
	#addLDProfiles()
	
	
	def addNamespace(self, namespace, polygenic=0):
		return self.addNamespaces([(namespace,polygenic)])[namespace]
	#addNamespace()
	
	
	def addNamespaces(self, namespaces):
		# namespaces=[ (namespace,polygenic), ... ]
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		ret = {}
		# use ABORT to avoid wasting autoincrements on existing rows,
		# and execute() to avoid bailing out of executemany() due to ABORT
		for n in namespaces:
			try:
				dbc.execute("INSERT OR ABORT INTO `db`.`namespace` (namespace,polygenic) VALUES (LOWER(?),?); SELECT LAST_INSERT_ROWID()", n)
			except apsw.ConstraintError:
				dbc.execute("SELECT namespace_id FROM `db`.`namespace` WHERE namespace = LOWER(?)", n[0:1])
			for row in dbc:
				ret[n[0]] = row[0]
		return ret
	#addNamespaces()
	
	
	def addRelationship(self, relationship):
		return self.addRelationships([(relationship,)])[relationship]
	#addRelationship()
	
	
	def addRelationships(self, relationships):
		# relationships=[ (relationship,), ... ]
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		ret = {}
		# use ABORT to avoid wasting autoincrements on existing rows,
		# and execute() to avoid bailing out of executemany() due to ABORT
		for r in relationships:
			try:
				dbc.execute("INSERT OR ABORT INTO `db`.`relationship` (relationship) VALUES (LOWER(?)); SELECT LAST_INSERT_ROWID()", r)
			except apsw.ConstraintError:
				dbc.execute("SELECT relationship_id FROM `db`.`relationship` WHERE relationship = LOWER(?)", r[0:1])
			for row in dbc:
				ret[r[0]] = row[0]
		return ret
	#addRelationships()
	
	
	def addRole(self, role, description=None, coding=None, exon=None):
		return self.addRoles([(role,description,coding,exon)])[role]
	#addRole()
	
	
	def addRoles(self, roles):
		# roles=[ (role,description,coding,exon), ... ]
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		ret = {}
		# use ABORT to avoid wasting autoincrements on existing rows,
		# and execute() to avoid bailing out of executemany() due to ABORT
		for r in roles:
			try:
				dbc.execute("INSERT OR ABORT INTO `db`.`role` (role,description,coding,exon) VALUES (LOWER(?),?,?,?); SELECT LAST_INSERT_ROWID()", r)
			except apsw.ConstraintError:
				dbc.execute("SELECT role_id FROM `db`.`role` WHERE role = LOWER(?)", r[0:1])
			for row in dbc:
				ret[r[0]] = row[0]
		return ret
	#addRoles()
	
	
	def addSource(self, source):
		return self.addSources([(source,)])[source]
	#addSource()
	
	
	def addSources(self, sources):
		# sources=[ (source,), ... ]
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		ret = {}
		# use ABORT to avoid wasting autoincrements on existing rows,
		# and execute() to avoid bailing out of executemany() due to ABORT
		for s in sources:
			try:
				dbc.execute("INSERT OR ABORT INTO `db`.`source` (source) VALUES (LOWER(?)); SELECT LAST_INSERT_ROWID()", s)
			except apsw.ConstraintError:
				dbc.execute("SELECT source_id FROM `db`.`source` WHERE source = LOWER(?)", s[0:1])
			for row in dbc:
				ret[s[0]] = row[0]
		return ret
	#addSources()
	
		
	def addType(self, type):
		return self.addTypes([(type,)])[type]
	#addType()
	
	
	def addTypes(self, types):
		# types=[ (type,), ... ]
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		ret = {}
		# use ABORT to avoid wasting autoincrements on existing rows,
		# and execute() to avoid bailing out of executemany() due to ABORT
		for t in types:
			try:
				dbc.execute("INSERT OR ABORT INTO `db`.`type` (type) VALUES (LOWER(?)); SELECT LAST_INSERT_ROWID()", t)
			except apsw.ConstraintError:
				dbc.execute("SELECT type_id FROM `db`.`type` WHERE type = LOWER(?)", t[0:1])
			for row in dbc:
				ret[t[0]] = row[0]
		return ret
	#addTypes()
	
	
	def deleteAll(self):
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		tables = [
			'snp_merge', 'snp_locus', 'snp_entrez_role',
			'biopolymer', 'biopolymer_name', 'biopolymer_name_name', 'biopolymer_region',
			'group', 'group_name', 'group_group', 'group_biopolymer', 'group_member_name',
			'chain', 'chain_data',
		]
		for table in tables:
			dbc.execute("DELETE FROM `db`.`%s` WHERE source_id = %d" % (table,self.getSourceID()))
	#deleteAll()
	
	
	##################################################
	# source metadata management
	
	
	def getSourceName(self):
		return self.__class__.__name__[7:]
	#getSourceName()
	
	
	def getSourceID(self):
		return self._sourceID
	#getSourceID()
	
	
	def setSourceBuilds(self, grch=None, ucschg=None):
		self._loki.testDatabaseUpdate()
		sql = "UPDATE `db`.`source` SET grch = ?, ucschg = ?, current_ucschg = ? WHERE source_id = ?"
		self._db.cursor().execute(sql, (grch, ucschg, ucschg, self.getSourceID()))
	#setSourceBuilds()
	
	
	##################################################
	# snp data management
	
	
	def addSNPMerges(self, snpMerges):
		# snpMerges=[ (rsMerged,rsCurrent), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('snp_merge')
		sql = "INSERT OR IGNORE INTO `db`.`snp_merge` (rsMerged,rsCurrent,source_id) VALUES (?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, snpMerges)
	#addSNPMerges()
	
	
	def addSNPLoci(self, snpLoci):
		# snpLoci=[ (rs,chr,pos,validated), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('snp_locus')
		sql = "INSERT OR IGNORE INTO `db`.`snp_locus` (rs,chr,pos,validated,source_id) VALUES (?,?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, snpLoci)
	#addSNPLoci()
	
	
	def addChromosomeSNPLoci(self, chromosome, snpLoci):
		# snpLoci=[ (rs,pos,validated), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('snp_locus')
		sql = "INSERT OR IGNORE INTO `db`.`snp_locus` (rs,chr,pos,validated,source_id) VALUES (?,%d,?,?,%d)" % (chromosome,self.getSourceID(),)
		self._db.cursor().executemany(sql, snpLoci)
	#addChromosomeSNPLoci()
	
	
	def addSNPEntrezRoles(self, snpRoles):
		# snpRoles=[ (rs,entrez_id,role_id), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('snp_entrez_role')
		sql = "INSERT OR IGNORE INTO `db`.`snp_entrez_role` (rs,entrez_id,role_id,source_id) VALUES (?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, snpRoles)
	#addSNPEntrezRoles()
	
	
	##################################################
	# biopolymer data management
	
	
	def addBiopolymers(self, biopolymers):
		# biopolymers=[ (type_id,label,description), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer')
		sql = "INSERT INTO `db`.`biopolymer` (type_id,label,description,source_id) VALUES (?,?,?,%d); SELECT last_insert_rowid()" % (self.getSourceID(),)
		return [ row[0] for row in self._db.cursor().executemany(sql, biopolymers) ]
	#addBiopolymers()
	
	
	def addTypedBiopolymers(self, typeID, biopolymers):
		# biopolymers=[ (label,description), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer')
		sql = "INSERT INTO `db`.`biopolymer` (type_id,label,description,source_id) VALUES (%d,?,?,%d); SELECT last_insert_rowid()" % (typeID,self.getSourceID(),)
		return [ row[0] for row in self._db.cursor().executemany(sql, biopolymers) ]
	#addTypedBiopolymers()
	
	
	def addBiopolymerNames(self, biopolymerNames):
		# biopolymerNames=[ (biopolymer_id,namespace_id,name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer_name')
		sql = "INSERT OR IGNORE INTO `db`.`biopolymer_name` (biopolymer_id,namespace_id,name,source_id) VALUES (?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, biopolymerNames)
	#addBiopolymerNames()
	
	
	def addBiopolymerNamespacedNames(self, namespaceID, biopolymerNames):
		# biopolymerNames=[ (biopolymer_id,name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer_name')
		sql = "INSERT OR IGNORE INTO `db`.`biopolymer_name` (biopolymer_id,namespace_id,name,source_id) VALUES (?,%d,?,%d)" % (namespaceID,self.getSourceID(),)
		self._db.cursor().executemany(sql, biopolymerNames)
	#addBiopolymerNamespacedNames()
	
	
	def addBiopolymerNameNames(self, biopolymerNameNames):
		# biopolymerNameNames=[ (old_namespace_id,old_name,old_type_id,new_namespace_id,new_name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer_name_name')
		sql = "INSERT OR IGNORE INTO `db`.`biopolymer_name_name` (namespace_id,name,type_id,new_namespace_id,new_name,source_id) VALUES (?,?,?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, biopolymerNameNames)
	#addBiopolymerNameNames()
	
	
	def addBiopolymerTypedNameNamespacedNames(self, oldTypeID, newNamespaceID, biopolymerNameNames):
		# biopolymerNameNames=[ (old_namespace_id,old_name,new_name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer_name_name')
		sql = "INSERT OR IGNORE INTO `db`.`biopolymer_name_name` (namespace_id,name,type_id,new_namespace_id,new_name,source_id) VALUES (?,?,%d,%d,?,%d)" % (oldTypeID,newNamespaceID,self.getSourceID(),)
		self._db.cursor().executemany(sql, biopolymerNameNames)
	#addBiopolymerTypedNameNamespacedNames()
	
	
	def addBiopolymerRegions(self, biopolymerRegions):
		# biopolymerRegions=[ (biopolymer_id,ldprofile_id,chr,posMin,posMax), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer_region')
		sql = "INSERT OR IGNORE INTO `db`.`biopolymer_region` (biopolymer_id,ldprofile_id,chr,posMin,posMax,source_id) VALUES (?,?,?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, biopolymerRegions)
	#addBiopolymerRegions()
	
	
	def addBiopolymerLDProfileRegions(self, ldprofileID, biopolymerRegions):
		# biopolymerRegions=[ (biopolymer_id,chr,posMin,posMax), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('biopolymer_region')
		sql = "INSERT OR IGNORE INTO `db`.`biopolymer_region` (biopolymer_id,ldprofile_id,chr,posMin,posMax,source_id) VALUES (?,%d,?,?,?,%d)" % (ldprofileID,self.getSourceID(),)
		self._db.cursor().executemany(sql, biopolymerRegions)
	#addBiopolymerLDProfileRegions()
	
	
	##################################################
	# group data management
	
	
	def addGroups(self, groups):
		# groups=[ (type_id,label,description), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group')
		sql = "INSERT INTO `db`.`group` (type_id,label,description,source_id) VALUES (?,?,?,%d); SELECT last_insert_rowid()" % (self.getSourceID(),)
		return [ row[0] for row in self._db.cursor().executemany(sql, groups) ]
	#addGroups()
	
	
	def addTypedGroups(self, typeID, groups):
		# groups=[ (label,description), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group')
		sql = "INSERT INTO `db`.`group` (type_id,label,description,source_id) VALUES (%d,?,?,%d); SELECT last_insert_rowid()" % (typeID,self.getSourceID(),)
		return [ row[0] for row in self._db.cursor().executemany(sql, groups) ]
	#addTypedGroups()
	
	
	def addGroupNames(self, groupNames):
		# groupNames=[ (group_id,namespace_id,name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_name')
		sql = "INSERT OR IGNORE INTO `db`.`group_name` (group_id,namespace_id,name,source_id) VALUES (?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, groupNames)
	#addGroupNames()
	
	
	def addGroupNamespacedNames(self, namespaceID, groupNames):
		# groupNames=[ (group_id,name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_name')
		sql = "INSERT OR IGNORE INTO `db`.`group_name` (group_id,namespace_id,name,source_id) VALUES (?,%d,?,%d)" % (namespaceID,self.getSourceID(),)
		self._db.cursor().executemany(sql, groupNames)
	#addGroupNamespacedNames()
	
	
	def addGroupRelationships(self, groupRels):
		# groupRels=[ (group_id,related_group_id,relationship_id,contains), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_group')
		# we SHOULD be able to do (?1,?2,?3) and (?2,?1,?3) with the same 3 bindings for each execution,
		# but apsw or SQLite appears to treat the compound statement separately, so we have to copy the bindings
		sql = "INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?1,?2,?3,1,(CASE WHEN ?4 IS NULL THEN NULL WHEN ?4 > 0 THEN 1 WHEN ?4 < 0 THEN -1 ELSE 0 END),%d)" % (self.getSourceID(),)
		sql += ";INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?2,?1,?3,-1,(CASE WHEN ?4 IS NULL THEN NULL WHEN ?4 > 0 THEN -1 WHEN ?4 < 0 THEN 1 ELSE 0 END),%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, (2*gr for gr in groupRels))
	#addGroupRelationships()
	
	
	def addGroupParentRelationships(self, groupRels):
		# groupRels=[ (group_id,related_group_id,relationship_id), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_group')
		sql = "INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?1,?2,?3,1,1,%d)" % (self.getSourceID(),)
		sql += ";INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?2,?1,?3,-1,-1,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, (2*gr for gr in groupRels))
	#addGroupParentRelationships()
	
	
	def addGroupChildRelationships(self, groupRels):
		# groupRels=[ (group_id,related_group_id,relationship_id), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_group')
		sql = "INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?1,?2,?3,1,-1,%d)" % (self.getSourceID(),)
		sql += ";INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?2,?1,?3,-1,1,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, (2*gr for gr in groupRels))
	#addGroupChildRelationships()
	
	
	def addGroupSiblingRelationships(self, groupRels):
		# groupRels=[ (group_id,related_group_id,relationship_id), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_group')
		sql = "INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?1,?2,?3,1,0,%d)" % (self.getSourceID(),)
		sql += ";INSERT OR IGNORE INTO `db`.`group_group` (group_id,related_group_id,relationship_id,direction,contains,source_id)"
		sql += " VALUES (?2,?1,?3,-1,0,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, (2*gr for gr in groupRels))
	#addGroupSiblingRelationships()
	
	
	def addGroupBiopolymers(self, groupBiopolymers):
		# groupBiopolymers=[ (group_id,biopolymer_id), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_biopolymer')
		sql = "INSERT OR IGNORE INTO `db`.`group_biopolymer` (group_id,biopolymer_id,specificity,implication,quality,source_id) VALUES (?,?,100,100,100,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, groupBiopolymers)
	#addGroupBiopolymers()
	
	
	def addGroupMemberNames(self, groupMemberNames):
		# groupMemberNames=[ (group_id,member,type_id,namespace_id,name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_member_name')
		sql = "INSERT OR IGNORE INTO `db`.`group_member_name` (group_id,member,type_id,namespace_id,name,source_id) VALUES (?,?,?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, groupMemberNames)
	#addGroupMemberNames()
	
	
	def addGroupMemberTypedNamespacedNames(self, typeID, namespaceID, groupMemberNames):
		# groupMemberNames=[ (group_id,member,name), ... ]
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('group_member_name')
		sql = "INSERT OR IGNORE INTO `db`.`group_member_name` (group_id,member,type_id,namespace_id,name,source_id) VALUES (?,?,%d,%d,?,%d)" % (typeID,namespaceID,self.getSourceID(),)
		self._db.cursor().executemany(sql, groupMemberNames)
	#addGroupMemberTypedNamespacedNames()
	
	
	##################################################
	# liftover data management
	
	
	def addChains(self, old_ucschg, new_ucschg, chain_list):
		# chain_list=[ (score,old_chr,old_start,old_end,new_chr,new_start,new_end,is_forward), ... ]
		"""
		Adds all of the chains described in chain_list and returns the
		ids of the added chains.  The chain_list must be an iterable
		container of objects that can be inserted into the chain table
		"""
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('chain')
		sql = "INSERT INTO `db`.`chain` (score,old_ucschg,old_chr,old_start,old_end,new_ucschg,new_chr,new_start,new_end,is_fwd,source_id)"
		sql += " VALUES (?,%d,?,?,?,%d,?,?,?,?,%d); SELECT last_insert_rowid()" % (old_ucschg,new_ucschg,self.getSourceID())
		return [ row[0] for row in self._db.cursor().executemany(sql, chain_list) ]
	#addChains()
	
	
	def addChainData(self, chain_data_list):
		"""
		Adds all of the chain data into the chain data table
		"""
		self._loki.testDatabaseUpdate()
		self.prepareTableForUpdate('chain_data')
		sql = "INSERT INTO `db`.`chain_data` (chain_id,old_start,old_end,new_start,source_id) VALUES (?,?,?,?,%d)" % (self.getSourceID(),)
		self._db.cursor().executemany(sql, chain_data_list)
	#addChainData()
	
	
	##################################################
	# source utility methods
	
	
	def zfile(self, fileName, splitChar="\n", chunkSize=1*1024*1024):
		dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		with open(fileName,'rb') as filePtr:
			text = ""
			loop = True
			while loop:
				data = filePtr.read(chunkSize)
				if data:
					text += dc.decompress(data)
					data = None
				else:
					text += dc.flush()
					loop = False
				if text:
					lines = text.split(splitChar)
					i,x = 0,len(lines)-1
					text = lines[x]
					while i < x:
						yield lines[i]
						i += 1
					lines = None
			#while data remains
			if text:
				yield text
		#with fileName
	#zfile()
	
	
	def findConnectedComponents(self, neighbors):
		f = set()
		c = list()
		for v in neighbors:
			if v not in f:
				f.add(v)
				c.append(self._findConnectedComponents_recurse(neighbors, v, f, {v}))
		return c
	#findConnectedComponents()
	
	
	def _findConnectedComponents_recurse(self, n, v, f, c):
		for u in n[v]:
			if u not in f:
				f.add(u)
				c.add(u)
				self._findConnectedComponents_recurse(n, v, f, c)
		return c
	#_findConnectedComponents_recurse()
	
	
	def findEdgeDisjointCliques(self, neighbors):
		# neighbors = {'a':{'b','c'}, 'b':{'a'}, 'c':{'a'}, ...}
		# 'a' not in neighbors['a']
		# 'b' in neighbors['a'] => 'a' in neighbors['b']
		
		# clone neighbors so we can modify the local copy
		n = { v:set(neighbors[v]) for v in neighbors }
		c = list()
		
		while True:
			# prune isolated vertices and extract hanging pairs
			for v in n.keys():
				try:
					if len(n[v]) == 0:
						del n[v]
					elif len(n[v]) == 1:
						u, = n[v]
						n[v].add(v)
						c.append(n[v])
						del n[v]
						n[u].remove(v)
						if len(n[u]) == 0:
							del n[u]
				except KeyError:
					pass
			#foreach vertex
			
			# if nothing remains, we're done
			if len(n) == 0:
				return c
			
			# find maximal cliques on the remaining graph
			cliques = self.findMaximalCliques(n)
			
			# add disjoint cliques to the solution and remove the covered edges from the graph
			cliques.sort(key=len, reverse=True)
			for clique in cliques:
				ok = True
				for v in clique:
					if len(n[v] & clique) != len(clique) - 1:
						ok = False
						break
				if ok:
					c.append(clique)
					for v in clique:
						n[v] -= clique
			#foreach clique
		#loop
	#findEdgeDisjointCliques()
	
	
	def findMaximalCliques(self, neighbors):
		# neighbors = {'a':{'b','c'}, 'b':{'a'}, 'c':{'a'}, ...}
		# 'a' not in neighbors['a']
		# 'b' in neighbors['a'] => 'a' in neighbors['b']
		#
		# this implementation of the Bron-Kerbosch algorithm incorporates the
		# top-level degeneracy ordering described in:
		#   Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time
		#   David Eppstein, Maarten Loeffler, Darren Strash
		
		# build vertex-degree and degree-vertices maps
		vd = dict()
		dv = list()
		for v in neighbors:
			d = len(neighbors[v])
			vd[v] = d
			while len(dv) <= d:
				dv.append(set())
			dv[d].add(v)
		#foreach vertex
		
		# compute degeneracy ordering
		o = list()
		while len(dv) > 0:
			for dvSet in dv:
				try:
					v = dvSet.pop()
				except KeyError:
					continue
				o.append(v)
				vd[v] = None
				for u in neighbors[v]:
					if vd[u]:
						dv[vd[u]].remove(u)
						vd[u] -= 1
						dv[vd[u]].add(u)
				while len(dv) > 0 and len(dv[-1]) == 0:
					dv.pop()
				break
			#for dvSet in dv (until dvSet is non-empty)
		#while dv remains
		vd = dv = None
		
		# run first recursion layer in degeneracy order
		p = set(o)
		x = set()
		c = list()
		for v in o:
			self._findMaximalCliques_recurse({v}, p & neighbors[v], x & neighbors[v], neighbors, c)
			p.remove(v)
			x.add(v)
		return c
	#findMaximalCliques()
	
	
	def _findMaximalCliques_recurse(self, r, p, x, n, c):
		if len(p) == 0:
			if len(x) == 0:
				return c.append(r)
		else:
			# cursory tests yield best performance by choosing the pivot
			# arbitrarily from x first if x is not empty, else p; also tried
			# picking from p always, picking the pivot with highest degree,
			# and picking the pivot earliest in degeneracy order
			u = iter(x).next() if (len(x) > 0) else iter(p).next()
			for v in (p - n[u]):
				self._findMaximalCliques_recurse(r | {v}, p & n[v], x & n[v], n, c)
				p.remove(v)
				x.add(v)
	#_findMaximalCliques_recurse()
	
	
	def downloadFilesFromFTP(self, remHost, remFiles):
		# remFiles=function(ftp) or {'filename.ext':'/path/on/remote/host/to/filename.ext',...}
		
		# connect to source server
		self.log("connecting to FTP server %s ..." % remHost)
		ftp = ftplib.FTP(remHost)
		ftp.login() # anonymous
		self.log(" OK\n")
		
		# if remFiles is callable, let it identify the files it wants
		if hasattr(remFiles, '__call__'):
			self.log("locating current files ...")
			remFiles = remFiles(ftp)
			self.log(" OK\n")
		
		# check local file sizes and times, and identify all needed remote paths
		remDirs = set()
		remSize = {}
		remTime = {}
		locSize = {}
		locTime = {}
		for (locPath, remFile) in remFiles.iteritems():
			remDirs.add(remFile[0:remFile.rfind('/')])
			
			remSize[remFile] = None
			remTime[remFile] = None
			locSize[locPath] = None
			locTime[locPath] = None
			if os.path.exists(locPath):
				stat = os.stat(locPath)
				locSize[locPath] = long(stat.st_size)
				locTime[locPath] = datetime.datetime.fromtimestamp(stat.st_mtime)
		
		# define FTP directory list parser
		# unfortunately the FTP protocol doesn't specify an easily parse-able
		# format, but most servers return "ls -l"-ish space-delimited columns
		# (permissions) (?) (user) (group) (size) (month) (day) (year-or-time) (filename)
		now = datetime.datetime.utcnow()
		def ftpDirCB(rem_dir, line):
			words = line.split()
			remFn = rem_dir + "/" + words[8]
			if len(words) >= 9 and remFn in remSize:
				remSize[remFn] = long(words[4])
				timestamp = ' '.join(words[5:8])
				try:
					time = datetime.datetime.strptime(timestamp,'%b %d %Y')
				except ValueError:
					try:
						time = datetime.datetime.strptime("%s %d" % (timestamp,now.year),'%b %d %H:%M %Y')
					except ValueError:
						try:
							time = datetime.datetime.strptime("%s %d" % (timestamp,now.year-1),'%b %d %H:%M %Y')
						except ValueError:
							time = now
					if (
							(time.year == now.year and time.month > now.month) or
							(time.year == now.year and time.month == now.month and time.day > now.day)
					):
						time = time.replace(year=now.year-1)
				remTime[remFn] = time
		
		# check remote file sizes and times
		self.log("identifying changed files ...")
		for remDir in remDirs:
			ftp.dir(remDir, lambda x: ftpDirCB(remDir, x))
		self.log(" OK\n")
		
		# download files as needed
		self.logPush("downloading changed files ...\n")
		for locPath in sorted(remFiles.keys()):
			if remSize[remFiles[locPath]] == locSize[locPath] and remTime[remFiles[locPath]] <= locTime[locPath]:
				self.log("%s: up to date\n" % locPath)
			else:
				self.log("%s: downloading ..." % locPath)
				#TODO: download to temp file, then rename?
				with open(locPath, 'wb') as locFile:
					#ftp.cwd(remFiles[locPath][0:remFiles[locPath].rfind('/')])
					ftp.retrbinary('RETR '+remFiles[locPath], locFile.write)
				#TODO: verify file size and retry a few times if necessary
				self.log(" OK\n")
			modTime = time.mktime(remTime[remFiles[locPath]].utctimetuple())
			os.utime(locPath, (modTime,modTime))
		
		# disconnect from source server
		try:
			ftp.quit()
		except Exception:
			ftp.close()
		self.logPop("... OK\n")
	#downloadFilesFromFTP()
	
	
	def downloadFilesFromHTTP(self, remHost, remFiles):
		# remFiles={'filename.ext':'/path/on/remote/host/to/filename.ext',...}
		
		# check local file sizes and times
		remSize = {}
		remTime = {}
		locSize = {}
		locTime = {}
		for locPath in remFiles:
			remSize[locPath] = None
			remTime[locPath] = None
			locSize[locPath] = None
			locTime[locPath] = None
			if os.path.exists(locPath):
				stat = os.stat(locPath)
				locSize[locPath] = long(stat.st_size)
				locTime[locPath] = datetime.datetime.fromtimestamp(stat.st_mtime)
		
		# check remote file sizes and times
		self.log("identifying changed files ...")
		for locPath in remFiles:
			try:
				http = httplib.HTTPConnection(remHost)
				http.request('HEAD', remFiles[locPath])
				response = http.getresponse()
			except IOError as e:
				self.log(" ERROR: %s" % e)
				return False
			
			content_length = response.getheader('content-length')
			if content_length:
				remSize[locPath] = long(content_length)
			
			last_modified = response.getheader('last-modified')
			if last_modified:
				try:
					remTime[locPath] = datetime.datetime.strptime(last_modified,'%a, %d %b %Y %H:%M:%S %Z')
				except ValueError:
					remTime[locPath] = datetime.datetime.utcnow()
			
			http.close()
		self.log(" OK\n")
		
		# download files as needed
		self.logPush("downloading changed files ...\n")
		for locPath in sorted(remFiles.keys()):
			if remSize[locPath] and remSize[locPath] == locSize[locPath] and remTime[locPath] and remTime[locPath] <= locTime[locPath]:
				self.log("%s: up to date\n" % locPath)
			else:
				self.log("%s: downloading ..." % locPath)
				#TODO: download to temp file, then rename?
				with open(locPath, 'wb') as locFile:
					http = httplib.HTTPConnection(remHost)
					http.request('GET', remFiles[locPath])
					response = http.getresponse()
					while True:
						data = response.read()
						if not data:
							break
						locFile.write(data)
					http.close()
				self.log(" OK\n")
			if remTime[locPath]:
				modTime = time.mktime(remTime[locPath].utctimetuple())
				os.utime(locPath, (modTime,modTime))
		self.logPop("... OK\n")
	#downloadFilesFromHTTP()
	
	
#Source
