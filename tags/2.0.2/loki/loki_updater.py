#!/usr/bin/env python

import collections
import hashlib
import os
import pkgutil

import loki_db
import loki_source
import loaders


class Updater(object):
	
	
	##################################################
	# constructor
	
	
	def __init__(self, lokidb, is_test=False):
		assert(isinstance(lokidb, loki_db.Database))
		self._is_test = is_test
		self._loki = lokidb
		self._db = lokidb._db
		self._sourceLoaders = None
		self._sourceClasses = dict()
		self._sourceObjects = dict()
		self._updating = False
		self._tablesUpdated = None
		self._tablesDeindexed = None
	#__init__()
	
	
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
	
	
	def flagTableUpdate(self, table):
		self._tablesUpdated.add(table)
	#flagTableUpdate()
	
	
	def prepareTableForUpdate(self, table):
		if self._updating:
			self.flagTableUpdate(table)
			if table not in self._tablesDeindexed:
				self._tablesDeindexed.add(table)
				self._loki.dropDatabaseIndecies(None, 'db', table)
	#prepareTableForUpdate()
	
	
	def prepareTableForQuery(self, table):
		if self._updating:
			if table in self._tablesDeindexed:
				self._tablesDeindexed.remove(table)
				self._loki.createDatabaseIndecies(None, 'db', table)
	#prepareTableForQuery()
	
	
	def findSourceModules(self):
		if self._sourceLoaders == None:
			self._sourceLoaders = {}
			loader_path = loaders.__path__
			if self._is_test:
				loader_path = [os.path.join(l, "test") for l in loaders.__path__]
			for srcImporter,srcModuleName,_ in pkgutil.iter_modules(loader_path):
				if srcModuleName.startswith('loki_source_'):
					self._sourceLoaders[srcModuleName[12:]] = srcImporter.find_module(srcModuleName)
	#findSourceModules()
	
	
	def getSourceModules(self):
		self.findSourceModules()
		return self._sourceLoaders.keys()
	#getSourceModules()
	
	
	def loadSourceModules(self, sources=None):
		self.findSourceModules()
		srcSet = set()
		for srcName in (set(sources) if sources else self._sourceLoaders.keys()):
			if srcName not in self._sourceClasses:
				if srcName not in self._sourceLoaders:
					self.log("WARNING: unknown source '%s'\n" % srcName)
					continue
				#if module not available
				srcModule = self._sourceLoaders[srcName].load_module('loki_source_%s' % srcName)
				srcClass = getattr(srcModule, 'Source_%s' % srcName)
				if not issubclass(srcClass, loki_source.Source):
					self.log("WARNING: invalid module for source '%s'\n" % srcName)
					continue
				self._sourceClasses[srcName] = srcClass
			#if module class not loaded
			srcSet.add(srcName)
		#foreach source
		return srcSet
	#loadSourceModules()
	
	
	def getSourceModuleVersions(self, sources=None):
		srcSet = self.loadSourceModules(sources)
		return { srcName : self._sourceClasses[srcName].getVersionString() for srcName in srcSet }
	#getSourceModuleVersions()
	
	
	def getSourceModuleOptions(self, sources=None):
		srcSet = self.loadSourceModules(sources)
		return { srcName : self._sourceClasses[srcName].getOptions() for srcName in srcSet }
	#getSourceModuleOptions()
	
	
	def attachSourceModules(self, sources=None):
		sources = self.loadSourceModules(sources)
		srcSet = set()
		for srcName in sources:
			if srcName not in self._sourceObjects:
				if srcName not in self._sourceClasses:
					raise Exception("loadSourceModules() reported false positive for '%s'" % srcName)
				self._sourceObjects[srcName] = self._sourceClasses[srcName](self._loki)
			#if module not instantiated
			srcSet.add(srcName)
		#foreach source
		return srcSet
	#attachSourceModules()
	
	
	def updateDatabase(self, sources=None, sourceOptions=None, cacheOnly=False):
		self._loki.testDatabaseUpdate()
		if self._updating:
			raise Exception('_updating set before updateDatabase()')
		
		# check for extraneous options
		srcSet = self.attachSourceModules(sources)
		srcOpts = sourceOptions or {}
		for srcName in srcOpts.keys():
			if srcName not in srcSet:
				self.log("WARNING: not updating from source '%s' for which options were supplied\n" % srcName)
		
		# update all specified sources
		iwd = os.path.abspath(os.getcwd())
		self._updating = True
		self._tablesUpdated = set()
		self._tablesDeindexed = set()
		with self._db:
			cursor = self._db.cursor()
			for srcName in sorted(srcSet):
				srcObj = self._sourceObjects[srcName]
				srcID = srcObj.getSourceID()
				
				# validate options, if any
				options = srcOpts.get(srcName, {})
				if options:
					self.logPush("validating %s options ...\n" % srcName)
					msg = srcObj.validateOptions(options)
					if msg != True:
						ok = False
						self.log("ERROR: %s\n" % msg)
						self.logPop()
						continue
					for opt,val in options.iteritems():
						self.log("%s = %s\n" % (opt,val))
					self.logPop("... OK\n")
				
				# store source options
				cursor.execute("DELETE FROM `db`.`source_option` WHERE source_id = ?", (srcID,))
				sql = "INSERT INTO `db`.`source_option` (source_id, option, value) VALUES (%d,?,?)" % srcID
				cursor.executemany(sql, options.iteritems())
				
				# switch to a temp subdirectory for this source
				path = os.path.join(iwd, srcName)
				if not os.path.exists(path):
					os.makedirs(path)
				os.chdir(path)
				
				# download files into a local cache
				if not cacheOnly:
					self.logPush("downloading %s data ...\n" % srcName)
					srcObj.download(options)
					self.logPop("... OK\n")
				
				# store source file metadata
				# all timestamps are assumed to be in UTC, but if a source
				# provides file timestamps with no TZ (like via FTP) we use them
				# as-is and assume they're supposed to be UTC
				self.log("analyzing %s data files ..." % srcName)
				cursor.execute("DELETE FROM `db`.`source_file` WHERE source_id = ?", (srcID,))
				sql = "INSERT INTO `db`.`source_file` (source_id, filename, size, modified, md5) VALUES (%d,?,?,DATETIME(?,'unixepoch'),?)" % srcID
				for filename in os.listdir('.'):
					stat = os.stat(filename)
					md5 = hashlib.md5()
					with open(filename,'rb') as f:
						chunk = f.read(8*1024*1024)
						while chunk:
							md5.update(chunk)
							chunk = f.read(8*1024*1024)
					cursor.execute(sql, (filename, stat.st_size, stat.st_mtime, md5.hexdigest()))
				self.log(" OK\n")
				
				# process new files
				self.logPush("processing %s data ...\n" % srcName)
				srcObj.update(options)
				self.logPop("... OK\n")
				
				# update source metadata
				sql = "UPDATE `db`.`source` SET updated = DATETIME('now'), version = ? WHERE source_id = ?"
				cursor.execute(sql, (srcObj.getVersionString(), srcID))
			#foreach source
			
			# cross-map GRCh/UCSChg build versions for all sources
			ucscGRC = collections.defaultdict(int)
			for row in self._db.cursor().execute("SELECT grch,ucschg FROM `db`.`grch_ucschg`"):
				ucscGRC[row[1]] = max(row[0], ucscGRC[row[1]])
				cursor.execute("UPDATE `db`.`source` SET grch = ? WHERE grch IS NULL AND ucschg = ?", (row[0],row[1]))
				cursor.execute("UPDATE `db`.`source` SET ucschg = ? WHERE ucschg IS NULL AND grch = ?", (row[1],row[0]))
			cursor.execute("UPDATE `db`.`source` SET current_ucschg = ucschg WHERE current_ucschg IS NULL")
			
			# check all sources' UCSChg build versions and set the latest as the target
			hgSources = collections.defaultdict(set)
			for row in cursor.execute("SELECT source_id, current_ucschg FROM `db`.`source` WHERE current_ucschg IS NOT NULL"):
				hgSources[row[1]].add(row[0])
			if hgSources:
				targetHG = max(hgSources)
				self.log("database genome build: GRCh%s / UCSChg%s\n" % (ucscGRC.get(targetHG,'?'), targetHG))
				targetUpdated = (int(self._loki.getDatabaseSetting('ucschg') or 0) != targetHG)
				self._loki.setDatabaseSetting('ucschg', targetHG)
			
			# liftOver sources with old build versions, if there are any
			if len(hgSources) > 1:
				locusSources = set(row[0] for row in cursor.execute("SELECT DISTINCT source_id FROM `db`.`snp_locus`"))
				regionSources = set(row[0] for row in cursor.execute("SELECT DISTINCT source_id FROM `db`.`biopolymer_region`"))
				chainsUpdated = ('grch_ucschg' in self._tablesUpdated or 'chain' in self._tablesUpdated or 'chain_data' in self._tablesUpdated)
				for oldHG in sorted(hgSources):
					if oldHG == targetHG:
						continue
					if not self._loki.hasLiftOverChains(oldHG, targetHG):
						self.log("ERROR: no chains available to lift hg%d to hg%d\n" % (oldHG, targetHG))
						continue
					
					if targetUpdated or chainsUpdated or 'snp_locus' in self._tablesUpdated:
						sourceIDs = hgSources[oldHG] & locusSources
						if sourceIDs:
							self.liftOverSNPLoci(oldHG, targetHG, sourceIDs)
					if targetUpdated or chainsUpdated or 'biopolymer_region' in self._tablesUpdated:
						sourceIDs = hgSources[oldHG] & regionSources
						if sourceIDs:
							self.liftOverRegions(oldHG, targetHG, sourceIDs)
					
					sql = "UPDATE `db`.`source` SET current_ucschg = %d WHERE source_id = ?" % targetHG
					cursor.executemany(sql, ((sourceID,) for sourceID in hgSources[oldHG]))
				#foreach old build
			#if any old builds
			
			# post-process as needed
			if 'snp_merge' in self._tablesUpdated or 'snp_locus' in self._tablesUpdated:
				self.updateMergedSNPLoci()
			if 'snp_merge' in self._tablesUpdated or 'snp_entrez_role' in self._tablesUpdated:
				self.updateMergedSNPEntrezRoles()
			if 'biopolymer_name' in self._tablesUpdated or 'biopolymer_name_name' in self._tablesUpdated:
				self.resolveBiopolymerNames()
			if 'biopolymer_name' in self._tablesUpdated or 'snp_entrez_role' in self._tablesUpdated:
				self.resolveSNPBiopolymerRoles()
			if 'biopolymer_name' in self._tablesUpdated or 'group_member_name' in self._tablesUpdated:
				self.resolveGroupMembers()
			if 'biopolymer_region' in self._tablesUpdated:
				self.updateBiopolymerZones()
			
			# reindex all remaining tables
			if self._tablesDeindexed:
				self.log("finalizing update ...")
				self._loki.createDatabaseIndecies(None, 'db', self._tablesDeindexed)
				self.log(" OK\n")
		#with db transaction
		self._updating = False
		self._tablesUpdated = None
		self._tablesDeindexed = None
		os.chdir(iwd)
	#updateDatabase()
	
	
	def liftOverSNPLoci(self, oldHG, newHG, sourceIDs):
		self.log("lifting over SNP loci from hg%d to hg%d ..." % (oldHG,newHG))
		self.prepareTableForUpdate('snp_locus')
		cursor = self._db.cursor()
		numLift = numNull = 0
		tally = dict()
		trash = set()
		
		# identify range of _ROWID_ in snp_locus
		# (two separate queries is faster because a simple MIN() or MAX() only peeks at the index;
		# SQLite isn't clever enough to do that for both at the same time, it does a table scan instead)
		firstRowID = min(row[0] for row in cursor.execute("SELECT MIN(_ROWID_) FROM `db`.`snp_locus`"))
		lastRowID = max(row[0] for row in cursor.execute("SELECT MAX(_ROWID_) FROM `db`.`snp_locus`"))
		
		# define a callback to store loci that can't be lifted over, for later deletion
		def errorCallback(region):
			trash.add( (region[0],) )
		
		# we can't SELECT and UPDATE the same table at the same time,
		# so read in batches of 2.5 million at a time based on _ROWID_
		minRowID = firstRowID
		maxRowID = minRowID + 2500000 - 1
		while minRowID <= lastRowID:
			sql = "SELECT _ROWID_, chr, pos FROM `db`.`snp_locus`"
			sql += " WHERE (_ROWID_ BETWEEN ? AND ?) AND source_id IN (%s)" % (','.join(str(i) for i in sourceIDs))
			oldLoci = list(cursor.execute(sql, (minRowID,maxRowID)))
			newLoci = self._loki.generateLiftOverLoci(oldHG, newHG, oldLoci, tally, errorCallback)
			sql = "UPDATE OR REPLACE `db`.`snp_locus` SET chr = ?2, pos = ?3 WHERE _ROWID_ = ?1"
			cursor.executemany(sql, newLoci)
			numLift += tally['lift']
			numNull += tally['null']
			if trash:
				cursor.executemany("DELETE FROM `db`.`snp_locus` WHERE _ROWID_ = ?", trash)
				trash.clear()
			minRowID = maxRowID + 1
			maxRowID = minRowID + 2500000 - 1
		#foreach batch
		
		self.log(" OK: %d loci lifted over, %d dropped\n" % (numLift,numNull))
	#liftOverSNPLoci()	
	
	
	def liftOverRegions(self, oldHG, newHG, sourceIDs):
		self.log("lifting over regions from hg%d to hg%d ..." % (oldHG,newHG))
		self.prepareTableForUpdate('biopolymer_region')
		cursor = self._db.cursor()
		numLift = numNull = 0
		tally = dict()
		trash = set()
		
		# identify range of _ROWID_ in biopolymer_region
		# (two separate queries is faster because a simple MIN() or MAX() only peeks at the index;
		# SQLite isn't clever enough to do that for both at the same time, it does a table scan instead)
		firstRowID = min(row[0] for row in cursor.execute("SELECT MIN(_ROWID_) FROM `db`.`biopolymer_region`"))
		lastRowID = max(row[0] for row in cursor.execute("SELECT MAX(_ROWID_) FROM `db`.`biopolymer_region`"))
		
		# define a callback to store regions that can't be lifted over, for later deletion
		def errorCallback(region):
			trash.add( (region[0],) )
		
		# we can't SELECT and UPDATE the same table at the same time,
		# so read in batches of 2.5 million at a time based on _ROWID_
		# (for regions this will probably be all of them in one go, but just in case)
		minRowID = firstRowID
		maxRowID = minRowID + 2500000 - 1
		while minRowID <= lastRowID:
			sql = "SELECT _ROWID_, chr, posMin, posMax FROM `db`.`biopolymer_region`"
			sql += " WHERE (_ROWID_ BETWEEN ? AND ?) AND source_id IN (%s)" % (','.join(str(i) for i in sourceIDs))
			oldRegions = list(cursor.execute(sql, (minRowID,maxRowID)))
			newRegions = self._loki.generateLiftOverRegions(oldHG, newHG, oldRegions, tally, errorCallback)
			sql = "UPDATE OR REPLACE `db`.`biopolymer_region` SET chr = ?2, posMin = ?3, posMax = ?4 WHERE _ROWID_ = ?1"
			cursor.executemany(sql, newRegions)
			numLift += tally['lift']
			numNull += tally['null']
			if trash:
				cursor.executemany("DELETE FROM `db`.`biopolymer_region` WHERE _ROWID_ = ?", trash)
				trash.clear()
			minRowID = maxRowID + 1
			maxRowID = minRowID + 2500000 - 1
		#foreach batch
		
		self.log(" OK: %d regions lifted over, %d dropped\n" % (numLift,numNull))
	#liftOverRegions()
	
	
	def updateMergedSNPLoci(self):
		self.log("checking for merged SNP loci ...")
		self._loki.testDatabaseUpdate()
		self.prepareTableForQuery('snp_locus')
		self.prepareTableForQuery('snp_merge')
		self._db.cursor().execute("""
INSERT OR IGNORE INTO `db`.`snp_locus` (rs, chr, pos, validated, source_id)
SELECT sm.rsCurrent, sl.chr, sl.pos, sl.validated, sl.source_id
FROM `db`.`snp_locus` AS sl
JOIN `db`.`snp_merge` AS sm
  ON sm.rsMerged = sl.rs
""")
		numCopied = self._db.changes()
		if numCopied:
			self.flagTableUpdate('snp_locus')
		self.log(" OK: %d loci copied\n" % (numCopied,))
	#updateMergedSNPLoci()
	
	
	def updateMergedSNPEntrezRoles(self):
		self.log("checking for merged SNP roles ...")
		self._loki.testDatabaseUpdate()
		self.prepareTableForQuery('snp_entrez_role')
		self.prepareTableForQuery('snp_merge')
		self._db.cursor().execute("""
INSERT OR IGNORE INTO `db`.`snp_entrez_role` (rs, entrez_id, role_id, source_id)
SELECT sm.rsCurrent, ser.entrez_id, ser.role_id, ser.source_id
FROM `db`.`snp_entrez_role` AS ser
JOIN `db`.`snp_merge` AS sm
  ON sm.rsMerged = ser.rs
""")
		numCopied = self._db.changes()
		if numCopied:
			self.flagTableUpdate('snp_entrez_role')
		self.log(" OK: %d roles copied\n" % (numCopied,))
	#updateMergedSNPEntrezRoles()
	
	
	def resolveBiopolymerNames(self):
		#TODO: iterative?
		self.log("resolving biopolymer names ...")
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		
		# calculate confidence scores for each possible name match
		dbc.execute("""
CREATE TEMP TABLE `temp`.`_biopolymer_name_name_score` (
  new_namespace_id INTERGER NOT NULL,
  new_name VARCHAR(256) NOT NULL,
  biopolymer_id INTEGER NOT NULL,
  polygenic TINYINT NOT NULL,
  implication INTEGER NOT NULL,
  PRIMARY KEY (new_namespace_id, new_name, biopolymer_id)
)
""")
		self.prepareTableForQuery('biopolymer_name_name')
		self.prepareTableForQuery('biopolymer_name')
		self.prepareTableForQuery('biopolymer')
		self.prepareTableForQuery('namespace')
		dbc.execute("""
INSERT INTO `temp`.`_biopolymer_name_name_score` (new_namespace_id, new_name, biopolymer_id, polygenic, implication)
/* calculate implication score for each possible match for each name */
SELECT
  bnn.new_namespace_id,
  bnn.new_name,
  bn.biopolymer_id,
  COALESCE(n.polygenic, 0) AS polygenic,
  COUNT(1) AS implication
FROM `db`.`biopolymer_name_name` AS bnn
JOIN `db`.`biopolymer_name` AS bn USING (name)
JOIN `db`.`biopolymer` AS b USING (biopolymer_id)
LEFT JOIN `db`.`namespace` AS n
  ON n.namespace_id = bnn.new_namespace_id
WHERE bnn.namespace_id IN (0, bn.namespace_id)
  AND bnn.type_id IN (0, b.type_id)
GROUP BY bnn.new_namespace_id, bnn.new_name, bn.biopolymer_id
""")
		
		# extrapolate new biopolymer_name records
		self.prepareTableForUpdate('biopolymer_name')
		dbc.execute("DELETE FROM `db`.`biopolymer_name` WHERE source_id = 0")
		dbc.execute("""
INSERT OR IGNORE INTO `db`.`biopolymer_name` (biopolymer_id, namespace_id, name, source_id)
/* identify specific match with the best score for each name */
SELECT
  biopolymer_id,
  new_namespace_id,
  new_name,
  0 AS source_id
FROM (
  /* identify names with only one best-score match */
  SELECT
    new_namespace_id,
    new_name,
    name_implication,
    SUM(CASE WHEN implication >= name_implication THEN 1 ELSE 0 END) AS match_implication
  FROM (
    /* identify best score for each name */
    SELECT
      new_namespace_id,
      new_name,
      MAX(implication) AS name_implication
    FROM `temp`.`_biopolymer_name_name_score`
    GROUP BY new_namespace_id, new_name
  )
  JOIN `temp`.`_biopolymer_name_name_score` USING (new_namespace_id, new_name)
  GROUP BY new_namespace_id, new_name
  HAVING polygenic > 0 OR match_implication = 1
)
JOIN `temp`.`_biopolymer_name_name_score` USING (new_namespace_id, new_name)
WHERE polygenic > 0 OR implication >= name_implication
""")
		
		# clean up
		dbc.execute("DROP TABLE `temp`.`_biopolymer_name_name_score`")
		numTotal = numUnrec = numMatch = 0
		self.prepareTableForQuery('biopolymer_name_name')
		self.prepareTableForQuery('biopolymer_name')
		self.prepareTableForQuery('biopolymer')
		numTotal = numUnrec = numMatch = 0
		for row in dbc.execute("""
SELECT COUNT(), SUM(CASE WHEN matches < 1 THEN 1 ELSE 0 END)
FROM (
  SELECT COUNT(DISTINCT b.biopolymer_id) AS matches
  FROM `db`.`biopolymer_name_name` AS bnn
  LEFT JOIN `db`.`biopolymer_name` AS bn
    ON bn.name = bnn.name
    AND bnn.namespace_id IN (0, bn.namespace_id)
  LEFT JOIN `db`.`biopolymer` AS b
    ON b.biopolymer_id = bn.biopolymer_id
    AND bnn.type_id IN (0, b.type_id)
  GROUP BY bnn.new_namespace_id, bnn.new_name
)
"""):
			numTotal = row[0] or 0
			numUnrec = row[1] or 0
		for row in dbc.execute("""
SELECT COUNT()
FROM (
  SELECT 1
  FROM `db`.`biopolymer_name`
  WHERE source_id = 0
  GROUP BY namespace_id, name
)
"""):
			numMatch = row[0] or 0
		numAmbig = numTotal - numUnrec - numMatch
		self.log(" OK: %d identifiers (%d ambiguous, %d unrecognized)\n" % (numMatch,numAmbig,numUnrec))
	#resolveBiopolymerNames()
	
	
	def resolveSNPBiopolymerRoles(self):
		self.log("resolving SNP roles ...")
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		
		typeID = self._loki.getTypeID('gene')
		namespaceID = self._loki.getNamespaceID('entrez_gid')
		if typeID and namespaceID:
			self.prepareTableForUpdate('snp_biopolymer_role')
			self.prepareTableForQuery('snp_entrez_role')
			self.prepareTableForQuery('biopolymer_name')
			dbc.execute("DELETE FROM `db`.`snp_biopolymer_role`")
			# we have to convert entrez_id to a string because the optimizer
			# won't use the index on biopolymer_name.name if the types don't match
			dbc.execute("""
INSERT OR IGNORE INTO `db`.`snp_biopolymer_role` (rs, biopolymer_id, role_id, source_id)
SELECT ser.rs, bn.biopolymer_id, ser.role_id, ser.source_id
FROM `db`.`snp_entrez_role` AS ser
JOIN `db`.`biopolymer_name` AS bn
  ON bn.namespace_id = ? AND bn.name = ''||ser.entrez_id
JOIN `db`.`biopolymer` AS b
  ON b.biopolymer_id = bn.biopolymer_id AND b.type_id = ?
""", (namespaceID,typeID))
		#if type[gene] and namespace[entrez_gid]
		
		#TODO: warning for unknown entrez_ids
		
		numTotal = numSNPs = numGenes = 0
		self.prepareTableForQuery('snp_biopolymer_role')
		for row in dbc.execute("SELECT COUNT(), COUNT(DISTINCT rs), COUNT(DISTINCT biopolymer_id) FROM `db`.`snp_biopolymer_role`"):
			numTotal = row[0]
			numSNPs = row[1]
			numGenes = row[2]
		self.log(" OK: %d roles (%d SNPs, %d genes)\n" % (numTotal,numSNPs,numGenes))
	#resolveSNPBiopolymerRoles()
	
	
	def resolveGroupMembers(self):
		self.log("resolving group members ...")
		self._loki.testDatabaseUpdate()
		dbc = self._db.cursor()
		
		# calculate confidence scores for each possible name match
		dbc.execute("""
CREATE TEMP TABLE `temp`.`_group_member_name_score` (
  group_id INTERGER NOT NULL,
  member INTEGER NOT NULL,
  biopolymer_id INTEGER NOT NULL,
  polynames INTEGER NOT NULL,
  implication INTEGER NOT NULL,
  quality INTEGER NOT NULL
)
""")
		self.prepareTableForQuery('group_member_name')
		self.prepareTableForQuery('biopolymer_name')
		self.prepareTableForQuery('biopolymer')
		self.prepareTableForQuery('namespace')
		dbc.execute("""
INSERT INTO `temp`.`_group_member_name_score` (group_id, member, biopolymer_id, polynames, implication, quality)
/* calculate implication and quality scores for each possible match for each member */
SELECT
  group_id,
  member,
  biopolymer_id,
  polynames,
  COUNT(DISTINCT gmn_rowid) AS implication,
  (CASE WHEN polynames > 0 THEN 1000 * COUNT(DISTINCT gmn_rowid) ELSE SUM(1000 / match_count) END) AS quality
FROM (
  /* count the number of possible matches for each name of each member */
  SELECT
    gmn._ROWID_ AS gmn_rowid,
    gmn.group_id,
    gmn.member,
    gmn.namespace_id,
    gmn.name,
    gmn.type_id,
    polynames,
    COUNT(DISTINCT bn.biopolymer_id) AS match_count
  FROM (
    /* count the number of matchable polyregion names for each member */
    SELECT
      gmn.group_id,
      gmn.member,
      COUNT(DISTINCT (CASE WHEN n.polygenic > 0 THEN gmn._ROWID_ ELSE NULL END)) AS polynames
    FROM `db`.`group_member_name` AS gmn
    JOIN `db`.`biopolymer_name` AS bn USING (name)
    JOIN `db`.`biopolymer` AS b USING (biopolymer_id)
    LEFT JOIN `db`.`namespace` AS n
      ON n.namespace_id = gmn.namespace_id
    WHERE gmn.namespace_id IN (0, bn.namespace_id)
      AND gmn.type_id IN (0, b.type_id)
    GROUP BY gmn.group_id, gmn.member
  )
  JOIN `db`.`group_member_name` AS gmn USING (group_id, member)
  JOIN `db`.`biopolymer_name` AS bn USING (name)
  JOIN `db`.`biopolymer` AS b USING (biopolymer_id)
  LEFT JOIN `db`.`namespace` AS n
    ON n.namespace_id = gmn.namespace_id
  WHERE gmn.namespace_id IN (0, bn.namespace_id)
    AND gmn.type_id IN (0, b.type_id)
    AND (n.polygenic > 0 OR polynames = 0)
  GROUP BY gmn.group_id, gmn.member, gmn.namespace_id, gmn.name
) AS gmn
JOIN `db`.`biopolymer_name` AS bn USING (name)
JOIN `db`.`biopolymer` AS b USING (biopolymer_id)
WHERE gmn.namespace_id IN (0, bn.namespace_id)
  AND gmn.type_id IN (0, b.type_id)
GROUP BY group_id, member, biopolymer_id
""")
		dbc.execute("CREATE INDEX `temp`.`_group_member_name_score__group_member_biopolymer` ON `_group_member_name_score` (group_id, member, biopolymer_id)")
		
		# generate group_biopolymer assignments with confidence scores
		self.prepareTableForUpdate('group_biopolymer')
		dbc.execute("DELETE FROM `db`.`group_biopolymer` WHERE source_id = 0")
		dbc.execute("""
/* group-biopolymer assignments with confidence scores */
INSERT INTO `db`.`group_biopolymer` (group_id, biopolymer_id, specificity, implication, quality, source_id)
SELECT
  group_id,
  biopolymer_id,
  MAX(specificity) AS specificity,
  MAX(implication) AS implication,
  MAX(quality) AS quality,
  0 AS source_id
FROM (
  /* identify specific matches with the best score for each member */
  SELECT
    group_id,
    member,
    biopolymer_id,
    (CASE
      WHEN polynames THEN 100 / member_variance
      ELSE 100 / match_basic
    END) AS specificity,
    (CASE
      WHEN polynames THEN 100 * implication / member_implication
      WHEN implication = member_implication THEN 100 / match_implication
      ELSE 0
    END) AS implication,
    (CASE
      WHEN polynames THEN 100 * quality / member_quality
      WHEN quality = member_quality THEN 100 / match_quality
      ELSE 0
    END) AS quality
  FROM (
    /* identify number of matches with the best score for each member */
    SELECT
      group_id,
      member,
      polynames,
      COUNT(DISTINCT implication) AS member_variance,
      member_implication,
      member_quality,
      COUNT() match_basic,
      SUM(CASE WHEN implication >= member_implication THEN 1 ELSE 0 END) AS match_implication,
      SUM(CASE WHEN quality >= member_quality THEN 1 ELSE 0 END) AS match_quality
    FROM (
      /* identify best scores for each member */
      SELECT
        group_id,
        member,
        polynames,
        MAX(implication) AS member_implication,
        MAX(quality) AS member_quality
      FROM `temp`.`_group_member_name_score`
      GROUP BY group_id, member, polynames
    )
    JOIN `temp`.`_group_member_name_score` USING (group_id, member, polynames)
    GROUP BY group_id, member, polynames
  )
  JOIN `temp`.`_group_member_name_score` USING (group_id, member, polynames)
  GROUP BY group_id, member, biopolymer_id
)
GROUP BY group_id, biopolymer_id
""")
		
		# generate group_biopolymer placeholders for unrecognized members
		self.prepareTableForUpdate('group_biopolymer')
		self.prepareTableForQuery('group_member_name')
		self.prepareTableForQuery('biopolymer_name')
		self.prepareTableForQuery('biopolymer')
		dbc.execute("""
INSERT INTO `db`.`group_biopolymer` (group_id, biopolymer_id, specificity, implication, quality, source_id)
SELECT
  group_id,
  0 AS biopolymer_id,
  COUNT() AS specificity,
  0 AS implication,
  0 AS quality,
  0 AS source_id
FROM (
  SELECT gmn.group_id
  FROM `db`.`group_member_name` AS gmn
  LEFT JOIN `db`.`biopolymer_name` AS bn
    ON bn.name = gmn.name
    AND gmn.namespace_id IN (0, bn.namespace_id)
  LEFT JOIN `db`.`biopolymer` AS b
    ON b.biopolymer_id = bn.biopolymer_id
    AND gmn.type_id IN (0, b.type_id)
  GROUP BY gmn.group_id, gmn.member
  HAVING MAX(b.biopolymer_id) IS NULL
)
GROUP BY group_id
""")
		
		# clean up
		dbc.execute("DROP TABLE `temp`.`_group_member_name_score`")
		numTotal = numSourced = numMatch = numAmbig = numUnrec = 0
		self.prepareTableForQuery('group_biopolymer')
		for row in dbc.execute("""
SELECT
  COALESCE(SUM(CASE WHEN biopolymer_id > 0 THEN 1 ELSE 0 END),0) AS total,
  COALESCE(SUM(CASE WHEN biopolymer_id > 0 AND source_id > 0 THEN 1 ELSE 0 END),0) AS sourced,
  COALESCE(SUM(CASE WHEN biopolymer_id > 0 AND source_id = 0 AND specificity >= 100 AND implication >= 100 AND quality >= 100 THEN 1 ELSE 0 END),0) AS definite,
  COALESCE(SUM(CASE WHEN biopolymer_id > 0 AND source_id = 0 AND (specificity < 100 OR implication < 100 OR quality < 100) THEN 1 ELSE 0 END),0) AS conditional,
  COALESCE(SUM(CASE WHEN biopolymer_id = 0 AND source_id = 0 THEN specificity ELSE 0 END),0) AS unmatched
FROM `db`.`group_biopolymer`
"""):
			numTotal = row[0]
			numSourced = row[1]
			numMatch = row[2]
			numAmbig = row[3]
			numUnrec = row[4]
		self.log(" OK: %d associations (%d explicit, %d definite, %d conditional, %d unrecognized)\n" % (numTotal,numSourced,numMatch,numAmbig,numUnrec))
	#resolveGroupMembers()
	
	
	def updateBiopolymerZones(self):
		self.log("calculating zone coverage ...")
		self._loki.testDatabaseUpdate()
		size = self._loki.getDatabaseSetting('zone_size')
		if not size:
			raise Exception("ERROR: could not determine database setting 'zone_size'")
		size = int(size)
		dbc = self._db.cursor()
		
		# make sure all regions are correctly oriented
		dbc.execute("UPDATE `db`.`biopolymer_region` SET posMin = posMax, posMax = posMin WHERE posMin > posMax")
		
		# define zone generator
		def _zones(size, regions):
			# regions=[ (id,chr,posMin,posMax),... ]
			# yields:[ (id,chr,zone),... ]
			for r in regions:
				for z in xrange(int(r[2]/size),int(r[3]/size)+1):
					yield (r[0],r[1],z)
		#_zones()
		
		# feed all regions through the zone generator
		self.prepareTableForUpdate('biopolymer_zone')
		self.prepareTableForQuery('biopolymer_region')
		dbc.execute("DELETE FROM `db`.`biopolymer_zone`")
		dbc.executemany(
			"INSERT OR IGNORE INTO `db`.`biopolymer_zone` (biopolymer_id,chr,zone) VALUES (?,?,?)",
			_zones(
				size,
				self._db.cursor().execute("SELECT biopolymer_id,chr,MIN(posMin),MAX(posMax) FROM `db`.`biopolymer_region` GROUP BY biopolymer_id, chr")
			)
		)
		
		# clean up
		self.prepareTableForQuery('biopolymer_zone')
		for row in dbc.execute("SELECT COUNT(), COUNT(DISTINCT biopolymer_id) FROM `db`.`biopolymer_zone`"):
			numTotal = row[0]
			numGenes = row[1]
		self.log(" OK: %d records (%d regions)\n" % (numTotal,numGenes))
	#updateBiopolymerZones()
	
	
#Updater
