#!/usr/bin/env python

import collections
import hashlib
import os
import pkgutil
import sys
import traceback

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
				#print "deindexing %s" % table #DEBUG
				self._tablesDeindexed.add(table)
				self._loki.dropDatabaseIndecies(None, 'db', table)
	#prepareTableForUpdate()
	
	
	def prepareTableForQuery(self, table):
		if self._updating:
			if table in self._tablesDeindexed:
				#print "reindexing %s" % table DEBUG
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
	
	
	def updateDatabase(self, sources=None, sourceOptions=None, cacheOnly=False, forceUpdate=False):
		if self._updating:
			raise Exception("_updating set before updateDatabase()")
		self._loki.testDatabaseWriteable()
		if self._loki.getDatabaseSetting('finalized',int):
			raise Exception("cannot update a finalized database")
		
		# check for extraneous options
		self.logPush("preparing for update ...\n")
		srcSet = self.attachSourceModules(sources)
		srcOpts = sourceOptions or {}
		for srcName in srcOpts.keys():
			if srcName not in srcSet:
				self.log("WARNING: not updating from source '%s' for which options were supplied\n" % srcName)
		logIndent = self.logPop("... OK\n")
		
		# update all specified sources
		iwd = os.path.abspath(os.getcwd())
		self._updating = True
		self._tablesUpdated = set()
		self._tablesDeindexed = set()
		srcErrors = set()
		cursor = self._db.cursor()
		cursor.execute("SAVEPOINT 'updateDatabase'")
		try:
			for srcName in sorted(srcSet):
				cursor.execute("SAVEPOINT 'updateDatabase_%s'" % (srcName,))
				try:
					srcObj = self._sourceObjects[srcName]
					srcID = srcObj.getSourceID()
					
					# validate options, if any
					options = srcOpts.get(srcName, {})
					if options:
						self.logPush("validating %s options ...\n" % srcName)
						msg = srcObj.validateOptions(options)
						if msg != True:
							raise Exception(msg)
						for opt,val in options.iteritems():
							self.log("%s = %s\n" % (opt,val))
						self.logPop("... OK\n")
					
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
					
					# calculate source file metadata
					# all timestamps are assumed to be in UTC, but if a source
					# provides file timestamps with no TZ (like via FTP) we use them
					# as-is and assume they're supposed to be UTC
					self.log("analyzing %s data files ..." % srcName)
					filehash = dict()
					for filename in os.listdir('.'):
						stat = os.stat(filename)
						md5 = hashlib.md5()
						with open(filename,'rb') as f:
							chunk = f.read(8*1024*1024)
							while chunk:
								md5.update(chunk)
								chunk = f.read(8*1024*1024)
						filehash[filename] = (filename, long(stat.st_size), long(stat.st_mtime), md5.hexdigest())
					self.log(" OK\n")
					
					# compare current loader version, options and file metadata to the last update
					skip = not forceUpdate
					last = '?'
					if skip:
						for row in cursor.execute("SELECT version, DATETIME(updated,'localtime') FROM `db`.`source` WHERE source_id = ?", (srcID,)):
							skip = skip and (row[0] == srcObj.getVersionString())
							last = row[1]
					if skip:
						n = 0
						for row in cursor.execute("SELECT option, value FROM `db`.`source_option` WHERE source_id = ?", (srcID,)):
							n += 1
							skip = skip and (row[0] in options) and (row[1] == options[row[0]])
						skip = skip and (n == len(options))
					if skip:
						n = 0
						for row in cursor.execute("SELECT filename, size, md5 FROM `db`.`source_file` WHERE source_id = ?", (srcID,)):
							n += 1
							skip = skip and (row[0] in filehash) and (row[1] == filehash[row[0]][1]) and (row[2] == filehash[row[0]][3])
						skip = skip and (n == len(filehash))
					
					# skip the update if the current loader and all source file versions match the last update
					if skip:
						self.log("skipping %s update, no data or software changes since %s\n" % (srcName,last))
					else:
						# process new files (or old files with a new loader)
						self.logPush("processing %s data ...\n" % srcName)
						
						cursor.execute("DELETE FROM `db`.`warning` WHERE source_id = ?", (srcID,))
						srcObj.update(options)
						cursor.execute("UPDATE `db`.`source` SET updated = DATETIME('now'), version = ? WHERE source_id = ?", (srcObj.getVersionString(), srcID))
						
						cursor.execute("DELETE FROM `db`.`source_option` WHERE source_id = ?", (srcID,))
						sql = "INSERT INTO `db`.`source_option` (source_id, option, value) VALUES (%d,?,?)" % srcID
						cursor.executemany(sql, options.iteritems())
						
						cursor.execute("DELETE FROM `db`.`source_file` WHERE source_id = ?", (srcID,))
						sql = "INSERT INTO `db`.`source_file` (source_id, filename, size, modified, md5) VALUES (%d,?,?,DATETIME(?,'unixepoch'),?)" % srcID
						cursor.executemany(sql, filehash.values())
						
						self.logPop("... OK\n")
					#if skip
				except:
					srcErrors.add(srcName)
					excType,excVal,excTrace = sys.exc_info()
					while self.logPop() > logIndent:
						pass
					self.logPush("ERROR: failed to update %s\n" % (srcName,))
					if excTrace:
						for line in traceback.format_list(traceback.extract_tb(excTrace)[-1:]):
							self.log(line)
					for line in traceback.format_exception_only(excType,excVal):
						self.log(line)
					self.logPop()
					cursor.execute("ROLLBACK TRANSACTION TO SAVEPOINT 'updateDatabase_%s'" % (srcName,))
				finally:
					cursor.execute("RELEASE SAVEPOINT 'updateDatabase_%s'" % (srcName,))
				#try/except/finally
			#foreach source
			
			# pull the latest GRCh/UCSChg conversions
			#   http://genome.ucsc.edu/FAQ/FAQreleases.html
			#   http://genome.ucsc.edu/goldenPath/releaseLog.html
			# TODO: find a better machine-readable source for this data
			if not cacheOnly:
				self.log("updating GRCh:UCSChg genome build identities ...")
				import urllib2
				import re
				response = urllib2.urlopen('http://genome.ucsc.edu/FAQ/FAQreleases.html')
				page = ""
				while True:
					data = response.read()
					if not data:
						break
					page += data
				rowHuman = False
				for tablerow in re.finditer(r'<tr>.*?</tr>', page, re.IGNORECASE | re.DOTALL):
					cols = tuple(match.group()[4:-5].strip().lower() for match in re.finditer(r'<td>.*?</td>', tablerow.group(), re.IGNORECASE | re.DOTALL))
					if cols and ((cols[0] == 'human') or (rowHuman and (cols[0] in ('','&nbsp;')))):
						rowHuman = True
						grch = ucschg = None
						try:
							if cols[1].startswith('hg'):
								ucschg = int(cols[1][2:])
							if cols[3].startswith('genome reference consortium grch'):
								grch = int(cols[3][32:])
							if cols[3].startswith('ncbi build '):
								grch = int(cols[3][11:])
						except:
							pass
						if grch and ucschg:
							cursor.execute("INSERT OR REPLACE INTO `db`.`grch_ucschg` (grch,ucschg) VALUES (?,?)", (grch,ucschg))
					else:
						rowHuman = False
				#foreach tablerow
				self.log(" OK\n")
			#if not cacheOnly
			
			# cross-map GRCh/UCSChg build versions for all sources
			ucscGRC = collections.defaultdict(int)
			for row in self._db.cursor().execute("SELECT grch,ucschg FROM `db`.`grch_ucschg`"):
				ucscGRC[row[1]] = max(row[0], ucscGRC[row[1]])
				cursor.execute("UPDATE `db`.`source` SET grch = ? WHERE grch IS NULL AND ucschg = ?", (row[0],row[1]))
				cursor.execute("UPDATE `db`.`source` SET ucschg = ? WHERE ucschg IS NULL AND grch = ?", (row[1],row[0]))
			cursor.execute("UPDATE `db`.`source` SET current_ucschg = ucschg WHERE current_ucschg IS NULL")
			
			# check for any source with an unrecognized GRCh build
			mismatch = False
			for row in cursor.execute("SELECT source, grch, ucschg FROM `db`.`source` WHERE (grch IS NULL) != (ucschg IS NULL)"):
				self.log("WARNING: unrecognized genome build for '%s' (NCBI GRCh%s, UCSC hg%s)\n" % (row[0],(row[1] or "?"),(row[2] or "?")))
				mismatch = True
			if mismatch:
				self.log("WARNING: database may contain incomparable genome positions!\n")
			
			# check all sources' UCSChg build versions and set the latest as the target
			hgSources = collections.defaultdict(set)
			for row in cursor.execute("SELECT source_id, current_ucschg FROM `db`.`source` WHERE current_ucschg IS NOT NULL"):
				hgSources[row[1]].add(row[0])
			if hgSources:
				targetHG = max(hgSources)
				self.log("database genome build: GRCh%s / UCSChg%s\n" % (ucscGRC.get(targetHG,'?'), targetHG))
				targetUpdated = (self._loki.getDatabaseSetting('ucschg',int) != targetHG)
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
			#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'snp_merge' in self._tablesUpdated:
				self.cleanupSNPMerges()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'snp_merge' in self._tablesUpdated or 'snp_locus' in self._tablesUpdated:
				self.updateMergedSNPLoci()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'snp_locus' in self._tablesUpdated:
				self.cleanupSNPLoci()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'snp_merge' in self._tablesUpdated or 'snp_entrez_role' in self._tablesUpdated:
				self.updateMergedSNPEntrezRoles()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'snp_entrez_role' in self._tablesUpdated:
				self.cleanupSNPEntrezRoles()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'snp_merge' in self._tablesUpdated or 'gwas' in self._tablesUpdated:
				self.updateMergedGWASAnnotations()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'biopolymer_name' in self._tablesUpdated or 'biopolymer_name_name' in self._tablesUpdated:
				self.resolveBiopolymerNames()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'biopolymer_name' in self._tablesUpdated or 'snp_entrez_role' in self._tablesUpdated:
				self.resolveSNPBiopolymerRoles()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'biopolymer_name' in self._tablesUpdated or 'group_member_name' in self._tablesUpdated:
				self.resolveGroupMembers()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			if 'biopolymer_region' in self._tablesUpdated:
				self.updateBiopolymerZones()
				#self.log("MEMORY: %d bytes (%d peak)\n" % self._loki.getDatabaseMemoryUsage()) #DEBUG
			
			# reindex all remaining tables
			self.log("finishing update ...")
			if self._tablesDeindexed:
				self._loki.createDatabaseIndecies(None, 'db', self._tablesDeindexed)
			if self._tablesUpdated:
				self._loki.setDatabaseSetting('optimized',0)
			self.log(" OK\n")
		except:
			excType,excVal,excTrace = sys.exc_info()
			while self.logPop() > logIndent:
				pass
			self.logPush("ERROR: failed to update the database\n")
			if excTrace:
				for line in traceback.format_list(traceback.extract_tb(excTrace)[-1:]):
					self.log(line)
			for line in traceback.format_exception_only(excType,excVal):
				self.log(line)
			self.logPop()
			cursor.execute("ROLLBACK TRANSACTION TO SAVEPOINT 'updateDatabase'")
		finally:
			cursor.execute("RELEASE SAVEPOINT 'updateDatabase'")
			self._updating = False
			self._tablesUpdated = None
			self._tablesDeindexed = None
			os.chdir(iwd)
		#try/except/finally
		
		# report and return
		if srcErrors:
			self.logPush("WARNING: data from these sources was not updated:\n")
			for srcName in sorted(srcErrors):
				self.log("%s\n" % srcName)
			self.logPop()
			return False
		return True
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
			sql = "SELECT _ROWID_, chr, pos, NULL FROM `db`.`snp_locus`"
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
			sql = "SELECT _ROWID_, chr, posMin, posMax, NULL FROM `db`.`biopolymer_region`"
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
	
	
	def cleanupSNPMerges(self):
		self.log("verifying SNP merge records ...")
		self.prepareTableForQuery('snp_merge')
		dbc = self._db.cursor()
		
		# for each set of ROWIDs which constitute a duplicated snp merge, cull all but one
		cull = set()
		sql = "SELECT GROUP_CONCAT(_ROWID_) FROM `db`.`snp_merge` GROUP BY rsMerged HAVING COUNT() > 1"
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		for row in dbc.execute(sql):
			cull.update( (long(i),) for i in row[0].split(',')[1:] )
		#last = None
		#for row in dbc.execute("SELECT _ROWID_, rsMerged FROM `db`.`snp_merge` ORDER BY rsMerged"):
		#	if last == row[1]:
		#		cull.add(row[0:1])
		#	last = row[1]
		if cull:
			self.flagTableUpdate('snp_merge')
			dbc.executemany("DELETE FROM `db`.`snp_merge` WHERE _ROWID_ = ?", cull)
		self.log(" OK: %d duplicate merges\n" % (len(cull),))
	#cleanupSNPMerges()
	
	
	def updateMergedSNPLoci(self):
		self.log("checking for merged SNP loci ...")
		self.prepareTableForQuery('snp_locus')
		self.prepareTableForQuery('snp_merge')
		dbc = self._db.cursor()
		sql = """
INSERT INTO `db`.`snp_locus` (rs, chr, pos, validated, source_id)
SELECT sm.rsCurrent, sl.chr, sl.pos, sl.validated, sl.source_id
FROM `db`.`snp_locus` AS sl
JOIN `db`.`snp_merge` AS sm
  ON sm.rsMerged = sl.rs
"""
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		dbc.execute(sql)
		numCopied = self._db.changes()
		if numCopied:
			self.flagTableUpdate('snp_locus')
		self.log(" OK: %d loci copied\n" % (numCopied,))
	#updateMergedSNPLoci()
	
	
	def cleanupSNPLoci(self):
		self.log("verifying SNP loci ...")
		self.prepareTableForQuery('snp_locus')
		dbc = self._db.cursor()
		# for each set of ROWIDs which constitute a duplicated snp-locus, cull all but one
		# but, make sure that if any of the originals were validated, the remaining one is also
		valid = set()
		cull = set()
		sql = "SELECT GROUP_CONCAT(_ROWID_), MAX(validated) FROM `db`.`snp_locus` GROUP BY rs, chr, pos HAVING COUNT() > 1"
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		for row in dbc.execute(sql):
			rowids = row[0].split(',')
			if row[1]:
				valid.add( (long(rowids[0]),) )
			cull.update( (long(i),) for i in rowids[1:] )
		#last = None
		#for row in dbc.execute("SELECT _ROWID_, rs||':'||chr||':'||pos, validated FROM `db`.`snp_locus` ORDER BY rs, chr, pos"):
		#	if last == row[1]:
		#		cull.add(row[0:1])
		#		if row[2]:
		#			valid.add( last.split(':') )
		#	last = row[1]
		if valid:
			dbc.executemany("UPDATE `db`.`snp_locus` SET validated = 1 WHERE _ROWID_ = ?", valid)
			#dbc.executemany("UPDATE `db`.`snp_locus` SET validated = 1 WHERE rs = ? AND chr = ? AND pos = ?", valid)
		if cull:
			self.flagTableUpdate('snp_locus')
			dbc.executemany("DELETE FROM `db`.`snp_locus` WHERE _ROWID_ = ?", cull)
		self.log(" OK: %d duplicate loci\n" % (len(cull),))
	#cleanupSNPLoci()
	
	
	def updateMergedSNPEntrezRoles(self):
		self.log("checking for merged SNP roles ...")
		self.prepareTableForQuery('snp_entrez_role')
		self.prepareTableForQuery('snp_merge')
		dbc = self._db.cursor()
		sql = """
INSERT OR IGNORE INTO `db`.`snp_entrez_role` (rs, entrez_id, role_id, source_id)
SELECT sm.rsCurrent, ser.entrez_id, ser.role_id, ser.source_id
FROM `db`.`snp_entrez_role` AS ser
JOIN `db`.`snp_merge` AS sm
  ON sm.rsMerged = ser.rs
"""
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		dbc.execute(sql)
		numCopied = self._db.changes()
		if numCopied:
			self.flagTableUpdate('snp_entrez_role')
		self.log(" OK: %d roles copied\n" % (numCopied,))
	#updateMergedSNPEntrezRoles()
	
	
	def cleanupSNPEntrezRoles(self):
		self.log("verifying SNP roles ...")
		self.prepareTableForQuery('snp_entrez_role')
		dbc = self._db.cursor()
		cull = set()
		sql = "SELECT GROUP_CONCAT(_ROWID_) FROM `db`.`snp_entrez_role` GROUP BY rs, entrez_id, role_id HAVING COUNT() > 1"
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		for row in dbc.execute(sql):
			cull.update( (long(i),) for i in row[0].split(',')[1:] )
		#last = None
		#for row in dbc.execute("SELECT _ROWID_, rs||':'||entrez_id||':'||role_id FROM `db`.`snp_entrez_role` ORDER BY rs, entrez_id, role_id"):
		#	if last == row[1]:
		#		cull.add(row[0:1])
		#	last = row[1]
		if cull:
			self.flagTableUpdate('snp_entrez_role')
			dbc.executemany("DELETE FROM `db`.`snp_entrez_role` WHERE _ROWID_ = ?", cull)
		self.log(" OK: %d duplicate roles\n" % (len(cull),))
	#cleanupSNPEntrezRoles()
	
	
	def updateMergedGWASAnnotations(self):
		self.log("checking for merged GWAS annotated SNPs ...")
		self.prepareTableForQuery('gwas')
		self.prepareTableForQuery('snp_merge')
		dbc = self._db.cursor()
		sql = """
INSERT INTO `db`.`gwas` (rs, chr, pos, trait, snps, orbeta, allele95ci, riskAfreq, pubmed_id, source_id)
SELECT sm.rsCurrent, w.chr, w.pos, w.trait, w.snps, w.orbeta, w.allele95ci, w.riskAfreq, w.pubmed_id, w.source_id
FROM `db`.`gwas` AS w
JOIN `db`.`snp_merge` AS sm
  ON sm.rsMerged = w.rs
"""
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		dbc.execute(sql)
		numCopied = self._db.changes()
		if numCopied:
			self.flagTableUpdate('gwas')
		self.log(" OK: %d annotations copied\n" % (numCopied,))
	#updateMergedGWASAnnotations()
	
	
	def resolveBiopolymerNames(self):
		#TODO: iterative?
		self.log("resolving biopolymer names ...")
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
		dbc = self._db.cursor()
		
		typeID = self._loki.getTypeID('gene')
		namespaceID = self._loki.getNamespaceID('entrez_gid')
		numUnrec = 0
		if typeID and namespaceID:
			self.prepareTableForUpdate('snp_biopolymer_role')
			self.prepareTableForQuery('snp_entrez_role')
			self.prepareTableForQuery('biopolymer_name')
			dbc.execute("DELETE FROM `db`.`snp_biopolymer_role`")
			# we have to convert entrez_id to a string because the optimizer
			# won't use the index on biopolymer_name.name if the types don't match
			dbc.execute("""
INSERT INTO `db`.`snp_biopolymer_role` (rs, biopolymer_id, role_id, source_id)
SELECT ser.rs, bn.biopolymer_id, ser.role_id, ser.source_id
FROM `db`.`snp_entrez_role` AS ser
JOIN `db`.`biopolymer_name` AS bn
  ON bn.namespace_id = ? AND bn.name = ''||ser.entrez_id
JOIN `db`.`biopolymer` AS b
  ON b.biopolymer_id = bn.biopolymer_id AND b.type_id = ?
""", (namespaceID,typeID))
			numUnrec = sum(row[0] for row in dbc.execute("""
SELECT COUNT() FROM (
SELECT 1
FROM `db`.`snp_entrez_role` AS ser
LEFT JOIN `db`.`biopolymer_name` AS bn
  ON bn.namespace_id = ? AND bn.name = ''||ser.entrez_id
LEFT JOIN `db`.`biopolymer` AS b
  ON b.biopolymer_id = bn.biopolymer_id AND b.type_id = ?
GROUP BY ser._ROWID_
HAVING MAX(b.biopolymer_id) IS NULL
)
""", (namespaceID,typeID)))
		#if type[gene] and namespace[entrez_gid]
		
		self.prepareTableForQuery('snp_biopolymer_role')
		cull = set()
		sql = "SELECT GROUP_CONCAT(_ROWID_) FROM `db`.`snp_biopolymer_role` GROUP BY rs, biopolymer_id, role_id HAVING COUNT() > 1"
		#for row in dbc.execute("EXPLAIN QUERY PLAN "+sql): #DEBUG
		#	print row
		for row in dbc.execute(sql):
			cull.update( (long(i),) for i in row[0].split(',')[1:] )
		#last = None
		#for row in dbc.execute("SELECT _ROWID_, rs||':'||biopolymer_id||':'||role_id FROM `db`.`snp_biopolymer_role` ORDER BY rs, biopolymer_id, role_id"):
		#	if last == row[1]:
		#		cull.add(row[0:1])
		#	last = row[1]
		if cull:
			self.flagTableUpdate('snp_biopolymer_role')
			dbc.executemany("DELETE FROM `db`.`snp_biopolymer_role` WHERE _ROWID_ = ?", cull)
		
		numTotal = numSNPs = numGenes = 0
		for row in dbc.execute("SELECT COUNT(), COUNT(DISTINCT rs), COUNT(DISTINCT biopolymer_id) FROM `db`.`snp_biopolymer_role`"):
			numTotal = row[0]
			numSNPs = row[1]
			numGenes = row[2]
		self.log(" OK: %d roles (%d SNPs, %d genes; %d unrecognized)\n" % (numTotal,numSNPs,numGenes,numUnrec))
	#resolveSNPBiopolymerRoles()
	
	
	def resolveGroupMembers(self):
		self.log("resolving group members ...")
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
		size = self._loki.getDatabaseSetting('zone_size',int)
		if not size:
			raise Exception("ERROR: could not determine database setting 'zone_size'")
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
