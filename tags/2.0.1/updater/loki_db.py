#!/usr/bin/env python

import apsw
import bisect
import itertools
import sys


class Database(object):
	
	
	##################################################
	# public class data
	
	
	ver_maj,ver_min,ver_rev,ver_dev,ver_date = 2,0,0,'a7','2012-08-30'
	
	# hardcode translations between chromosome numbers and textual tags
	chr_num = {}
	chr_name = {}
	cnum = 0
	for cname in ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','XY','MT'):
		cnum += 1
		chr_num[cnum] = cnum
		chr_num['%s' % cnum] = cnum
		chr_num[cname] = cnum
		chr_name[cnum] = cname
		chr_name['%s' % cnum] = cname
		chr_name[cname] = cname
	chr_num['M'] = chr_num['MT']
	chr_name['M'] = chr_name['MT']
	
	
	##################################################
	# private class data
	
	
	_schema = {
		'db': {
			##################################################
			# configuration tables
			
			
			'setting': {
				'table': """
(
  setting VARCHAR(32) PRIMARY KEY NOT NULL,
  value VARCHAR(256)
)
""",
				'data': [
					('zone_size','100000'),
					('finalized','0'),
				],
				'index': {}
			}, #.db.setting
			
			
			##################################################
			# metadata tables
			
			
			'ldprofile': {
				'table': """
(
  ldprofile_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ldprofile VARCHAR(32) UNIQUE NOT NULL,
  comment VARCHAR(128),
  description VARCHAR(128)
)
""",
				'index': {}
			}, #.db.ldprofile
			
			
			'namespace': {
				'table': """
(
  namespace_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  namespace VARCHAR(32) UNIQUE NOT NULL,
  polygenic TINYINT NOT NULL DEFAULT 0
)
""",
				'index': {}
			}, #.db.namespace
			
			
			'relationship': {
				'table': """
(
  relationship_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  relationship VARCHAR(32) UNIQUE NOT NULL
)
""",
				'index': {}
			}, #.db.relationship
			
			
			'role': {
				'table': """
(
  role_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  role VARCHAR(32) UNIQUE NOT NULL,
  description VARCHAR(128),
  coding TINYINT,
  exon TINYINT
)
""",
				'index': {}
			}, #.db.role
			
			
			'source': {
				'table': """
(
  source_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  source VARCHAR(32) UNIQUE NOT NULL,
  updated DATETIME,
  version VARCHAR(32),
  grch INTEGER,
  ucschg INTEGER,
  current_ucschg INTEGER
)
""",
				'index': {}
			}, #.db.source
			
			
			'source_file': {
				'table': """
(
  source_id TINYINT NOT NULL,
  filename VARCHAR(256) NOT NULL,
  size BIGINT,
  modified DATETIME,
  md5 VARCHAR(64),
  PRIMARY KEY (source_id, filename)
)
""",
				'index': {}
			}, #.db.source_file
			
			
			'type': {
				'table': """
(
  type_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  type VARCHAR(32) UNIQUE NOT NULL
)
""",
				'index': {}
			}, #.db.type
			
			
			'warning': {
				'table': """
(
  warning_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  source_id TINYINT NOT NULL,
  warning VARCHAR(8192)
)
""",
				'index': {
					'warning__source': '(source_id)',
				}
			}, #.db.warning
			
			
			##################################################
			# snp tables
			
			
			'snp_merge': {
				'table': """
(
  rsMerged INTEGER PRIMARY KEY NOT NULL,
  rsCurrent INTEGER NOT NULL,
  source_id TINYINT NOT NULL
)
""",
				'index': {}
			}, #.db.snp_merge
			
			
			'snp_locus': {
				'table': """
(
  rs INTEGER NOT NULL,
  chr TINYINT NOT NULL,
  pos BIGINT NOT NULL,
  validated TINYINT NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (rs,chr,pos)
)
""",
				'index': {
					'snp_locus__chr_pos_rs': '(chr,pos,rs)',
					# a (validated,...) index would be nice but adds >1GB to the file size :/
					#'snp_locus__valid_chr_pos_rs': '(validated,chr,pos,rs)',
				}
			}, #.db.snp_locus
			
			
			'snp_entrez_role': {
				'table': """
(
  rs INTEGER NOT NULL,
  entrez_id INTEGER NOT NULL,
  role_id INTEGER NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (entrez_id,rs,role_id)
)
""",
				'index': {}
			}, #.db.snp_entrez_role
			
			
			'snp_biopolymer_role': {
				'table': """
(
  rs INTEGER NOT NULL,
  biopolymer_id INTEGER NOT NULL,
  role_id INTEGER NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (rs,biopolymer_id,role_id)
)
""",
				'index': {
					'snp_biopolymer_role__biopolymer_rs_role': '(biopolymer_id,rs,role_id)',
				}
			}, #.db.snp_biopolymer_role
			
			
			##################################################
			# biopolymer tables
			
			
			'biopolymer': {
				'table': """
(
  biopolymer_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  type_id TINYINT NOT NULL,
  label VARCHAR(64) NOT NULL,
  description VARCHAR(256),
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'biopolymer__type': '(type_id)',
					'biopolymer__label_type': '(label,type_id)',
				}
			}, #.db.biopolymer
			
			
			'biopolymer_name': {
				'table': """
(
  biopolymer_id INTEGER NOT NULL,
  namespace_id INTEGER NOT NULL,
  name VARCHAR(256) NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (biopolymer_id,namespace_id,name)
)
""",
				'index': {
					'biopolymer_name__name_namespace_biopolymer': '(name,namespace_id,biopolymer_id)',
				}
			}, #.db.biopolymer_name
			
			
			'biopolymer_name_name': {
				# PRIMARY KEY column order satisfies the need to GROUP BY new_namespace_id, new_name
				'table': """
(
  namespace_id INTEGER NOT NULL,
  name VARCHAR(256) NOT NULL,
  type_id TINYINT NOT NULL,
  new_namespace_id INTEGER NOT NULL,
  new_name VARCHAR(256) NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (new_namespace_id,new_name,type_id,namespace_id,name)
)
""",
				'index': {}
			}, #.db.biopolymer_name_name
			
			
			'biopolymer_region': {
				'table': """
(
  biopolymer_id INTEGER NOT NULL,
  ldprofile_id INTEGER NOT NULL,
  chr TINYINT NOT NULL,
  posMin BIGINT NOT NULL,
  posMax BIGINT NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (biopolymer_id,ldprofile_id,chr,posMin,posMax)
)
""",
				'index': {
					'biopolymer_region__ldprofile_chr_min': '(ldprofile_id,chr,posMin)',
					'biopolymer_region__ldprofile_chr_max': '(ldprofile_id,chr,posMax)',
				}
			}, #.db.biopolymer_region
			
			
			'biopolymer_zone': {
				'table': """
(
  biopolymer_id INTEGER NOT NULL,
  chr TINYINT NOT NULL,
  zone INTEGER NOT NULL,
  PRIMARY KEY (biopolymer_id,chr,zone)
)
""",
				'index': {
					'biopolymer_zone__zone': '(chr,zone,biopolymer_id)',
				}
			}, #.db.biopolymer_zone
			
			
			##################################################
			# group tables
			
			
			'group': {
				'table': """
(
  group_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  type_id TINYINT NOT NULL,
  label VARCHAR(64) NOT NULL,
  description VARCHAR(256),
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'group__type': '(type_id)',
					'group__label_type': '(label,type_id)',
				}
			}, #.db.group
			
			
			'group_name': {
				'table': """
(
  group_id INTEGER NOT NULL,
  namespace_id INTEGER NOT NULL,
  name VARCHAR(256) NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (group_id,namespace_id,name)
)
""",
				'index': {
					'group_name__name_namespace_group': '(name,namespace_id,group_id)',
					'group_name__source_name': '(source_id,name)',
				}
			}, #.db.group_name
			
			
			'group_group': {
				'table': """
(
  group_id INTEGER NOT NULL,
  related_group_id INTEGER NOT NULL,
  relationship_id SMALLINT NOT NULL,
  direction TINYINT NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (group_id,related_group_id,relationship_id,direction)
)
""",
				'index': {
					'group_group__related': '(related_group_id,group_id)',
				}
			}, #.db.group_group
			
			
			'group_biopolymer': {
				'table': """
(
  group_id INTEGER NOT NULL,
  biopolymer_id INTEGER NOT NULL,
  specificity TINYINT NOT NULL,
  implication TINYINT NOT NULL,
  quality TINYINT NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (group_id,biopolymer_id,source_id)
)
""",
				'index': {
					'group_biopolymer__biopolymer': '(biopolymer_id,group_id)',
				}
			}, #.db.group_biopolymer
			
			
			'group_member_name': {
				'table': """
(
  group_id INTEGER NOT NULL,
  member INTEGER NOT NULL,
  type_id TINYINT NOT NULL,
  namespace_id INTEGER NOT NULL,
  name VARCHAR(256) NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (group_id,member,type_id,namespace_id,name)
)
""",
				'index': {}
			}, #.db.group_member_name
			
			
			##################################################
			# liftover tables
			
			
			'grch_ucschg': {
				'table': """
(
  grch INTEGER PRIMARY KEY,
  ucschg INTEGER NOT NULL
)
""",
				# hardcode translations between GRCh/NCBI and UCSC 'hg' genome build numbers
				# TODO: find a source for this so we can transition to new builds without a code update, i.e.
				#         http://genome.ucsc.edu/FAQ/FAQreleases.html
				#         http://genome.ucsc.edu/goldenPath/releaseLog.html
				#       these aren't meant to be machine-readable but might be good enough if they're kept updated
				'data': [
					(34,16),
					(35,17),
					(36,18),
					(37,19),
				],
				'index': {}
			}, #.db.grch_ucschg
			
			
			'chain': {
				'table': """
(
  chain_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  old_ucschg INTEGER NOT NULL,
  old_chr TINYINT NOT NULL,
  old_start BIGINT NOT NULL,
  old_end BIGINT NOT NULL,
  new_ucschg INTEGER NOT NULL,
  new_chr TINYINT NOT NULL,
  new_start BIGINT NOT NULL,
  new_end BIGINT NOT NULL,
  score BIGINT NOT NULL,
  is_fwd TINYINT NOT NULL,
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'chain__oldhg_newhg_chr': '(old_ucschg,new_ucschg,old_chr)',
				}
			}, #.db.chain
			
			
			'chain_data': {
				'table': """
(
  chain_id INTEGER NOT NULL,
  old_start BIGINT NOT NULL,
  old_end BIGINT NOT NULL,
  new_start BIGINT NOT NULL,
  source_id TINYINT NOT NULL,
  PRIMARY KEY (chain_id,old_start)
)
""",
				'index': {
					'chain_data__end': '(chain_id,old_end)',
				}
			}, #.db.chain_data
			
		}, #.db
		
	} #_schema{}
	
	
	##################################################
	# class interrogation
	
	
	@classmethod
	def getVersionString(cls):
		return "%d.%d.%d%s%s (%s)" % (cls.ver_maj, cls.ver_min, cls.ver_rev, ("-" if cls.ver_dev else ""), (cls.ver_dev or ""), cls.ver_date)
	#getVersionString()
	
	
	@classmethod
	def checkMinimumVersion(cls, major=None, minor=None, revision=None, development=None):
		if (major == None) or (cls.ver_maj > major):
			return True
		if cls.ver_maj < major:
			return False
		
		if (minor == None) or (cls.ver_min > minor):
			return True
		if cls.ver_min < minor:
			return False
		
		if (revision == None) or (cls.ver_rev > revision):
			return True
		if cls.ver_rev < revision:
			return False
		
		if (development == None) or (not cls.ver_dev) or (cls.ver_dev > development):
			return True
		if cls.ver_dev < development:
			return False
		
		return True
	#checkMinimumVersion()
	
	
	@classmethod
	def getDatabaseDriverName(cls):
		return "SQLite"
	#getDatabaseDriverName()
	
	
	@classmethod
	def getDatabaseDriverVersion(cls):
		return apsw.sqlitelibversion()
	#getDatabaseDriverVersion()
	
	
	@classmethod
	def getDatabaseInterfaceName(cls):
		return "APSW"
	#getDatabaseInterfaceName()
	
	
	@classmethod
	def getDatabaseInterfaceVersion(cls):
		return apsw.apswversion()
	#getDatabaseInterfaceVersion()
	
	
	##################################################
	# constructor
	
	
	def __init__(self, dbFile=None, testing=False, updating=False):
		# initialize instance properties
		self._is_test = testing
		self._updating = updating
		self._verbose = False
		self._logger = None
		self._logFile = sys.stderr
		self._logIndent = 0
		self._logHanging = False
		self._db = apsw.Connection('')
		self._dbFile = None
		self._dbNew = None
		self._updater = None
		self._liftOverCache = dict() # { (from,to) : [] }
		
		self.configureDatabase()
		self.attachDatabaseFile(dbFile)
	#__init__()
	
	
	##################################################
	# context manager
	
	
	def __enter__(self):
		return self._db.__enter__()
	#__enter__()
	
	
	def __exit__(self, excType, excVal, traceback):
		return self._db.__exit__(excType, excVal, traceback)
	#__exit__()
	
	
	##################################################
	# logging
	
	def _checkTesting(self):
		now_test = self.getDatabaseSetting("testing")
		if now_test is None or bool(int(now_test)) == bool(self._is_test):
			self.setDatabaseSetting("testing", bool(self._is_test))
			return True
		else:
			return False
	# setTesting(is_test)
	
	def getVerbose(self):
		return self._verbose
	#getVerbose()
	
	
	def setVerbose(self, verbose=True):
		self._verbose = verbose
	#setVerbose()
	
	
	def setLogger(self, logger=None):
		self._logger = logger
	#setLogger()
	
	
	def log(self, message=""):
		if self._logger:
			self._logger.log(message)
		elif self._verbose:
			if (self._logIndent > 0) and (not self._logHanging):
				self._logFile.write(self._logIndent * "  ")
				self._logHanging = True
			self._logFile.write(message)
			if (message == "") or (message[-1] != "\n"):
				self._logHanging = True
				self._logFile.flush()
			else:
				self._logHanging = False
		#if _logger
	#log()
	
	
	def logPush(self, message=None):
		if self._logger:
			self._logger.logPush(message)
		else:
			if message:
				self.log(message)
			if self._logHanging:
				self.log("\n")
			self._logIndent += 1
		#if _logger
	#logPush()
	
	
	def logPop(self, message=None):
		if self._logger:
			self._logger.logPop(message)
		else:
			if self._logHanging:
				self.log("\n")
			self._logIndent = max(0, self._logIndent - 1)
			if message:
				self.log(message)
		#if _logger
	#logPop()
	
	
	##################################################
	# database management
	
	
	def configureDatabase(self, db=None, pageSize=4096, cacheSize=-65536, sync='OFF'):
		cursor = self._db.cursor()
		db = ("%s." % db) if db else None
		
		# linux VFS doesn't usually report actual disk cluster size,
		# so sqlite ends up using 1KB pages by default
		cursor.execute("PRAGMA %spage_size = %d" % (db,pageSize))
		
		# cache_size is pages if positive, kilobytes if negative
		cursor.execute("PRAGMA %scache_size = %d" % (db,cacheSize))
		
		# for typical read-only usage, synchronization behavior is moot anyway,
		# and while updating we're not that worried about a power failure
		# corrupting the database file since the user could just start the
		# update over from the beginning; so, we'll take the performance gain
		cursor.execute("PRAGMA %ssynchronous = %s" % (db,sync))
	#configureDatabase()
	
		
	def attachTempDatabase(self, db):
		cursor = self._db.cursor()
		
		# detach the current db, if any
		try:
			cursor.execute("DETACH DATABASE `%s`" % db)
		except apsw.SQLError as e:
			if not str(e).startswith('SQLError: no such database: '):
				raise e
		
		# attach a new temp db
		cursor.execute("ATTACH DATABASE '' AS `%s`" % db)
		self.configureDatabase(db)
	#attachTempDatabase()
	
		
	def attachDatabaseFile(self, dbFile, readOnly=False, quiet=False):
		cursor = self._db.cursor()
		
		# detach the current db file, if any
		if self._dbFile and not quiet:
			self.log("unloading knowledge database file '%s' ..." % self._dbFile)
		try:
			cursor.execute("DETACH DATABASE `db`")
		except apsw.SQLError as e:
			if not str(e).startswith('SQLError: no such database: '):
				raise e
		if self._dbFile and not quiet:
			self.log(" OK\n")
		
		# reset db info
		self._dbFile = None
		self._dbNew = None
		
		# attach the new db file, if any
		if dbFile:
			if not quiet:
				self.logPush("loading knowledge database file '%s' ..." % dbFile)
			cursor.execute("ATTACH DATABASE ? AS `db`", (dbFile,))
			self.configureDatabase('db')
			self._dbFile = dbFile
			self._dbNew = (0 == max(row[0] for row in cursor.execute("SELECT COUNT(1) FROM `db`.`sqlite_master`")))
			
			# establish or audit database schema
			err_msg = ""
			with self._db:
				if self._dbNew:
					self.createDatabaseObjects(None, 'db')
					ok = True
				else:
					ok = self.auditDatabaseObjects(None, 'db')
					if not ok:
						err_msg = "Audit of database failed"
				
				if ok and self._updating:
					ok = self._checkTesting()
					if not ok:
						err_msg = "Testing settings do not match loaded database"
			
			if ok:
				if not quiet:
					self.logPop("... OK\n")
			else:
				self._dbFile = None
				self._dbNew = None
				cursor.execute("DETACH DATABASE `db`")
				if not quiet:
					self.logPop("... ERROR (" + err_msg + ")\n")
		#if new dbFile
	#attachDatabaseFile()
	
	
	def detachDatabaseFile(self, quiet=False):
		return self.attachDatabaseFile(None, quiet=quiet)
	#detachDatabaseFile()
	
	
	def testDatabaseUpdate(self):
		if self._dbFile == None:
			raise Exception("ERROR: no knowledge database file is loaded")
		if int(self.getDatabaseSetting('finalized') or 0):
			raise Exception("ERROR: knowledge database has been finalized")
		try:
			if self._db.readonly('db'):
				raise Exception("ERROR: knowledge database file may not be modified")
		except AttributeError: # apsw.Connection.readonly() added in 3.7.11
			try:
				self._db.cursor().execute("UPDATE `db`.`setting` SET value = value")
			except apsw.ReadOnlyError:
				raise Exception("ERROR: knowledge database file may not be modified")
		return True
	#testDatabaseUpdate()
	
	
	def createDatabaseObjects(self, schema, dbName, tblList=None, doTables=True, idxList=None, doIndecies=True):
		cursor = self._db.cursor()
		schema = schema or self._schema[dbName]
		dbType = "TEMP " if (dbName == "temp") else ""
		if tblList and isinstance(tblList, str):
			tblList = (tblList,)
		if idxList and isinstance(idxList, str):
			idxList = (idxList,)
		for tblName in (tblList or schema.keys()):
			if doTables:
				cursor.execute("CREATE %sTABLE IF NOT EXISTS `%s`.`%s` %s" % (dbType, dbName, tblName, schema[tblName]['table']))
				if 'data' in schema[tblName] and schema[tblName]['data']:
					sql = "INSERT OR IGNORE INTO `%s`.`%s` VALUES (%s)" % (dbName, tblName, ("?,"*len(schema[tblName]['data'][0]))[:-1])
					# TODO: change how 'data' is defined so it can be tested without having to try inserting
					try:
						cursor.executemany(sql, schema[tblName]['data'])
					except apsw.ReadOnlyError:
						pass
			if doIndecies:
				for idxName in (idxList or schema[tblName]['index'].keys()):
					if idxName not in schema[tblName]['index']:
						raise Exception("ERROR: no definition for index '%s' on table '%s'" % (idxName,tblName))
					cursor.execute("CREATE INDEX IF NOT EXISTS `%s`.`%s` ON `%s` %s" % (dbName, idxName, tblName, schema[tblName]['index'][idxName]))
				#foreach idxName in idxList
				cursor.execute("ANALYZE `%s`.`%s`" % (dbName,tblName))
		#foreach tblName in tblList
		
		# this shouldn't be necessary since we don't manually modify the sqlite_stat* tables
		#if doIndecies:
		#	cursor.execute("ANALYZE `%s`.`sqlite_master`" % (dbName,))
	#createDatabaseObjects()
	
	
	def createDatabaseTables(self, schema, dbName, tblList, doIndecies=False):
		return self.createDatabaseObjects(schema, dbName, tblList, True, None, doIndecies)
	#createDatabaseTables()
	
	
	def createDatabaseIndecies(self, schema, dbName, tblList, doTables=False, idxList=None):
		return self.createDatabaseObjects(schema, dbName, tblList, doTables, idxList, True)
	#createDatabaseIndecies()
	
	
	def dropDatabaseObjects(self, schema, dbName, tblList=None, doTables=True, idxList=None, doIndecies=True):
		cursor = self._db.cursor()
		schema = schema or self._schema[dbName]
		if tblList and isinstance(tblList, str):
			tblList = (tblList,)
		if idxList and isinstance(idxList, str):
			idxList = (idxList,)
		for tblName in (tblList or schema.keys()):
			if doTables:
				if dbName == 'db':
					self.testDatabaseUpdate()
				cursor.execute("DROP TABLE IF EXISTS `%s`.`%s`" % (dbName, tblName))
			elif doIndecies:
				for idxName in (idxList or schema[tblName]['index'].keys()):
					cursor.execute("DROP INDEX IF EXISTS `%s`.`%s`" % (dbName, idxName))
				#foreach idxName in idxList
		#foreach tblName in tblList
	#dropDatabaseObjects()
	
	
	def dropDatabaseTables(self, schema, dbName, tblList):
		return self.dropDatabaseObjects(schema, dbName, tblList, True, None, True)
	#dropDatabaseTables()
	
	
	def dropDatabaseIndecies(self, schema, dbName, tblList, idxList=None):
		return self.dropDatabaseObjects(schema, dbName, tblList, False, idxList, True)
	#dropDatabaseIndecies()
	
	
	def auditDatabaseObjects(self, schema, dbName, tblList=None, doTables=True, idxList=None, doIndecies=True, doRepair=True):
		# fetch current schema
		cursor = self._db.cursor()
		current = dict()
		dbMaster = "`sqlite_temp_master`" if (dbName == "temp") else ("`%s`.`sqlite_master`" % dbName)
		for row in cursor.execute("SELECT tbl_name,type,name,sql FROM %s WHERE type IN ('table','index')" % dbMaster):
			tblName,objType,idxName,objDef = row
			if tblName not in current:
				current[tblName] = {'table':None, 'index':{}}
			if objType == 'table':
				current[tblName]['table'] = objDef
			elif objType == 'index':
				current[tblName]['index'][idxName] = objDef
		tblEmpty = dict()
		for tblName in current:
			tblEmpty[tblName] = True
			for row in cursor.execute("SELECT 1 FROM `%s`.`%s` LIMIT 1" % (dbName,tblName)):
				tblEmpty[tblName] = False
		# audit requested objects
		schema = schema or self._schema[dbName]
		if tblList and isinstance(tblList, str):
			tblList = (tblList,)
		if idxList and isinstance(idxList, str):
			idxList = (idxList,)
		ok = True
		for tblName in (tblList or schema.keys()):
			if doTables:
				if tblName in current:
					if current[tblName]['table'].rstrip() == ("CREATE TABLE `%s` %s" % (tblName, schema[tblName]['table'].rstrip())):
						if 'data' in schema[tblName] and schema[tblName]['data']:
							sql = "INSERT OR IGNORE INTO `%s`.`%s` VALUES (%s)" % (dbName, tblName, ("?,"*len(schema[tblName]['data'][0]))[:-1])
							# TODO: change how 'data' is defined so it can be tested without having to try inserting
							try:
								cursor.executemany(sql, schema[tblName]['data'])
							except apsw.ReadOnlyError:
								pass
					elif doRepair and tblEmpty[tblName]:
						self.log("WARNING: table '%s' schema mismatch -- repairing ..." % tblName)
						self.dropDatabaseTables(schema, dbName, tblName)
						self.createDatabaseTables(schema, dbName, tblName, doIndecies)
						self.log(" OK\n")
					elif doRepair:
						self.log("ERROR: table '%s' schema mismatch -- cannot repair\n" % tblName)
						ok = False
					else:
						self.log("ERROR: table '%s' schema mismatch\n" % tblName)
						ok = False
					#if definition match
				elif doRepair:
					self.log("WARNING: table '%s' is missing -- repairing ..." % tblName)
					self.createDatabaseTables(schema, dbName, tblName, doIndecies)
					self.log(" OK\n")
				else:
					self.log("ERROR: table '%s' is missing\n" % tblName)
					ok = False
				#if tblName in current
			#if doTables
			if doIndecies:
				for idxName in (idxList or schema[tblName]['index'].keys()):
					if (tblName not in current) and not (doTables and doRepair):
						self.log("ERROR: table '%s' is missing for index '%s'\n" % (tblName, idxName))
						ok = False
					elif tblName in current and idxName in current[tblName]['index']:
						if current[tblName]['index'][idxName].rstrip() == ("CREATE INDEX `%s` ON `%s` %s" % (idxName, tblName, schema[tblName]['index'][idxName].rstrip())):
							pass
						elif doRepair:
							self.log("WARNING: index '%s' on table '%s' schema mismatch -- repairing ..." % (idxName, tblName))
							self.dropDatabaseIndecies(schema, dbName, tblName, idxName)
							self.createDatabaseIndecies(schema, dbName, tblName, False, idxName)
							self.log(" OK\n")
						else:
							self.log("ERROR: index '%s' on table '%s' schema mismatch\n" % (idxName, tblName))
							ok = False
						#if definition match
					elif doRepair:
						self.log("WARNING: index '%s' on table '%s' is missing -- repairing ..." % (idxName, tblName))
						self.createDatabaseIndecies(schema, dbName, tblName, False, idxName)
						self.log(" OK\n")
					else:
						self.log("ERROR: index '%s' on table '%s' is missing\n" % (idxName, tblName))
						ok = False
					#if tblName,idxName in current
				#foreach idxName in idxList
			#if doIndecies
		#foreach tblName in tblList
		return ok
	#auditDatabaseObjects()
	
		
	def defragmentDatabase(self):
		# unfortunately sqlite's VACUUM doesn't work on attached databases,
		# so we have to detach, make a new direct connection, then re-attach
		if self._dbFile:
			dbFile = self._dbFile
			self.detachDatabaseFile(quiet=True)
			db = apsw.Connection(dbFile)
			db.cursor().execute("VACUUM")
			db.close()
			self.attachDatabaseFile(dbFile, quiet=True)
	#defragmentDatabase()
	
	
	def finalizeDatabase(self):
		self.log("discarding intermediate data ...")
		self.testDatabaseUpdate()
		self.dropDatabaseTables(None, 'db', ('snp_entrez_role','biopolymer_name_name','group_member_name'))
		self.createDatabaseTables(None, 'db', ('snp_entrez_role','biopolymer_name_name','group_member_name'), True)
		self.log(" OK\n")
		self.log("updating optimizer statistics ...")
		self._db.cursor().execute("ANALYZE `db`")
		self.log(" OK\n")
		self.log("compacting knowledge database file (this may take several hours!) ...")
		self.defragmentDatabase()
		self.log(" OK\n")
		self.setDatabaseSetting('finalized', 1)
	#finalizeDatabase()
	
	
	def getDatabaseSetting(self, setting):
		value = None
		if self._dbFile:
			for row in self._db.cursor().execute("SELECT value FROM `db`.`setting` WHERE setting = ?", (setting,)):
				value = row[0]
		return value
	#getDatabaseSetting()
	
	
	def setDatabaseSetting(self, setting, value):
		self.testDatabaseUpdate()
		self._db.cursor().execute("INSERT OR REPLACE INTO `db`.`setting` (setting, value) VALUES (?, ?)", (setting,value))
	#setDatabaseSetting()
	
	
	def getSourceModules(self):
		if not self._updater:
			import loki_updater
			self._updater = loki_updater.Updater(self, self._is_test)
		return self._updater.getSourceModules()
	#getSourceModules()
	
	
	def getSourceModuleVersions(self, sources=None):
		if not self._updater:
			import loki_updater
			self._updater = loki_updater.Updater(self, self._is_test)
		return self._updater.getSourceModuleVersions(sources)
	#getSourceModuleVersions()
	
	
	def getSourceModuleOptions(self, sources=None):
		if not self._updater:
			import loki_updater
			self._updater = loki_updater.Updater(self, self._is_test)
		return self._updater.getSourceModuleOptions(sources)
	#getSourceModuleOptions()
	
	
	def updateDatabase(self, sources=None, sourceOptions=None, cacheOnly=False):
		self.testDatabaseUpdate()
		if not self._updater:
			import loki_updater
			self._updater = loki_updater.Updater(self, self._is_test)
		return self._updater.updateDatabase(sources, sourceOptions, cacheOnly)
	#updateDatabase()
	
	
	def prepareTableForUpdate(self, table):
		if self._updater:
			return self._updater.prepareTableForUpdate(table)
		return None
	#prepareTableForUpdate()
	
	
	def prepareTableForQuery(self, table):
		if self._updater:
			return self._updater.prepareTableForQuery(table)
		return None
	#prepareTableForQuery()
	
	
	##################################################
	# metadata retrieval
	
	
	def getLDProfileID(self, ldprofile):
		return self.getLDProfileIDs([ldprofile])[ldprofile]
	#getLDProfileID()
	
	
	def getLDProfileIDs(self, ldprofiles):
		if not self._dbFile:
			return { l:None for l in ldprofiles }
		sql = "SELECT i.ldprofile, l.ldprofile_id FROM (SELECT ? AS ldprofile) AS i LEFT JOIN `db`.`ldprofile` AS l ON l.ldprofile = LOWER(i.ldprofile)"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(ldprofiles)) }
		return ret
	#getLDProfileIDs()
	
	
	def getNamespaceID(self, namespace):
		return self.getNamespaceIDs([namespace])[namespace]
	#getNamespaceID()
	
	
	def getNamespaceIDs(self, namespaces):
		if not self._dbFile:
			return { n:None for n in namespaces }
		sql = "SELECT i.namespace, n.namespace_id FROM (SELECT ? AS namespace) AS i LEFT JOIN `db`.`namespace` AS n ON n.namespace = LOWER(i.namespace)"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(namespaces)) }
		return ret
	#getNamespaceIDs()
	
	
	def getRelationshipID(self, relationship):
		return self.getRelationshipIDs([relationship])[relationship]
	#getRelationshipID()
	
	
	def getRelationshipIDs(self, relationships):
		if not self._dbFile:
			return { r:None for r in relationships }
		sql = "SELECT i.relationship, r.relationship_id FROM (SELECT ? AS relationship) AS i LEFT JOIN `db`.`relationship` AS r ON r.relationship = LOWER(i.relationship)"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(relationships)) }
		return ret
	#getRelationshipIDs()
	
	
	def getRoleID(self, role):
		return self.getRoleIDs([role])[role]
	#getRoleID()
	
	
	def getRoleIDs(self, roles):
		if not self._dbFile:
			return { r:None for r in roles }
		sql = "SELECT i.role, role_id FROM (SELECT ? AS role) AS i LEFT JOIN `db`.`role` AS r ON r.role = LOWER(i.role)"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(roles)) }
		return ret
	#getRoleIDs()
	
	
	def getSourceID(self, source):
		return self.getSourceIDs([source])[source]
	#getSourceID()
	
	
	def getSourceIDs(self, sources):
		if not self._dbFile:
			return { s:None for s in sources }
		sql = "SELECT i.source, s.source_id FROM (SELECT ? AS source) AS i LEFT JOIN `db`.`source` AS s ON s.source = LOWER(i.source)"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(sources)) }
		return ret
	#getSourceIDs()
	
	
	def getTypeID(self, type):
		return self.getTypeIDs([type])[type]
	#getTypeID()
	
	
	def getTypeIDs(self, types):
		if not self._dbFile:
			return { t:None for t in types }
		sql = "SELECT i.type, t.type_id FROM (SELECT ? AS type) AS i LEFT JOIN `db`.`type` AS t ON t.type = LOWER(i.type)"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(types)) }
		return ret
	#getTypeIDs()
	
	
	##################################################
	# snp data retrieval
	
	
	def generateCurrentRSesByRS(self, rses, tally=None):
		# rses=[ rs, ... ]
		# tally=dict()
		# yield:[ (rsInput,rsCurrent), ... ]
		sql = """
SELECT i.rsMerged, COALESCE(sm.rsCurrent, i.rsMerged) AS rsCurrent
FROM (SELECT ? AS rsMerged) AS i
LEFT JOIN `db`.`snp_merge` AS sm USING (rsMerged)
"""
		with self._db:
			if tally != None:
				numMerge = numMatch = 0
				for row in self._db.cursor().executemany(sql, itertools.izip(rses)):
					if row[1] != row[0]:
						numMerge += 1
					else:
						numMatch += 1
					yield row
				tally['merge'] = numMerge
				tally['match'] = numMatch
			else:
				for row in self._db.cursor().executemany(sql, itertools.izip(rses)):
					yield row
	#generateCurrentRSesByRS()
	
	
	def generateSNPLociByRS(self, rses, minMatch=1, maxMatch=1, validated=None, tally=None):
		# rses=[ rs, ... ]
		# tally=dict()
		# yield:[ (rs,chr,pos), ... ]
		sql = """
SELECT i.rs, sl.chr, sl.pos
FROM (SELECT ? AS rs) AS i
LEFT JOIN `db`.`snp_locus` AS sl
  ON sl.rs = i.rs
"""
		if validated != None:
			sql += "  AND sl.validated = %d" % (1 if validated else 0)
		key = matches = None
		numNull = numAmbig = numMatch = 0
		with self._db:
			for row in itertools.chain(self._db.cursor().executemany(sql, itertools.izip(rses)), [(None,None,None)]):
				if key != row[0]:
					if key:
						if tally != None:
							if not matches:
								numNull += 1
							elif len(matches) > 1:
								numAmbig += 1
							else:
								numMatch += 1
						if (minMatch or 0) < 1 and not matches:
							yield (key,None,None)
						elif (minMatch or 0) <= len(matches) <= (maxMatch or len(matches)):
							for match in matches:
								yield match
					key = row[0]
					matches = set()
				if row[1] and row[2]:
					matches.add(row)
			#foreach row
		if tally != None:
			tally['null'] = numNull
			tally['ambig'] = numAmbig
			tally['match'] = numMatch
	#generateSNPLociByRS()
	
	
	##################################################
	# biopolymer data retrieval
	
	
	def generateBiopolymersByID(self, ids):
		# ids=[ id, ... ]
		# yield:[ (id,type_id,label,description), ... ]
		sql = "SELECT biopolymer_id, type_id, label, description FROM `db`.`biopolymer` WHERE biopolymer_id = ?"
		return self._db.cursor().executemany(sql, itertools.izip(ids))
	#generateBiopolymersByID()
	
	
	def generateBiopolymerIDsByName(self, names, minMatch=1, maxMatch=1, tally=None, namespaceID=None, typeID=None):
		# names=[ name, ... ]
		# tally=dict()
		# namespaceID=0 means to search names using any namespace
		# namespaceID=None means to search primary labels directly
		# yields (name,id)
		
		sql = """
SELECT i.name, b.biopolymer_id
FROM (SELECT ? AS name) AS i
"""
		
		if namespaceID != None:
			sql += """
LEFT JOIN `db`.`biopolymer_name` AS bn
  ON bn.name = i.name
"""
			if namespaceID:
				sql += """
  AND bn.namespace_id = %d
""" % namespaceID
		
		sql += """
LEFT JOIN `db`.`biopolymer` AS b
"""
		
		if namespaceID != None:
			sql += """
  ON b.biopolymer_id = bn.biopolymer_id
"""
		else:
			sql += """
  ON b.label = i.name
"""
		
		if typeID:
			sql += """
  AND b.type_id = %d
""" % typeID
		
		key = matches = None
		numNull = numAmbig = numMatch = 0
		with self._db:
			for row in itertools.chain(self._db.cursor().executemany(sql, itertools.izip(names)), [(None,None)]):
				if key != row[0]:
					if key:
						if tally != None:
							if not matches:
								numNull += 1
							elif len(matches) > 1:
								numAmbig += 1
							else:
								numMatch += 1
						if minMatch < 1 and not matches:
							yield (key,None)
						elif (minMatch or 0) <= len(matches) <= (maxMatch or len(matches)):
							for match in matches:
								yield match
					key = row[0]
					matches = set()
				if row[1]:
					matches.add(row)
			#foreach row
		if tally != None:
			tally['null'] = numNull
			tally['ambig'] = numAmbig
			tally['match'] = numMatch
	#generateBiopolymerIDsByName()
	
	
	def generateBiopolymerNameStats(self, namespaceID=None, typeID=None):
		sql = """
SELECT
  `namespace`,
  COUNT() AS `names`,
  SUM(CASE WHEN matches = 1 THEN 1 ELSE 0 END) AS `unique`,
  SUM(CASE WHEN matches > 1 THEN 1 ELSE 0 END) AS `ambiguous`
FROM (
  SELECT bn.namespace_id, bn.name, COUNT(DISTINCT bn.biopolymer_id) AS matches
  FROM `db`.`biopolymer_name` AS bn
"""
		
		if typeID:
			sql += """
  JOIN `db`.`biopolymer` AS b
    ON b.biopolymer_id = bn.biopolymer_id AND b.type_id = %d
""" % typeID
		
		if namespaceID:
			sql += """
  WHERE bn.namespace_id = %d
""" % namespaceID
		
		sql += """
  GROUP BY bn.namespace_id, bn.name
)
JOIN `db`.`namespace` AS n USING (namespace_id)
GROUP BY namespace_id
"""
		
		for row in self._db.cursor().execute(sql):
			yield row
	#generateBiopolymerNameStats()
	
	
	##################################################
	# group data retrieval
	
	
	def generateGroupsByID(self, ids):
		# ids=[ id, ... ]
		# yield:[ (id,type_id,label,description), ... ]
		sql = "SELECT group_id, type_id, label, description FROM `db`.`group` WHERE group_id = ?"
		return self._db.cursor().executemany(sql, itertools.izip(ids))
	#generateGroupsByID()
	
	
	def generateGroupIDsByName(self, names, minMatch=1, maxMatch=1, tally=None, namespaceID=None, typeID=None):
		# names=[ name, ... ]
		# tally=dict()
		# namespaceID=0 means to search names using any namespace
		# namespaceID=None means to search primary labels directly
		# yields (name,id)
		
		sql = """
SELECT i.name, g.group_id
FROM (SELECT ? AS name) AS i
"""
		
		if namespaceID != None:
			sql += """
LEFT JOIN `db`.`group_name` AS gn
  ON gn.name = i.name
"""
			if namespaceID:
				sql += """
  AND gn.namespace_id = %d
""" % namespaceID
		
		sql += """
LEFT JOIN `db`.`group` AS g
"""
		
		if namespaceID != None:
			sql += """
  ON g.group_id = gn.group_id
"""
		else:
			sql += """
  ON g.label = i.name
"""
		
		if typeID:
			sql += """
  AND g.type_id = %d
""" % typeID
		
		key = matches = None
		numNull = numAmbig = numMatch = 0
		with self._db:
			for row in itertools.chain(self._db.cursor().executemany(sql, itertools.izip(names)), [(None,None)]):
				if key != row[0]:
					if key:
						if tally != None:
							if not matches:
								numNull += 1
							elif len(matches) > 1:
								numAmbig += 1
							else:
								numMatch += 1
						if minMatch < 1 and not matches:
							yield (key,None)
						elif (minMatch or 0) <= len(matches) <= (maxMatch or len(matches)):
							for match in matches:
								yield match
					key = row[0]
					matches = set()
				if row[1]:
					matches.add(row)
			#foreach row
		if tally != None:
			tally['null'] = numNull
			tally['ambig'] = numAmbig
			tally['match'] = numMatch
	#generateGroupIDsByName()
	
	
	def generateGroupNameStats(self, namespaceID=None, typeID=None):
		sql = """
SELECT
  `namespace`,
  COUNT() AS `names`,
  SUM(CASE WHEN matches = 1 THEN 1 ELSE 0 END) AS `unique`,
  SUM(CASE WHEN matches > 1 THEN 1 ELSE 0 END) AS `ambiguous`
FROM (
  SELECT gn.namespace_id, gn.name, COUNT(DISTINCT gn.group_id) AS matches
  FROM `db`.`group_name` AS gn
"""
		
		if typeID:
			sql += """
  JOIN `db`.`group` AS g
    ON g.group_id = gn.group_id AND g.type_id = %d
""" % typeID
		
		if namespaceID:
			sql += """
  WHERE gn.namespace_id = %d
""" % namespaceID
		
		sql += """
  GROUP BY gn.namespace_id, gn.name
)
JOIN `db`.`namespace` AS n USING (namespace_id)
GROUP BY namespace_id
"""
		
		for row in self._db.cursor().execute(sql):
			yield row
	#generateGroupNameStats()
	
	
	##################################################
	# liftover
	# 
	# originally from UCSC
	# reimplemented? in C++ for Biofilter 1.0 by Eric Torstenson
	# ported to Python by John Wallace
	
	
	def _generateApplicableLiftOverChains(self, oldHG, newHG, chrom, start, end):
		conv = (oldHG,newHG)
		if conv in self._liftOverCache:
			chains = self._liftOverCache[conv]
		else:
			chains = {'data':{}, 'keys':{}}
			sql = """
SELECT chain_id,
  c.old_chr, c.score, c.old_start, c.old_end, c.new_start, c.is_fwd, c.new_chr,
  cd.old_start, cd.old_end, cd.new_start
FROM `db`.`chain` AS c
JOIN `db`.`chain_data` AS cd USING (chain_id)
WHERE c.old_ucschg=? AND c.new_ucschg=?
ORDER BY c.old_chr, score DESC, cd.old_start
"""
			for row in self._db.cursor().execute(sql, conv):
				chain = (row[2], row[3], row[4], row[5], row[6], row[7], row[0])
				chr = row[1]
				
				if chr not in chains['data']:
					chains['data'][chr] = {chain: []}
					chains['keys'][chr] = [chain]
				elif chain not in chains['data'][chr]:
					chains['data'][chr][chain] = []
					chains['keys'][chr].append(chain)
				
				chains['data'][chr][chain].append( (row[8],row[9],row[10]) )
			#foreach row
			
			# Sort the chains by score
			for k in chains['keys']:
				chains['keys'][k].sort(reverse=True)
			
			self._liftOverCache[conv] = chains
		#if chains are cached
		
		for c in chains['keys'].get(chrom, []):
			# if the region overlaps the chain...
			if start <= c[2] and end >= c[1]:
				data = chains['data'][chrom][c]
				idx = bisect.bisect(data, (start, sys.maxint, sys.maxint))
				if idx:
					idx = idx-1
				
				if idx < len(data) - 1 and start == data[idx + 1]:
					idx = idx + 1
				
				while idx < len(data) and data[idx][0] < end:
					yield (c[-1], data[idx][0], data[idx][1], data[idx][2], c[4], c[5])
					idx = idx + 1
		#foreach chain
	#_generateApplicableLiftOverChains()
	
	
	def _liftOverRegionUsingChains(self, label, start, end, first_seg, end_seg, total_mapped_sz):
		"""
		Map a region given the 1st and last segment as well as the total mapped size
		"""
		mapped_reg = None
		
		# The front and end differences are the distances from the
		# beginning of the segment.
		
		# The front difference should be >= 0 and <= size of 1st segment
		front_diff = max(0, min(start - first_seg[1], first_seg[2] - first_seg[1]))
		
		# The end difference should be similar, but w/ last
		end_diff = max(0, min(end - end_seg[1], end_seg[2] - end_seg[1]))
		
		# Now, if we are moving forward, we add the difference
		# to the new_start, backward, we subtract
		# Also, at this point, if backward, swap start/end
		if first_seg[4]:
			new_start = first_seg[3] + front_diff
			new_end = end_seg[3] + end_diff
		else:
			new_start = end_seg[3] - end_diff
			new_end = first_seg[3] - front_diff
		
		# old_startHere, detect if we have mapped a sufficient fraction 
		# of the region.  liftOver uses a default of 95%
		mapped_size = total_mapped_sz - front_diff - (end_seg[2] - end_seg[1]) + end_diff
		
		if mapped_size / float(end - start) >= 0.95: # TODO: configurable threshold?
			mapped_reg = (label, first_seg[5], new_start, new_end)
		
		return mapped_reg
	#_liftOverRegionUsingChains()
	
	
	def generateLiftOverRegions(self, oldHG, newHG, regions, tally=None, errorCallback=None):
		# regions=[ (label,chr,posMin,posMax), ... ]
		oldHG = int(oldHG)
		newHG = int(newHG)
		numNull = numLift = 0
		for region in regions:
			label,chrom,start,end = region
			
			# We need to actually lift regions to detect dropped sections
			is_region = True
			
			# If the start and end are swapped, reverse them, please
			if start > end:
				(start, end) = (end, start)
			elif start == end:
				is_region = False
				end = start + 1	
			
			mapped_reg = None
			curr_chain = None
			total_mapped_sz = 0
			first_seg = None
			end_seg = None
			for seg in self._generateApplicableLiftOverChains(oldHG, newHG, chrom, start, end):
				if curr_chain is None:
					curr_chain = seg[0]
					first_seg = seg
					end_seg = seg
					total_mapped_sz = seg[2] - seg[1]
				elif seg[0] != curr_chain:
					mapped_reg = self._liftOverRegionUsingChains(label, start, end, first_seg, end_seg, total_mapped_sz)
					if mapped_reg:
						break
					curr_chain = seg[0] #TODO ?
					first_seg = seg
					end_seg = seg
					total_mapped_sz = seg[2] - seg[1]
				else:
					end_seg = seg
					total_mapped_sz = total_mapped_sz + seg[2] - seg[1]
			
			if not mapped_reg and first_seg is not None:
				mapped_reg = self._liftOverRegionUsingChains(label, start, end, first_seg, end_seg, total_mapped_sz)
			
			if mapped_reg:
				numLift += 1
				if not is_region:
					mapped_reg = (mapped_reg[0], mapped_reg[1], mapped_reg[2], mapped_reg[2])
				yield mapped_reg
			else:
				numNull += 1
				if errorCallback:
					errorCallback(region)
		#foreach region
		
		if tally != None:
			tally['null'] = numNull
			tally['lift'] = numLift
	#generateLiftOverRegions()
	
	
	def generateLiftOverLoci(oldHG, newHG, loci, tally=None, errorCallback=None):
		# loci=[ (label,chr,pos), ... ]
		regions = ((l[0],l[1],l[2],l[2]) for l in loci)
		for r in self.generateLiftOverRegions(oldHG, newHG, regions, tally, errorCallback):
			yield r[0:3]
	#generateLiftOverLoci()
	
	
#Database


#TODO: temporary liftover testing code!
if __name__ == "__main__":
	inputFile = file(sys.argv[1])
	loki = Database(sys.argv[2])
	outputFile = file(sys.argv[3],'w')
	unmapFile = file(sys.argv[4],'w')
	oldHG = int(sys.argv[5]) if (len(sys.argv) > 5) else 18
	newHG = int(sys.argv[6]) if (len(sys.argv) > 6) else 19
	
	def generateInput():
		for line in inputFile:
			chrom,start,end = line.split()
			if chrom[:3].upper() in ('CHM','CHR'):
				chrom = chrom[3:]
			yield (None, loki.chr_num.get(chrom,-1), int(start), int(end))
	
	def errorCallback(region):
		print >> unmapFile, "chr"+loki.chr_name.get(region[1],'?'), region[2], region[3]
	
	for region in loki.generateLiftOverRegions(oldHG, newHG, generateInput(), errorCallback=errorCallback):
		print >> outputFile, "chr"+loki.chr_name.get(region[1],'?'), region[2], region[3]
