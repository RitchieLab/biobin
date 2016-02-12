#!/usr/bin/env python

import apsw
import bisect
import itertools
import sys


class Database(object):
	
	
	##################################################
	# class interrogation
	
	
	@classmethod
	def getVersionTuple(cls):
		# tuple = (major,minor,revision,dev,build,date)
		# dev must be in ('a','b','rc','release') for lexicographic comparison
		return (2,3,0,'a',1,'2015-12-14')
	#getVersionTuple()
	
	
	@classmethod
	def getVersionString(cls):
		v = list(cls.getVersionTuple())
		# tuple = (major,minor,revision,dev,build,date)
		# dev must be > 'rc' for releases for lexicographic comparison,
		# but we don't need to actually print 'release' in the version string
		v[3] = '' if v[3] > 'rc' else v[3]
		return "%d.%d.%d%s%s (%s)" % tuple(v)
	#getVersionString()
	
	
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
	# public class data
	
	
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
					('schema','3'),
					('ucschg',None),
					('zone_size','100000'),
					('optimized','0'),
					('finalized','0'),
				],
				'index': {}
			}, #.db.setting
			
			
			##################################################
			# metadata tables
			
			
			'grch_ucschg': {
				'table': """
(
  grch INTEGER PRIMARY KEY,
  ucschg INTEGER NOT NULL
)
""",
				# translations known at time of writing are still provided,
				# but additional translations will also be fetched at update
				'data': [
					(34,16),
					(35,17),
					(36,18),
					(37,19),
					(38,38),
				],
				'index': {}
			}, #.db.grch_ucschg
			
			
			'ldprofile': {
				'table': """
(
  ldprofile_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ldprofile VARCHAR(32) UNIQUE NOT NULL,
  description VARCHAR(128),
  metric VARCHAR(32),
  value DOUBLE
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
			
			
			'source_option': {
				'table': """
(
  source_id TINYINT NOT NULL,
  option VARCHAR(32) NOT NULL,
  value VARCHAR(64),
  PRIMARY KEY (source_id, option)
)
""",
				'index': {}
			}, #.db.source_option
			
			
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
  rsMerged INTEGER NOT NULL,
  rsCurrent INTEGER NOT NULL,
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'snp_merge__merge_current': '(rsMerged,rsCurrent)',
				}
			}, #.db.snp_merge
			
			
			'snp_locus': { # all coordinates in LOKI are 1-based closed intervals
				'table': """
(
  rs INTEGER NOT NULL,
  chr TINYINT NOT NULL,
  pos BIGINT NOT NULL,
  validated TINYINT NOT NULL,
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'snp_locus__rs_chr_pos': '(rs,chr,pos)',
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
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'snp_entrez_role__rs_entrez_role': '(rs,entrez_id,role_id)',
				}
			}, #.db.snp_entrez_role
			
			
			'snp_biopolymer_role': {
				'table': """
(
  rs INTEGER NOT NULL,
  biopolymer_id INTEGER NOT NULL,
  role_id INTEGER NOT NULL,
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'snp_biopolymer_role__rs_biopolymer_role': '(rs,biopolymer_id,role_id)',
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
			
			
			'biopolymer_region': { # all coordinates in LOKI are 1-based closed intervals
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
  contains TINYINT,
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
			# gwas tables
			
			
			'gwas': { # all coordinates in LOKI are 1-based closed intervals
				'table': """
(
  gwas_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  rs INTEGER,
  chr TINYINT,
  pos BIGINT,
  trait VARCHAR(256) NOT NULL,
  snps VARCHAR(256),
  orbeta VARCHAR(8),
  allele95ci VARCHAR(16),
  riskAfreq VARCAHR(16),
  pubmed_id INTEGER,
  source_id TINYINT NOT NULL
)
""",
				'index': {
					'gwas__rs': '(rs)',
					'gwas__chr_pos': '(chr,pos)',
				}
			}, #.db.gwas
			
			
			##################################################
			# liftover tables
			
			
			'chain': { # all coordinates in LOKI are 1-based closed intervals
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
			
			
			'chain_data': { # all coordinates in LOKI are 1-based closed intervals
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
	# constructor
	
	
	def __init__(self, dbFile=None, testing=False, updating=False, tempMem=False):
		# initialize instance properties
		self._is_test = testing
		self._updating = updating
		self._verbose = True
		self._logger = None
		self._logFile = sys.stderr
		self._logIndent = 0
		self._logHanging = False
		self._db = apsw.Connection('')
		self._dbFile = None
		self._dbNew = None
		self._updater = None
		self._liftOverCache = dict() # { (from,to) : [] }
		
		self.configureDatabase(tempMem=tempMem)
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
			return self._logger.log(message)
		if self._verbose:
			if (self._logIndent > 0) and (not self._logHanging):
				self._logFile.write(self._logIndent * "  ")
				self._logHanging = True
			self._logFile.write(message)
			if (message == "") or (message[-1] != "\n"):
				self._logHanging = True
				self._logFile.flush()
			else:
				self._logHanging = False
		return self._logIndent
	#log()
	
	
	def logPush(self, message=None):
		if self._logger:
			return self._logger.logPush(message)
		if message:
			self.log(message)
		if self._logHanging:
			self.log("\n")
		self._logIndent += 1
		return self._logIndent
	#logPush()
	
	
	def logPop(self, message=None):
		if self._logger:
			return self._logger.logPop(message)
		if self._logHanging:
			self.log("\n")
		self._logIndent = max(0, self._logIndent - 1)
		if message:
			self.log(message)
		return self._logIndent
	#logPop()
	
	
	##################################################
	# database management
	
	
	def getDatabaseMemoryUsage(self, resetPeak=False):
		return (apsw.memoryused(), apsw.memoryhighwater(resetPeak))
	#getDatabaseMemoryUsage()
	
	
	def getDatabaseMemoryLimit(self):
		return apsw.softheaplimit(-1)
	#getDatabaseMemoryLimit()
	
	
	def setDatabaseMemoryLimit(self, limit=0):
		apsw.softheaplimit(limit)
	#setDatabaseMemoryLimit()
	
	
	def configureDatabase(self, db=None, tempMem=False):
		cursor = self._db.cursor()
		db = ("%s." % db) if db else ""
		
		# linux VFS doesn't usually report actual disk cluster size,
		# so sqlite ends up using 1KB pages by default; we prefer 4KB
		cursor.execute("PRAGMA %spage_size = 4096" % (db,))
		
		# cache_size is pages if positive, kibibytes if negative;
		# seems to only affect write performance
		cursor.execute("PRAGMA %scache_size = -32768" % (db,))
		
		# for typical read-only usage, synchronization behavior is moot anyway,
		# and while updating we're not that worried about a power failure
		# corrupting the database file since the user could just start the
		# update over from the beginning; so, we'll take the performance gain
		cursor.execute("PRAGMA %ssynchronous = OFF" % (db,))
		
		# the journal isn't that big, so keeping it in memory is faster; the
		# cost is that a system crash will corrupt the database rather than
		# leaving it recoverable with the on-disk journal (a program crash
		# should be fine since sqlite will rollback transactions before exiting)
		cursor.execute("PRAGMA %sjournal_mode = MEMORY" % (db,))
		
		# the temp store is used for all of sqlite's internal scratch space
		# needs, such as the TEMP database, indexing, etc; keeping it in memory
		# is much faster, but it can get quite large
		if tempMem and not db:
			cursor.execute("PRAGMA temp_store = MEMORY")
		
		# we want EXCLUSIVE while updating since the data shouldn't be read
		# until ready and we want the performance gain; for normal read usage,
		# NORMAL is better so multiple users can share a database file
		cursor.execute("PRAGMA %slocking_mode = %s" % (db,("EXCLUSIVE" if self._updating else "NORMAL")))
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
	
	
	def attachDatabaseFile(self, dbFile, quiet=False):
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
			self._dbFile = dbFile
			self._dbNew = (0 == max(row[0] for row in cursor.execute("SELECT COUNT(1) FROM `db`.`sqlite_master`")))
			self.configureDatabase('db')
			
			# establish or audit database schema
			err_msg = ""
			with self._db:
				if self._dbNew:
					self.createDatabaseObjects(None, 'db')
					ok = True
				else:
					self.updateDatabaseSchema()
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
	
	
	def testDatabaseWriteable(self):
		if self._dbFile == None:
			raise Exception("ERROR: no knowledge database file is loaded")
		try:
			if self._db.readonly('db'):
				raise Exception("ERROR: knowledge database file cannot be modified")
		except AttributeError: # apsw.Connection.readonly() added in 3.7.11
			try:
				self._db.cursor().execute("UPDATE `db`.`setting` SET value = value")
			except apsw.ReadOnlyError:
				raise Exception("ERROR: knowledge database file cannot be modified")
		return True
	#testDatabaseWriteable()
	
	
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
	
	
	def updateDatabaseSchema(self):
		cursor = self._db.cursor()
		
		if self.getDatabaseSetting('schema',int) < 2:
			self.logPush("updating database schema to version 2 ...\n")
			updateMap = {
				'snp_merge'           : 'rsMerged,rsCurrent,source_id',
				'snp_locus'           : 'rs,chr,pos,validated,source_id',
				'snp_entrez_role'     : 'rs,entrez_id,role_id,source_id',
				'snp_biopolymer_role' : 'rs,biopolymer_id,role_id,source_id',
			}
			for tblName,tblColumns in updateMap.iteritems():
				self.log("%s ..." % (tblName,))
				cursor.execute("ALTER TABLE `db`.`%s` RENAME TO `___old_%s___`" % (tblName,tblName))
				self.createDatabaseTables(None, 'db', tblName)
				cursor.execute("INSERT INTO `db`.`%s` (%s) SELECT %s FROM `db`.`___old_%s___`" % (tblName,tblColumns,tblColumns,tblName))
				cursor.execute("DROP TABLE `db`.`___old_%s___`" % (tblName,))
				self.createDatabaseIndecies(None, 'db', tblName)
				self.log(" OK\n")
			self.setDatabaseSetting('schema', 2)
			self.logPop("... OK\n")
		#schema<2
		
		if self.getDatabaseSetting('schema',int) < 3:
			self.log("updating database schema to version 3 ...")
			self.setDatabaseSetting('optimized', self.getDatabaseSetting('finalized',int))
			self.setDatabaseSetting('schema', 3)
			self.log(" OK\n")
		#schema<3
	#updateDatabaseSchema()
	
	
	def auditDatabaseObjects(self, schema, dbName, tblList=None, doTables=True, idxList=None, doIndecies=True, doRepair=True):
		# fetch current schema
		cursor = self._db.cursor()
		current = dict()
		dbMaster = "`sqlite_temp_master`" if (dbName == "temp") else ("`%s`.`sqlite_master`" % (dbName,))
		sql = "SELECT tbl_name,type,name,COALESCE(sql,'') FROM %s WHERE type IN ('table','index')" % (dbMaster,)
		for row in cursor.execute(sql):
			tblName,objType,idxName,objDef = row
			if tblName not in current:
				current[tblName] = {'table':None, 'index':{}}
			if objType == 'table':
				current[tblName]['table'] = " ".join(objDef.strip().split())
			elif objType == 'index':
				current[tblName]['index'][idxName] = " ".join(objDef.strip().split())
		tblEmpty = dict()
		sql = None
		for tblName in current:
			tblEmpty[tblName] = True
			sql = "SELECT 1 FROM `%s`.`%s` LIMIT 1" % (dbName,tblName)
			for row in cursor.execute(sql):
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
					if current[tblName]['table'] == ("CREATE TABLE `%s` %s" % (tblName, " ".join(schema[tblName]['table'].strip().split()))):
						if 'data' in schema[tblName] and schema[tblName]['data']:
							sql = u"INSERT OR IGNORE INTO `%s`.`%s` VALUES (%s)" % (dbName, tblName, ("?,"*len(schema[tblName]['data'][0]))[:-1])
							# TODO: change how 'data' is defined so it can be tested without having to try inserting
							try:
								cursor.executemany(sql, schema[tblName]['data'])
							except apsw.ReadOnlyError:
								pass
					elif doRepair and tblEmpty[tblName]:
						self.log("WARNING: table '%s' schema mismatch -- repairing ..." % tblName)
						self.dropDatabaseTables(schema, dbName, tblName)
						self.createDatabaseTables(schema, dbName, tblName)
						current[tblName]['index'] = dict()
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
						if current[tblName]['index'][idxName] == ("CREATE INDEX `%s` ON `%s` %s" % (idxName, tblName, " ".join(schema[tblName]['index'][idxName].strip().split()))):
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
	
	
	def finalizeDatabase(self):
		self.log("discarding intermediate data ...")
		self.dropDatabaseTables(None, 'db', ('snp_entrez_role','biopolymer_name_name','group_member_name'))
		self.createDatabaseTables(None, 'db', ('snp_entrez_role','biopolymer_name_name','group_member_name'), True)
		self.log(" OK\n")
		self.setDatabaseSetting('finalized', 1)
		self.setDatabaseSetting('optimized', 0)
	#finalizeDatabase()
	
	
	def optimizeDatabase(self):
		self.log("updating optimizer statistics ...")
		self._db.cursor().execute("ANALYZE `db`")
		self.log(" OK\n")
		self.log("compacting knowledge database file ...")
		self.defragmentDatabase()
		self.log(" OK\n")
		self.setDatabaseSetting('optimized', 1)
	#optimizeDatabase()
	
	
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
	
	
	def getDatabaseSetting(self, setting, type=None):
		value = None
		if self._dbFile:
			for row in self._db.cursor().execute("SELECT value FROM `db`.`setting` WHERE setting = ?", (setting,)):
				value = row[0]
		if type:
			value = type(value) if (value != None) else type()
		return value
	#getDatabaseSetting()
	
	
	def setDatabaseSetting(self, setting, value):
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
	
	
	def updateDatabase(self, sources=None, sourceOptions=None, cacheOnly=False, forceUpdate=False):
		if self.getDatabaseSetting('finalized',int):
			raise Exception("ERROR: cannot update a finalized database")
		if not self._updater:
			import loki_updater
			self._updater = loki_updater.Updater(self, self._is_test)
		return self._updater.updateDatabase(sources, sourceOptions, cacheOnly, forceUpdate)
	#updateDatabase()
	
	
	def prepareTableForUpdate(self, table):
		if self.getDatabaseSetting('finalized',int):
			raise Exception("ERROR: cannot update a finalized database")
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
	
	
	def generateGRChByUCSChg(self, ucschg):
		return (row[0] for row in self._db.cursor().execute("SELECT grch FROM grch_ucschg WHERE ucschg = ?", (ucschg,)))
	#generateGRChByUCSChg()
	
	
	def getUCSChgByGRCh(self, grch):
		ucschg = None
		for row in self._db.cursor().execute("SELECT ucschg FROM grch_ucschg WHERE grch = ?", (grch,)):
			ucschg = row[0]
		return ucschg
	#getUCSChgByGRCh()
	
	
	def getLDProfileID(self, ldprofile):
		return self.getLDProfileIDs([ldprofile])[ldprofile]
	#getLDProfileID()
	
	
	def getLDProfileIDs(self, ldprofiles):
		if not self._dbFile:
			return { l:None for l in ldprofiles }
		sql = "SELECT i.ldprofile, l.ldprofile_id FROM (SELECT ? AS ldprofile) AS i LEFT JOIN `db`.`ldprofile` AS l ON LOWER(TRIM(l.ldprofile)) = LOWER(TRIM(i.ldprofile))"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(ldprofiles)) }
		return ret
	#getLDProfileIDs()
	
	
	def getLDProfiles(self, ldprofiles=None):
		if not self._dbFile:
			return { l:None for l in (ldprofiles or list()) }
		with self._db:
			if ldprofiles:
				sql = "SELECT i.ldprofile, l.ldprofile_id, l.description, l.metric, l.value FROM (SELECT ? AS ldprofile) AS i LEFT JOIN `db`.`ldprofile` AS l ON LOWER(TRIM(l.ldprofile)) = LOWER(TRIM(i.ldprofile))"
				ret = { row[0]:row[1:] for row in self._db.cursor().executemany(sql, itertools.izip(ldprofiles)) }
			else:
				sql = "SELECT l.ldprofile, l.ldprofile_id, l.description, l.metric, l.value FROM `db`.`ldprofile` AS l"
				ret = { row[0]:row[1:] for row in self._db.cursor().execute(sql) }
		return ret
	#getLDProfiles()
	
	
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
	
	
	def getSourceIDs(self, sources=None):
		if not self._dbFile:
			return { s:None for s in (sources or list()) }
		if sources:
			sql = "SELECT i.source, s.source_id FROM (SELECT ? AS source) AS i LEFT JOIN `db`.`source` AS s ON s.source = LOWER(i.source)"
			with self._db:
				ret = { row[0]:row[1] for row in self._db.cursor().executemany(sql, itertools.izip(sources)) }
		else:
			sql = "SELECT source, source_id FROM `db`.`source`"
			with self._db:
				ret = { row[0]:row[1] for row in self._db.cursor().execute(sql) }
		return ret
	#getSourceIDs()
	
	
	def getSourceIDVersion(self, sourceID):
		sql = "SELECT version FROM `db`.`source` WHERE source_id = ?"
		ret = None
		with self._db:
			for row in self._db.cursor().execute(sql, (sourceID,)):
				ret = row[0]
		return ret
	#getSourceIDVersion()
	
	
	def getSourceIDOptions(self, sourceID):
		sql = "SELECT option, value FROM `db`.`source_option` WHERE source_id = ?"
		with self._db:
			ret = { row[0]:row[1] for row in self._db.cursor().execute(sql, (sourceID,)) }
		return ret
	#getSourceIDOptions()
	
	
	def getSourceIDFiles(self, sourceID):
		sql = "SELECT filename, COALESCE(modified,''), COALESCE(size,''), COALESCE(md5,'') FROM `db`.`source_file` WHERE source_id = ?"
		with self._db:
			ret = { row[0]:tuple(row[1:]) for row in self._db.cursor().execute(sql, (sourceID,)) }
		return ret
	#getSourceIDFiles()
	
	
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
	
	
	def generateCurrentRSesByRSes(self, rses, tally=None):
		# rses=[ (rsInput,extra), ... ]
		# tally=dict()
		# yield:[ (rsInput,extra,rsCurrent), ... ]
		sql = """
SELECT i.rsMerged, i.extra, COALESCE(sm.rsCurrent, i.rsMerged) AS rsCurrent
FROM (SELECT ? AS rsMerged, ? AS extra) AS i
LEFT JOIN `db`.`snp_merge` AS sm USING (rsMerged)
"""
		with self._db:
			if tally != None:
				numMerge = numMatch = 0
				for row in self._db.cursor().executemany(sql, rses):
					if row[2] != row[0]:
						numMerge += 1
					else:
						numMatch += 1
					yield row
				tally['merge'] = numMerge
				tally['match'] = numMatch
			else:
				for row in self._db.cursor().executemany(sql, rses):
					yield row
	#generateCurrentRSesByRSes()
	
	
	def generateSNPLociByRSes(self, rses, minMatch=1, maxMatch=1, validated=None, tally=None, errorCallback=None):
		# rses=[ (rs,extra), ... ]
		# tally=dict()
		# yield:[ (rs,extra,chr,pos), ... ]
		sql = """
SELECT i.rs, i.extra, sl.chr, sl.pos
FROM (SELECT ? AS rs, ? AS extra) AS i
LEFT JOIN `db`.`snp_locus` AS sl
  ON sl.rs = i.rs
ORDER BY sl.chr, sl.pos
"""
		if validated != None:
			sql += "  AND sl.validated = %d" % (1 if validated else 0)
		
		minMatch = int(minMatch) if (minMatch != None) else 0
		maxMatch = int(maxMatch) if (maxMatch != None) else None
		tag = matches = None
		n = numZero = numOne = numMany = 0
		with self._db:
			for row in itertools.chain(self._db.cursor().executemany(sql, rses), [(None,None,None,None)]):
				if tag != row[0:2]:
					if tag:
						if not matches:
							numZero += 1
						elif len(matches) == 1:
							numOne += 1
						else:
							numMany += 1
						
						if minMatch <= len(matches) <= (maxMatch if (maxMatch != None) else len(matches)):
							for match in (matches or [tag+(None,None)]):
								yield match
						elif errorCallback:
							errorCallback("\t".join((t or "") for t in tag), "%s match%s at index %d" % ((len(matches) or "no"),("" if len(matches) == 1 else "es"),n))
					tag = row[0:2]
					matches = list()
					n += 1
				if row[2] and row[3]:
					matches.append(row)
			#foreach row
		if tally != None:
			tally['zero'] = numZero
			tally['one']  = numOne
			tally['many'] = numMany
	#generateSNPLociByRSes()
	
	
	##################################################
	# biopolymer data retrieval
	
	
	def generateBiopolymersByIDs(self, ids):
		# ids=[ (id,extra), ... ]
		# yield:[ (id,extra,type_id,label,description), ... ]
		sql = "SELECT biopolymer_id, ?2 AS extra, type_id, label, description FROM `db`.`biopolymer` WHERE biopolymer_id = ?1"
		return self._db.cursor().executemany(sql, ids)
	#generateBiopolymersByIDs()
	
	
	def _lookupBiopolymerIDs(self, typeID, identifiers, minMatch, maxMatch, tally, errorCallback):
		# typeID=int or Falseish for any
		# identifiers=[ (namespace,name,extra), ... ]
		#   namespace='' or '*' for any, '-' for labels, '=' for biopolymer_id
		# minMatch=int or Falseish for none
		# maxMatch=int or Falseish for none
		# tally=dict() or None
		# errorCallback=callable(position,input,error)
		# yields (namespace,name,extra,id)
		
		sql = """
SELECT i.namespace, i.identifier, i.extra, COALESCE(bID.biopolymer_id,bLabel.biopolymer_id,bName.biopolymer_id) AS biopolymer_id
FROM (SELECT ?1 AS namespace, ?2 AS identifier, ?3 AS extra) AS i
LEFT JOIN `db`.`biopolymer` AS bID
  ON i.namespace = '='
  AND bID.biopolymer_id = 1*i.identifier
  AND ( ({0} IS NULL) OR (bID.type_id = {0}) )
LEFT JOIN `db`.`biopolymer` AS bLabel
  ON i.namespace = '-'
  AND bLabel.label = i.identifier
  AND ( ({0} IS NULL) OR (bLabel.type_id = {0}) )
LEFT JOIN `db`.`namespace` AS n
  ON i.namespace NOT IN ('=','-')
  AND n.namespace = COALESCE(NULLIF(NULLIF(LOWER(TRIM(i.namespace)),''),'*'),n.namespace)
LEFT JOIN `db`.`biopolymer_name` AS bn
  ON i.namespace NOT IN ('=','-')
  AND bn.name = i.identifier
  AND bn.namespace_id = n.namespace_id
LEFT JOIN `db`.`biopolymer` AS bName
  ON i.namespace NOT IN ('=','-')
  AND bName.biopolymer_id = bn.biopolymer_id
  AND ( ({0} IS NULL) OR (bName.type_id = {0}) )
""".format(int(typeID) if typeID else "NULL")
		
		minMatch = int(minMatch) if (minMatch != None) else 0
		maxMatch = int(maxMatch) if (maxMatch != None) else None
		tag = matches = None
		n = numZero = numOne = numMany = 0
		with self._db:
			for row in itertools.chain(self._db.cursor().executemany(sql, identifiers), [(None,None,None,None)]):
				if tag != row[0:3]:
					if tag:
						if not matches:
							numZero += 1
						elif len(matches) == 1:
							numOne += 1
						else:
							numMany += 1
						
						if minMatch <= len(matches) <= (maxMatch if (maxMatch != None) else len(matches)):
							for match in (matches or [tag+(None,)]):
								yield match
						elif errorCallback:
							errorCallback("\t".join((t or "") for t in tag), "%s match%s at index %d" % ((len(matches) or "no"),("" if len(matches) == 1 else "es"),n))
					tag = row[0:3]
					matches = set()
					n += 1
				if row[3]:
					matches.add(row)
			#foreach row
		if tally != None:
			tally['zero'] = numZero
			tally['one']  = numOne
			tally['many'] = numMany
	#_lookupBiopolymerIDs()
	
	
	def generateBiopolymerIDsByIdentifiers(self, identifiers, minMatch=1, maxMatch=1, tally=None, errorCallback=None):
		# identifiers=[ (namespace,name,extra), ... ]
		return self._lookupBiopolymerIDs(None, identifiers, minMatch, maxMatch, tally, errorCallback)
	#generateBiopolymerIDsByIdentifiers()
	
	
	def generateTypedBiopolymerIDsByIdentifiers(self, typeID, identifiers, minMatch=1, maxMatch=1, tally=None, errorCallback=None):
		# identifiers=[ (namespace,name,extra), ... ]
		return self._lookupBiopolymerIDs(typeID, identifiers, minMatch, maxMatch, tally, errorCallback)
	#generateTypedBiopolymerIDsByIdentifiers()
	
	
	def _searchBiopolymerIDs(self, typeID, texts):
		# texts=[ (text,extra), ... ]
		# yields (extra,label,id)
		
		sql = """
SELECT ?2 AS extra, b.label, b.biopolymer_id
FROM `db`.`biopolymer` AS b
LEFT JOIN `db`.`biopolymer_name` AS bn USING (biopolymer_id)
WHERE
  (
    b.label LIKE '%'||?1||'%'
    OR b.description LIKE '%'||?1||'%'
    OR bn.name LIKE '%'||?1||'%'
  )
"""
		
		if typeID:
			sql += """
  AND b.type_id = %d
""" % typeID
		#if typeID
		
		sql += """
GROUP BY b.biopolymer_id
"""
		
		return self._db.cursor().executemany(sql, texts)
	#_searchBiopolymerIDs()
	
	
	def generateBiopolymerIDsBySearch(self, searches):
		# searches=[ (text,extra), ... ]
		return self._searchBiopolymerIDs(None, searches)
	#generateBiopolymerIDsBySearch()
	
	
	def generateTypedBiopolymerIDsBySearch(self, typeID, searches):
		# searches=[ (text,extra), ... ]
		return self._searchBiopolymerIDs(typeID, searches)
	#generateTypedBiopolymerIDsBySearch()
	
	
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
	
	
	def generateGroupsByIDs(self, ids):
		# ids=[ (id,extra), ... ]
		# yield:[ (id,extra,type_id,label,description), ... ]
		sql = "SELECT group_id, ?2 AS extra, type_id, label, description FROM `db`.`group` WHERE group_id = ?1"
		return self._db.cursor().executemany(sql, ids)
	#generateGroupsByIDs()
	
	
	def _lookupGroupIDs(self, typeID, identifiers, minMatch, maxMatch, tally, errorCallback):
		# typeID=int or Falseish for any
		# identifiers=[ (namespace,name,extra), ... ]
		#   namespace='' or '*' for any, '-' for labels, '=' for group_id
		# minMatch=int or Falseish for none
		# maxMatch=int or Falseish for none
		# tally=dict() or None
		# errorCallback=callable(input,error)
		# yields (namespace,name,extra,id)
		
		sql = """
SELECT i.namespace, i.identifier, i.extra, COALESCE(gID.group_id,gLabel.group_id,gName.group_id) AS group_id
FROM (SELECT ?1 AS namespace, ?2 AS identifier, ?3 AS extra) AS i
LEFT JOIN `db`.`group` AS gID
  ON i.namespace = '='
  AND gID.group_id = 1*i.identifier
  AND ( ({0} IS NULL) OR (gID.type_id = {0}) )
LEFT JOIN `db`.`group` AS gLabel
  ON i.namespace = '-'
  AND gLabel.label = i.identifier
  AND ( ({0} IS NULL) OR (gLabel.type_id = {0}) )
LEFT JOIN `db`.`namespace` AS n
  ON i.namespace NOT IN ('=','-')
  AND n.namespace = COALESCE(NULLIF(NULLIF(LOWER(TRIM(i.namespace)),''),'*'),n.namespace)
LEFT JOIN `db`.`group_name` AS gn
  ON i.namespace NOT IN ('=','-')
  AND gn.name = i.identifier
  AND gn.namespace_id = n.namespace_id
LEFT JOIN `db`.`group` AS gName
  ON i.namespace NOT IN ('=','-')
  AND gName.group_id = gn.group_id
  AND ( ({0} IS NULL) OR (gName.type_id = {0}) )
""".format(int(typeID) if typeID else "NULL")
		
		minMatch = int(minMatch) if (minMatch != None) else 0
		maxMatch = int(maxMatch) if (maxMatch != None) else None
		tag = matches = None
		n = numZero = numOne = numMany = 0
		with self._db:
			for row in itertools.chain(self._db.cursor().executemany(sql, identifiers), [(None,None,None,None)]):
				if tag != row[0:3]:
					if tag:
						if not matches:
							numZero += 1
						elif len(matches) == 1:
							numOne += 1
						else:
							numMany += 1
						
						if minMatch <= len(matches) <= (maxMatch if (maxMatch != None) else len(matches)):
							for match in (matches or [tag+(None,)]):
								yield match
						elif errorCallback:
							errorCallback("\t".join((t or "") for t in tag), "%s match%s at index %d" % ((len(matches) or "no"),("" if len(matches) == 1 else "es"),n))
					tag = row[0:3]
					matches = set()
					n += 1
				if row[3]:
					matches.add(row)
			#foreach row
		if tally != None:
			tally['zero'] = numZero
			tally['one']  = numOne
			tally['many'] = numMany
	#_lookupGroupIDs()
	
	
	def generateGroupIDsByIdentifiers(self, identifiers, minMatch=1, maxMatch=1, tally=None, errorCallback=None):
		# identifiers=[ (namespace,name,extra), ... ]
		return self._lookupGroupIDs(None, identifiers, minMatch, maxMatch, tally, errorCallback)
	#generateGroupIDsByIdentifiers()
	
	
	def generateTypedGroupIDsByIdentifiers(self, typeID, identifiers, minMatch=1, maxMatch=1, tally=None, errorCallback=None):
		# identifiers=[ (namespace,name,extra), ... ]
		return self._lookupGroupIDs(typeID, identifiers, minMatch, maxMatch, tally, errorCallback)
	#generateTypedGroupIDsByIdentifiers()
	
	
	def _searchGroupIDs(self, typeID, texts):
		# texts=[ (text,extra), ... ]
		# yields (extra,label,id)
		
		sql = """
SELECT ?2 AS extra, g.label, g.group_id
FROM `db`.`group` AS g
LEFT JOIN `db`.`group_name` AS gn USING (group_id)
WHERE
  (
    g.label LIKE '%'||?1||'%'
    OR g.description LIKE '%'||?1||'%'
    OR gn.name LIKE '%'||?1||'%'
  )
"""
		
		if typeID:
			sql += """
  AND g.type_id = %d
""" % typeID
		#if typeID
		
		sql += """
GROUP BY g.group_id
"""
		
		return self._db.cursor().executemany(sql, texts)
	#_searchGroupIDs()
	
	
	def generateGroupIDsBySearch(self, searches):
		# searches=[ (text,extra), ... ]
		return self._searchGroupIDs(None, searches)
	#generateGroupIDsBySearch()
	
	
	def generateTypedGroupIDsBySearch(self, typeID, searches):
		# searches=[ (text,extra), ... ]
		return self._searchGroupIDs(typeID, searches)
	#generateTypedGroupIDsBySearch()
	
	
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
	# reimplemented again in Python by John Wallace
	
	
	def hasLiftOverChains(self, oldHG, newHG):
		sql = "SELECT COUNT() FROM `db`.`chain` WHERE old_ucschg = ? AND new_ucschg = ?"
		return max(row[0] for row in self._db.cursor().execute(sql, (oldHG, newHG)))
	#hasLiftOverChains()
	
	
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
			# if the region overlaps the chain... (1-based, closed intervals)
			if start <= c[2] and end >= c[1]:
				data = chains['data'][chrom][c]
				idx = bisect.bisect(data, (start, sys.maxint, sys.maxint)) - 1
				while (idx < 0) or (data[idx][1] < start):
					idx = idx + 1
				while (idx < len(data)) and (data[idx][0] <= end):
					yield (c[-1], data[idx][0], data[idx][1], data[idx][2], c[4], c[5])
					idx = idx + 1
		#foreach chain
	#_generateApplicableLiftOverChains()
	
	
	def _liftOverRegionUsingChains(self, label, start, end, extra, first_seg, end_seg, total_mapped_sz):
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
		mapped_size = total_mapped_sz - front_diff - (end_seg[2] - end_seg[1] + 1) + end_diff + 1
		
		if mapped_size / float(end - start + 1) >= 0.95: # TODO: configurable threshold?
			mapped_reg = (label, first_seg[5], new_start, new_end, extra)
		
		return mapped_reg
	#_liftOverRegionUsingChains()
	
	
	def generateLiftOverRegions(self, oldHG, newHG, regions, tally=None, errorCallback=None):
		# regions=[ (label,chr,posMin,posMax,extra), ... ]
		oldHG = int(oldHG)
		newHG = int(newHG)
		numNull = numLift = 0
		for region in regions:
			label,chrom,start,end,extra = region
			
			if start > end:
				start,end = end,start
			is_region = (start != end)
			
			# find and apply chains
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
					total_mapped_sz = seg[2] - seg[1] + 1
				elif seg[0] != curr_chain:
					mapped_reg = self._liftOverRegionUsingChains(label, start, end, extra, first_seg, end_seg, total_mapped_sz)
					if mapped_reg:
						break
					curr_chain = seg[0]
					first_seg = seg
					end_seg = seg
					total_mapped_sz = seg[2] - seg[1] + 1
				else:
					end_seg = seg
					total_mapped_sz = total_mapped_sz + seg[2] - seg[1] + 1
			
			if not mapped_reg and first_seg is not None:
				mapped_reg = self._liftOverRegionUsingChains(label, start, end, extra, first_seg, end_seg, total_mapped_sz)
			
			if mapped_reg:
				numLift += 1
				if not is_region:
					mapped_reg = (mapped_reg[0], mapped_reg[1], mapped_reg[2], mapped_reg[2], extra)
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
	
	
	def generateLiftOverLoci(self, oldHG, newHG, loci, tally=None, errorCallback=None):
		# loci=[ (label,chr,pos,extra), ... ]
		regions = ((l[0],l[1],l[2],l[2],l[3]) for l in loci)
		newloci = ((r[0],r[1],r[2],r[4]) for r in self.generateLiftOverRegions(oldHG, newHG, regions, tally, errorCallback))
		return newloci
	#generateLiftOverLoci()
	
	
#Database


# TODO: find a better place for this liftover testing code
"""
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
"""
