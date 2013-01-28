#!/usr/bin/env python

import argparse
import os
import posixpath
import shutil
import sys
import tarfile
import tempfile


# Add one level up to the beginning of the python path to allow for running from non-installed location
if __name__ == "__main__":
	sys.path.insert(0,os.path.abspath(os.path.join(__file__, "..", "..")))

from loki import loki_db

if __name__ == "__main__":
	version = "LOKI version %s" % (loki_db.Database.getVersionString())
	
	# define arguments
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=version,
	)
	parser.add_argument('--version', action='version',
			version=version+"\n%s version %s\n%s version %s" % (
				loki_db.Database.getDatabaseDriverName(), loki_db.Database.getDatabaseDriverVersion(),
				loki_db.Database.getDatabaseInterfaceName(), loki_db.Database.getDatabaseInterfaceVersion()
			)
	)
	parser.add_argument('-k', '--knowledge', type=str, metavar='file', action='store', default=None,
			help="the knowledge database file to use"
	)
	parser.add_argument('-a', '--archive', type=str, metavar='file', action='store', default=None,
			help="create (or re-use and update) a compressed archive of downloaded source data files"
	)
	parser.add_argument('--from-archive', type=str, metavar='file', action='store', default=None,
			help="an input source data archive to re-use but not update"
	)
	parser.add_argument('--to-archive', type=str, metavar='file', action='store', default=None,
			help="an output source data archive to create (or replace) but not re-use"
	)
	parser.add_argument('-d', '--temp-directory', type=str, metavar='dir', action='store', default=None,
			help="a directory to use for temporary storage of downloaded or archived source data files (default: platform dependent)"
	)
	parser.add_argument('-l', '--list-sources', type=str, metavar='source', nargs='*', action='append', default=None,
			help="list versions and options for the specified source loaders, or if none or '+' are specified, list all available sources"
	)
	parser.add_argument('-c', '--cache-only', action='store_true',
			help="do not download any new source data files, only use what's available in the provided archive"
	)
	parser.add_argument('-u', '--update', type=str, metavar='source', nargs='*', action='append', default=None,
			help="update the knowledge database file by downloading and processing new data from the specified sources, "
			+"or if none or '+' are specified, from all available sources"
	)
	parser.add_argument('-U', '--update-except', type=str, metavar='source', nargs='*', action='append', default=None,
			help="update the knowledge database file by downloading and processing new data from all available sources EXCEPT those specified"
	)
	parser.add_argument('-o', '--option', type=str, metavar=('source','optionstring'), nargs=2, action='append', default=None,
			help="additional option(s) to pass to the specified source loader module, in the format 'option=value[,option2=value2[,...]]'"
	)
	parser.add_argument('-f', '--finalize', action='store_true',
			help="finalize the knowledge database file"
	)
	parser.add_argument('-v', '--verbose', action='store_true',
			help="print warnings and log messages"
	)
	parser.add_argument('-t', '--test-data', action='store_true',
			help="Load testing data only"
	)
	
	# if no arguments, print usage and exit
	if len(sys.argv) < 2:
		print version
		print
		parser.print_usage()
		print
		print "Use -h for details."
		sys.exit(2)
	
	# parse arguments
	args = parser.parse_args()
		
	# instantiate database and load knowledge file, if any
	db = loki_db.Database(testing=args.test_data, updating=((args.update != None) or (args.update_except != None)))
	db.setVerbose(args.verbose)
	db.attachDatabaseFile(args.knowledge)
	
	# list sources?
	if args.list_sources != None:
		srcSet = set()
		for srcList in args.list_sources:
			srcSet |= set(srcList)
		if (not srcSet) or ('+' in srcSet):
			print "available source loaders:"
			srcSet = set()
		else:
			print "source loader options:"
		moduleVersions = db.getSourceModuleVersions(srcSet)
		moduleOptions = db.getSourceModuleOptions(srcSet)
		for srcName in sorted(moduleOptions.keys()):
			print "  %s : %s" % (srcName,moduleVersions[srcName])
			if moduleOptions[srcName]:
				for srcOption in sorted(moduleOptions[srcName].keys()):
					print "    %s = %s" % (srcOption,moduleOptions[srcName][srcOption])
			elif srcSet:
				print "    <no options>"
	
	# pass options?
	userOptions = {}
	if args.option != None:
		for optList in args.option:
			srcName = optList[0]
			if srcName not in userOptions:
				userOptions[srcName] = {}
			for optString in optList[1].split(','):
				opt,val = optString.split('=',1)
				userOptions[srcName][opt] = val
	userOptions = userOptions or None
	
	# parse requested update sources
	srcSet = None
	if args.update != None:
		srcSet = set()
		for srcList in args.update:
			srcSet |= set(srcList)
	notSet = None
	if args.update_except != None:
		notSet = set()
		for srcList in args.update_except:
			notSet |= set(srcList)
	
	# update?
	if (srcSet != None) or (notSet != None):
		if srcSet and '+' in srcSet:
			srcSet = set()
		srcSet = (srcSet or set(db.getSourceModules())) - (notSet or set())
		
		# create temp directory and unpack input archive, if any
		startDir = os.getcwd()
		fromArchive = args.from_archive or args.archive
		toArchive = args.to_archive or args.archive
		cacheDir = os.path.abspath(tempfile.mkdtemp(prefix='loki_update_cache.', dir=args.temp_directory))
		if args.temp_directory:
			print "using temporary directory '%s'" % cacheDir
		
		# try/finally to make sure we clean up the tempdir at the end
		try:
			if fromArchive:
				if os.path.exists(fromArchive) and tarfile.is_tarfile(fromArchive):
					print "unpacking archived source data files from '%s' ..." % fromArchive
					with tarfile.open(name=fromArchive, mode='r:*') as archive:
						# the archive should only contain directories named after sources,
						# so we can filter members by their normalized top-level directory
						for member in archive:
							srcName = posixpath.normpath(member.name).split('/',1)[0]
							if (not srcName) or srcName.startswith('.'):
								continue
							# if we're not writing an output archive, we only have to extract
							# the directories for the sources we need
							if (not toArchive) and (srcName not in srcSet):
								continue
							archive.extractall(cacheDir, [member])
					#with archive
					print "... OK"
				else:
					print "source data archive '%s' not found, starting fresh" % fromArchive
			#if fromArchive
			
			os.chdir(cacheDir)
			db.updateDatabase(srcSet, userOptions, args.cache_only)
			os.chdir(startDir)
			
			# create output archive, if requested
			if toArchive:
				print "archiving source data files in '%s' ..." % toArchive
				with tarfile.open(name=toArchive, mode='w:gz') as archive:
					for filename in sorted(os.listdir(cacheDir)):
						archive.add(os.path.join(cacheDir, filename), arcname=filename)
				print "... OK"
		finally:
			# clean up temp directory
			def onerror(func, path, exc):
				print "WARNING: unable to remove temporary file '%s': %s\n" % (path,exc)
			shutil.rmtree(cacheDir, onerror=onerror)
	#update

	# finalize?
	if args.finalize:
		db.finalizeDatabase()
	
#__main__
