#!/usr/bin/env python

import subprocess
import os
import sys
import itertools
import gzip
import re
import shutil
from ftplib import FTP
from optparse import OptionParser

def downloadFile(ftp_obj, fn):
	"""
	Downloads and checks the file, retrying as necessary
	"""
	print "Downloading", fn
	# TODO: retry on errors
	n_tries = 0
	max_tries = 10
	
	wd = ftp_obj.pwd()	
	
	fs_re = re.compile(r"(?:\S*\s*){4,4}(\d*)")
	retry = True
	
	while retry and n_tries < max_tries:
		n_tries = n_tries + 1
		ls_result = []
		try:
			f = file(fn, 'wb')
			ftp_obj.retrbinary("RETR " + fn, lambda x: f.write(x))
			f.close()
			ftp_obj.dir(fn, lambda x: ls_result.append(x))
			# check to make sure we got the whole file
			ftp_sz = int(fs_re.match(ls_result[0]).groups()[0])
			f_inf = os.stat(fn)
			# If the file sizes match, then we probably got it right.
			# I'm not going to worry about solar flares flipping a bit!
			if ftp_sz == f_inf.st_size:
				retry = False
			
		except ftplib.all_errors, e:
			f.close()
			ftp_obj.connect(ftp_obj.host)
			ftp_obj.login()
			ftp_obj.cwd(wd)
	
	if n_tries == max_tries:
		raise Exception("Could not download " + fn)


def downloadFiles(population):
	"""
	Downloads and unzips the HapMap files from the given population
	(population must be a 3-letter code used by HapMap)
	"""
	f_srv = FTP("ftp.ncbi.nlm.nih.gov")
	f_srv.login()
	f_srv.cwd("hapmap/ld_data/latest")
	
	# TODO: make this a little more restrictive, please
	file_list = f_srv.nlst("*" + population + "*")
	for fn in file_list:
		downloadFile(f_srv, fn)
	
	f_srv.close()
	
	fn_list = []
	# OK, now we have all the files, unzip them
	for fn in file_list:
		if False:
			print "Extracting", fn
			fz = gzip.open(fn)
			f = file(os.path.splitext(fn)[0],'w')
			for l in fz:
				print >> f, l,
			f.close()
			fz.close()
			os.remove(fz.name)
			fn_list.append(f.name)
		else:
			fn_list.append(fn)
	
	# And return a list of all the files that we've downloaded
	return fn_list		


def genPops(popList, dprimes, rsquared, opts):
	"""
	Generates the populations using LD-spline
	"""
	valid_pops = {"ASW" : "ASW", "A" : "ASW",
	              "CEU" : "CEU", "C" : "CEU",
	              "CHB" : "CHB", "H" : "CHB",
	              "CHD" : "CHD", "D" : "CHD",
	              "GIH" : "GIH", "G" : "GIH",
	              "JPT" : "JPT", "J" : "JPT",
	              "LWK" : "LWK", "L" : "LWK",
	              "MEX" : "MEX", "M" : "MEX",
	              "MKK" : "MKK", "K" : "MKK",
	              "TSI" : "TSI", "T" : "TSI",
	              "YRI" : "YRI", "Y" : "YRI"}
	              
	datadir = ".__data"
	              
	if not os.path.isdir(datadir):
		# If "data" exists and is not a directory, we'll fail here
		os.mkdir(datadir)
		
	os.chdir(datadir)
	
	# download the chain file
	f_srv = FTP("hgdownload.cse.ucsc.edu")
	f_srv.login("anonymous")
	f_srv.cwd("goldenPath/hg18/liftOver")
	chainfn = "hg18ToHg19.over.chain.gz"
	downloadFile(f_srv, chainfn)
	f_srv.close()
	

	
	# Now, write the configuration
	cfg_f = file("ldspline.cfg", "w")
	print >> cfg_f, "rs", " ".join((str(v) for v in rsquared))
	print >> cfg_f, "dp", " ".join((str(v) for v in dprimes))
	
	for pop in popList:
		# Find the population extension
		pop_ext = valid_pops[pop]	
		
		fn_list = downloadFiles(pop_ext)
		
		idx_file = file(pop_ext, 'w')
		
		# Pattern to extract the chromosome from the filename
		fn_pattern = re.compile("ld_chr(.*)_.*")
		for fn in fn_list:
			chr_list = fn_pattern.findall(fn)
			if chr_list:
				print >> idx_file, chr_list[0], "\t", fn
		
		idx_file.close()			
		
		# Now, get the splines
		print "Building LD Splines"
		os.system(opts.ldspline + " load " + idx_file.name)
		if not opts.keep:
			for fn in fn_list:
				os.remove(fn)
		
		# Export the splines
		print "Lifting splines to build 37"
		os.system(opts.ldspline + " export-lomap " + pop_ext)
		
		# Lift over the splines (36->37)
		os.system(opts.liftover + " " + pop_ext + ".bim " + chainfn + " " + pop_ext + ".new " + pop_ext + ".unmapped")
		
		# Re-import the splines
		os.system(opts.ldspline + " import-lomap " + pop_ext + " 37")
		if not opts.keep:
			os.remove(pop_ext + ".ldspline")
		
		# print this population's LDspline location 
		print >> cfg_f, pop_ext, os.path.join(os.getcwd(), pop_ext+"-b37.ldspline"), pop_ext + " population from HapMap"

	cfg_f.close()
	
	# OK, we're done w/ all the populations, time to give it to biofilter!
	print "Creating populations in biofilter"
	os.chdir("..")
	os.system(opts.poploader + " --DB " + opts.db + " --config " + os.path.join(datadir,cfg_f.name))
	
	#.... and we're done, so clean up after yourself, please!
	if not opts.keep:
		shutil.rmtree(datadir)

if __name__ == "__main__":
	parser = OptionParser()
	
	parser.add_option("-r", "--rsquared", dest="rsq", action="append",
		help="Comma-separated list of R-squared value cutoffs")
	parser.add_option("-p", "--populations", dest="pops", action="append",
		help="Comma-separated list of HapMap populations to add")
	parser.add_option("-d", "--dprime", dest="dprime", action="append",
		help="Comma-separated list of D-prime cutoffs")
	parser.add_option("-l", "--liftover", dest="liftover", action="store",
		help="Location of the liftOver executable (default is in path)",
		default="liftOver")
	parser.add_option("-s", "--ldspline", dest="ldspline", action="store",
		help="Location of the ldspline executable (default is in path)",
		default="ldspline")
	parser.add_option("-o", "--poploader", dest="poploader", action="store",
		help="Location of the pop_loader executable (default is in path)",
		default="pop_loader")
	parser.add_option("-b", "--db", dest="db", action="store",
		help="Location of the LOKI databse (default 'knowledge.bio')",
		default="knowledge.bio")
	parser.add_option("-k", "--keep-data", dest="keep", action="store_true",
		help="Do not delete data after finished working with it",
		default=False)
		
	(opts, args) = parser.parse_args();
	rsq = []
	try:
		if opts.rsq is not None:
			rsq = [float(f) for f in itertools.chain(*(a.split(",") for a in opts.rsq))]
	except ValueError, e:
		print "Error: Invalid R-Squared value"
		sys.exit(1)
		
	dprime = []
	try:
		if opts.dprime is not None:
			dprime = [float(f) for f in itertools.chain(*(a.split(",") for a in opts.dprime))]
	except ValueError, e:
		print "Error: Invalid R-Squared value"	
		sys.exit(2)
				
	if not rsq and not dprime:
		print "Error: You must provide R-squared or D-prime cutoffs"
		sys.exit(3)
	
	pops = []
	if opts.pops is not None:
		pops = [s.upper() for s in itertools.chain(*(a.split(",") for a in opts.pops))]
		
	if not pops:
		print "Error: You must supply at least one population"
		sys.exit(5)	
	
	# Try to find all of the child programs needed now
	exe_error = False
	null_f = file(os.devnull, 'w')
	try:
		if os.path.split(opts.ldspline)[0]:
			opts.ldspline = os.path.abspath(opts.ldspline)
		subprocess.call(opts.ldspline, stdout=null_f, stderr=subprocess.STDOUT)
	except OSError, e:
		print "Error: could not find ldspline executable"
		exe_error = True
	try:
		if os.path.split(opts.poploader)[0]:
			opts.poploader = os.path.abspath(opts.poploader)
		subprocess.call(opts.poploader, stdout=null_f, stderr=subprocess.STDOUT)
	except OSError, e:
		print "Error: could not find pop_loader executable"
		exe_error = True
	try:
		if os.path.split(opts.liftover)[0]:
			opts.liftover = os.path.abspath(opts.liftover)
		subprocess.call(opts.liftover, stdout=null_f, stderr=subprocess.STDOUT)
	except OSError, e:
		print "Error: could not find liftOver executable"
		exe_error = True
	
	null_f.close()	
	
	if exe_error:
		sys.exit(6)
		
	try:
		genPops(pops, dprime, rsq, opts)
	except KeyError, e:
		print "Error: Invalid population"
		sys.exit(4)
	except Exception, e:
		print e
		raise

