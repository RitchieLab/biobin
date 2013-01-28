#!/usr/bin/env python

import itertools
import os
import re
from loki import loki_source


class Source_chainfiles(loki_source.Source):
	"""
	A loader that loads all of the chainfiles into LOKI
	"""
	
	
	##################################################
	# private class data
	
	
	_reDir = re.compile('^hg[0-9]+$', re.IGNORECASE)
	_reFile = re.compile('^hg([0-9]+)tohg([0-9]+)\.over\.chain\.gz$', re.IGNORECASE)
	
	
	##################################################
	# source interface
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-08-30)'
	#getVersionString()
	
	
	def download(self, options):
		# define a callback to search for all available hgX liftover chain files
		def remFilesCallback(ftp):
			remFiles = {}
			ftp.cwd('/goldenPath')
			dirs = [ d for d in ftp.nlst() if self._reDir.match(d) ]
			for d in dirs:
				ftp.cwd('/goldenPath/%s/liftOver' % d)
				files = [ f for f in ftp.nlst() if self._reFile.match(f) ]
				for f in files:
					remFiles[f] = '/goldenPath/%s/liftOver/%s' % (d,f)
			
			return remFiles
		#remFilesCallback
		
		self.downloadFilesFromFTP("hgdownload.cse.ucsc.edu", remFilesCallback)
	#download()
	
	
	def update(self, options):
		"""
		Parse all of the chain files and insert them into the database
		"""	
		
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		for fn in os.listdir('.'):
			match = self._reFile.match(fn)
			if not match:
				continue
			old_ucschg = int(match.group(1))
			new_ucschg = int(match.group(2))
			self.log("parsing chains for hg%d -> hg%d ..." % (old_ucschg,new_ucschg))
			f = self.zfile(fn)
			
			is_hdr = True
			is_valid = True
			chain_hdrs = []
			chain_data = []
			curr_data = []
			for line in f:
				if is_hdr:
					if line:
						try:
							chain_hdrs.append(self._parseChain(line))
						except:
							is_valid = False
						is_hdr = False
				elif line:
					if is_valid:
						curr_data.append(line)
				else:
					if is_valid:
						chain_data.append(self._parseData(chain_hdrs[-1], '\n'.join(curr_data)))
					is_valid = True
					curr_data = []
					is_hdr = True
			
			hdr_ids = self.addChains(old_ucschg, new_ucschg, chain_hdrs)
			
			# Now, I want to take my list of IDs and my list of list of 
			# tuples and convert them into a list of tuples suitable for
			# entering in the chain_data table
			chain_id_data = zip(hdr_ids, chain_data)
			chain_data_itr = (tuple(itertools.chain((chn[0],),seg)) for chn in chain_id_data for seg in chn[1])
			
			self.addChainData(chain_data_itr)
			
			self.log("OK\n")
		# for fn in dir
		
	#update()
	
	def _parseChain(self, chain_hdr):
		"""
		Parses the chain header to extract the information required 
		for insertion into the database
		"""
		
		# get the 1st line
		hdr = chain_hdr.strip().split('\n')[0].strip()
		# Parse the first line
		wds = hdr.split()
		
		if wds[0] != "chain":
			raise Exception("Not a valid chain file")
		
		if wds[2][3:] not in self._loki.chr_num:
			raise Exception("Could not find chromosome: " + wds[2][3:] + "->" + wds[7][3:])
			
		is_fwd = (wds[9] == "+")
		if is_fwd:
			new_start = int(wds[10])
			new_end = int(wds[11])
		else:
			# NOTE: If we're going backward, this will mean that 
			# end < start
			new_start = int(wds[8]) - int(wds[10])
			new_end = int(wds[8]) - int(wds[11])
		
		
		# I want a tuple of (score, old_chr, old_start, old_end,
		# new_chr, new_start, new_end, is_forward)
		return (int(wds[1]), 
			self._loki.chr_num[wds[2][3:]], int(wds[5]), int(wds[6]),
			self._loki.chr_num.get(wds[7][3:],-1), new_start, new_end,
			int(is_fwd))
		
	def _parseData(self, chain_tuple, chain_data):
		"""
		Parses the chain data into a more readily usable and iterable 
		form (the data of the chain is everything after the 1st line)
		"""
		_data = [ tuple([int(v) for v in l.split()]) for l in chain_data.split('\n')[:-1] ]
		
		curr_pos = chain_tuple[2]
		new_pos = chain_tuple[5]
			
		_data_txform = []
		for l in _data:
			_data_txform.append((curr_pos, curr_pos + l[0], new_pos))
			curr_pos = curr_pos + l[0] + l[1]
			if chain_tuple[7]:
				new_pos = new_pos + l[0] + l[2]
			else:
				new_pos = new_pos - l[0] - l[2]
			
		_data_txform.append((curr_pos, curr_pos + int(chain_data.split()[-1]), new_pos))
		
		return _data_txform
	
#class Source_chainfiles
	
