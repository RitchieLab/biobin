#!/usr/bin/env python

import sys
import bisect

class liftOver(object):
	"""
	A class to do your heavy lifting for you
	"""
	
	def __init__(self, db, old_ucschg, new_ucschg, cached=False):
		# db is a loki_db.Database object
		self._db = db
		self._old_ucschg = old_ucschg
		self._new_ucschg = new_ucschg
		self._cached = cached
		self._minFrac = 0.95
		if self._cached:
			self._cached_data = {}
			self._cached_keys = {}
			self._chainData = self._initChains()
	
	def _initChains(self):
		"""
		Constructs the chain object that resides in memory
		"""
		for row in self._db._db.cursor().execute("SELECT chain_id, old_chr, score, chain.old_start, " + 
			"chain.old_end, chain.new_start, is_fwd, new_chr, " + 
			"chain_data.old_start, chain_data.old_end, chain_data.new_start " + 
			"FROM db.chain INNER JOIN db.chain_data USING (chain_id) " +
			"WHERE old_ucschg=? AND new_ucschg=?" + 
			"ORDER BY old_chr, score DESC, chain_data.old_start",
			(self._old_ucschg,self._new_ucschg)):
				
			chain = (row[2], row[3], row[4], row[5], row[6], row[7], row[0])
			chr = row[1]
					
			if chr not in self._cached_data:
				self._cached_data[chr] = {chain: []}
				self._cached_keys[chr] = [chain]
			elif chain not in self._cached_data[chr]:
				self._cached_data[chr][chain] = []
				self._cached_keys[chr].append(chain)
			
			self._cached_data[chr][chain].append((row[8],row[9],row[10]))
		
		# Sort the chains by score
		for k in self._cached_keys:
			self._cached_keys[k].sort(reverse=True)

			
				
	def _findChains(self, chrom, start, end):
		"""
		Finds all of the chain segments, either from the  
		"""
		
		if not self._cached:
			for row in self._db._db.cursor().execute(
			"SELECT chain.chain_id, chain_data.old_start, chain_data.old_end, chain_data.new_start, is_fwd, new_chr " +
			"FROM chain INNER JOIN chain_data ON chain.chain_id = chain_data.chain_id " +
			"WHERE old_ucschg=? AND new_ucschg=? AND old_chr=? AND chain.old_end>=? AND chain.old_start<? AND chain_data.old_end>=? AND chain_data.old_start<? " +
			"ORDER BY score DESC",
			(self._old_ucschg, self._new_ucschg, chrom, start, end, start, end)):
				yield row
		else:
			for c in self._cached_keys.get(chrom, []):
				# if the region overlaps the chain...
				if start <= c[2] and end >= c[1]:
					data = self._cached_data[chrom][c]
					idx = bisect.bisect(data, (start, sys.maxint, sys.maxint))
					if idx:
						idx = idx-1

					if idx < len(data) - 1 and start == data[idx + 1]:
						idx = idx + 1
					
					while idx < len(data) and data[idx][0] < end:
						yield (c[-1], data[idx][0], data[idx][1], data[idx][2], c[4], c[5])
						idx = idx + 1
					
					
	def liftRegion(self, chrom, start, end):
		"""
		Lift a region from a given chromosome on a specified assembly
		to a region on the new assembly
		"""
		
		# We need to actually lift regions to detect dropped sections
		is_region = True
		
		# If the start and end are swapped, reverse them, please
		if start > end:
			(start, end) = (end, start)
		elif start == end:
			is_region = False
			end = start + 1	
		
		ch_list = self._findChains(chrom, start, end)
		
		# This will be a tuple of (start, end) of the mapped region
		# If the function returns "None", then it was unable to map
		# the region into the new assembly
		mapped_reg = None
		
		curr_chain = None
		
		total_mapped_sz = 0
		first_seg = None
		end_seg = None
		for seg in ch_list:
			if curr_chain is None:
				curr_chain = seg[0]
				first_seg = seg
				end_seg = seg
				total_mapped_sz = seg[2] - seg[1]
			elif seg[0] != curr_chain:
				mapped_reg = self._mapRegion((start, end), first_seg, end_seg, total_mapped_sz)
				if not mapped_reg:
					first_seg = seg
					end_seg = seg
					total_mapped_sz = seg[2] - seg[1]
				else:
					break
			else:
				end_seg = seg
				total_mapped_sz = total_mapped_sz + seg[2] - seg[1]
				
		if not mapped_reg and first_seg is not None:
			mapped_reg = self._mapRegion((start, end), first_seg, end_seg, total_mapped_sz)
		
		if mapped_reg and not is_region:
			mapped_reg = (mapped_reg[0], mapped_reg[1], mapped_reg[1]) #bug?
		
		return mapped_reg
		
				
				

	def _mapRegion(self, region, first_seg, end_seg, total_mapped_sz):
		"""
		Map a region given the 1st and last segment as well as the total mapped size
		"""
		mapped_reg = None
		
		# The front and end differences are the distances from the
		# beginning of the segment.
		
		# The front difference should be >= 0 and <= size of 1st segment
		front_diff = max(0, min(region[0] - first_seg[1], first_seg[2] - first_seg[1]))
		
		# The end difference should be similar, but w/ last
		end_diff = max(0, min(region[1] - end_seg[1], end_seg[2] - end_seg[1]))
		
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
		
		if mapped_size / float(region[1] - region[0]) >= self._minFrac:
			mapped_reg = (first_seg[5], new_start, new_end)
			
		return mapped_reg

if __name__ == "__main__":
	import loki_db
	db = loki_db.Database(sys.argv[2])
	
	old = int(sys.argv[5]) if (len(sys.argv) > 5) else 18
	new = int(sys.argv[6]) if (len(sys.argv) > 6) else 19
	lo = liftOver(db, old, new, False)
	f = file(sys.argv[1])
	m = file(sys.argv[3],'w')
	u = file(sys.argv[4],'w')
	
	for l in f:
		wds = l.split()
		chrm = lo._db.chr_num.get(wds[0][3:],-1)
		n = lo.liftRegion(chrm, int(wds[1]), int(wds[2]))
		if n:
			print >> m, "chr" + lo._db.chr_list[n[0]-1], n[1], n[2]
		else:
			print >> u, l
