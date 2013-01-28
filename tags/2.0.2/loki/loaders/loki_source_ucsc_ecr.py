#!/usr/bin/env python

import itertools
from loki import loki_source


class Source_ucsc_ecr(loki_source.Source):
	"""
	A class to load the pairwise alignments between species as ECRs from the 
	UCSC inter-species alignments
	"""

	
	_remhost = "hgdownload.cse.ucsc.edu"
	_remPath = "goldenPath/hg19/phastCons46way/"
	_comparisons = {"vertebrate" : "", "placentalMammals" : "placental." , "primates" : "primates." }
	_min_sz = 100
	_min_pct = 0.7
	_max_gap = 50
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a2 (2012-09-13)'
	#getVersionString()
	
	
	@classmethod
	def getOptions(cls):
		return {
			'size': 'An integer defining the minimum length of an ECR (default: 100)',
			'identity' : 'A float defining the minimum identity of an ECR (default: 0.7)',
			'gap' : 'An integer defining the maximum gap length below the identity threshold (default: 50)'
		}
	#getOptions()
	
	
	def validateOptions(self, options):
		"""
		Validate the options
		"""
		for o,v in options.iteritems():
			try:
				if o == 'size':
					self._min_sz = int(v)
				elif o == 'identity':
					self._min_pct = float(v)
				elif o == 'gap':
					self._max_gap = int(v)
				else:
					return "unknown option '%s'" % o
			except ValueError:
				return "Cannot parse '%s' parameter value - given '%s'" % (o,v)
		
		return True
	
	def download(self, options):
		"""
		Download the files
		"""
		self._chr_list = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
		file_dict = dict(((sp + ".chr" + ch + ".phastCons.txt.gz", self._remPath + sp + "/chr" + ch + ".phastCons46way." + v + "wigFix.gz") for (sp, v) in self._comparisons.iteritems() for ch in self._chr_list))
		file_dict.update(dict(((sp + ".chrMT.phastCons.txt.gz", self._remPath + sp + "/chrM.phastCons46way." + v + "wigFix.gz") for (sp, v) in self._comparisons.iteritems())))
	
		self.downloadFilesFromFTP(self._remhost,file_dict)

	def update(self, options):
		"""
		Load the data from all of the files
		"""
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# Add a namespace
		ecr_ns = self.addNamespace("ucsc_ecr")
		
		# Add a type of "ecr"
		ecr_typeid = self.addType("ecr")
		
		# Add a type of "ecr_group"
		ecr_group_typeid = self.addType("ecr_group")
		
		# Make sure the '' ldprofile exists
		ecr_ldprofile_id = self.addLDProfile('', 'no LD adjustment')
		
		# Add a containment relationship
		rel_id = self.addRelationship("contains")			
		
		for sp in self._comparisons:
			self.logPush("processing ECRs for " + sp + " ...")
			desc = "ECRs for " + sp
			label = "ecr_" + sp
			
			# Add the group for this species (or comparison)
			ecr_gid = self.addTypedGroups(ecr_group_typeid, [(label, desc)])[0]
			self.addGroupNamespacedNames(ecr_ns, [(ecr_gid, label)])
			
			chr_grp_ids = []
			for ch in self._chr_list + ("MT",):
				ch_id = self._loki.chr_num[ch]
				self.log("processing Chromosome " + ch + " ...")
				f = self.zfile(sp + ".chr" + ch + ".phastCons.txt.gz")
				curr_band = 1
				num_regions = 0
				desc = "ECRs for " + sp + " on Chromosome " + ch
				chr_grp_ids.append(self.addTypedGroups(ecr_group_typeid, [("ecr_%s_chr%s" % (sp, ch), desc)])[0])
				self.addGroupNamespacedNames(ecr_ns, [(chr_grp_ids[-1], "ecr_%s_chr%s" % (sp, ch))])
				band_grps = []
				grp_rid = {}
				for regions in self.getRegions(f):
					label = "ecr_%s_chr%s_band%d" % (sp, ch, curr_band)
					desc = "ECRs for " + sp + " on Chromosome " + ch + ", Band %d" % (curr_band,)
					num_regions += len(regions)
					
					if regions:
						band_grps.append((label, desc))
																
					# Add the region itself
					reg_ids = self.addTypedBiopolymers(ecr_typeid, ((self.getRegionName(sp, ch, r), '') for r in regions))
					# Add the name of the region
					self.addBiopolymerNamespacedNames(ecr_ns, zip(reg_ids, (self.getRegionName(sp, ch, r) for r in regions)))
					# Add the region Boundaries
					# This gives a generator that yields [(region_id, (chrom_id, start, stop)) ... ]
					region_bound_gen = zip(((i,) for i in reg_ids), ((ch_id, r[0], r[1]) for r in regions))
					self.addBiopolymerLDProfileRegions(ecr_ldprofile_id, (tuple(itertools.chain(*c)) for c in region_bound_gen))			
	
					if regions:
						grp_rid[band_grps[-1]] = reg_ids
						#Add the region to the group
						#self.addGroupBiopolymers(((band_gids[-1], r_id) for r_id in reg_ids))
						
					curr_band += 1

				
				band_gids = self.addTypedGroups(ecr_group_typeid, band_grps)
				self.addGroupNamespacedNames(ecr_ns, zip(band_gids, (r[0] for r in band_grps)))
				gid_rid = []
				for i in xrange(len(band_gids)):
					gid_rid.extend(((band_gids[i], rid) for rid in grp_rid[band_grps[i]]))
				
				self.addGroupBiopolymers(gid_rid)

				self.addGroupRelationships(((chr_grp_ids[-1], b, rel_id, 1) for b in band_gids))
				
				self.log("OK (%d regions found in %d bands)\n" % (num_regions, curr_band - 1))
			
			self.addGroupRelationships(((ecr_gid, c, rel_id, 1) for c in chr_grp_ids))
			
			self.logPop("... OK\n")

				
				
	def getRegionName(self, species, ch, region):
		"""
		Returns a string representation of the name
		"""
		return species + ":chr" + ch + ":" + str(region[0]) + "-" + str(region[1])
				
	
	def getRegions(self, f):
		"""
		Yields the regions that meets the thresholds with a given maximum gap
		"""
		running_sum = 0
		n_pos = 0
		curr_gap = 0
		curr_pos = 1
		curr_start = 1
		curr_end = 0
		step = 1
		
		line = f.next()
		curr_band = []
		
		for l in f:
			try:
				p = float(l)
				if p >= self._min_pct:
					#If this is the 1st time we crossed the threshold, start the counters
					if curr_gap != 0 and running_sum / float(n_pos) < self._min_pct:
						if curr_end- curr_start >= self._min_sz:
							curr_band.append((curr_start, curr_end))
						
						# Restart the region tracking
						running_sum = 0
						n_pos = 0
					
					if n_pos == 0:
						curr_start = curr_pos
											
					curr_end = curr_pos
					running_sum += p
					n_pos += 1
					curr_gap = 0
				# If this is true, we're searching a gap				
				elif n_pos != 0:
					#print curr_gap, curr_end - curr_start
					# If we are in an acceptable gap, don't add on to the end
					if curr_gap < self._max_gap:
						running_sum += p
						n_pos += 1
						curr_gap += 1
					else:
						# If it's big enough, add it (we ran off the end of the gap)
						if curr_end - curr_start > self._min_sz:
							curr_band.append((curr_start, curr_end))
						n_pos = 0
						curr_start = 0
						curr_end = 0
						running_sum = 0
						curr_gap = 0
				# otherwise, just keep on trucking		
				curr_pos += step
			except ValueError:
				# At this point, we have a format line
				d = dict((v.split('=',2) for v in l.split() if v.find('=') != -1))
				
				# If this is moving us to a different place, we have to restart,
				# o/w just keep on trucking
				if int(d['start']) != curr_pos or int(d['step']) != step:
					
					if curr_end - curr_start > self._min_sz and running_sum / float(n_pos) >= self._min_pct and curr_gap < self._max_gap:
						curr_band.append((curr_start, curr_end))
					
					yield curr_band
					curr_band = []	
	
					running_sum = 0
					n_pos = 0
					curr_gap = 0
					curr_pos = int(d['start'])
					curr_start = int(d['start'])
					step = int(d['step'])
					curr_end = 0
				
		
		# Check on the last region...
		if curr_end - curr_start > self._min_sz and running_sum / float(n_pos) >= self._min_pct and curr_gap < self._max_gap:
			 curr_band.append((curr_start, curr_end))
		
		yield curr_band
