import loki_source

class Source_test_gene (loki_source.Source):
	
	def download(self, options):
		"""
		No download needed
		"""
		pass
		
	def update(self, options):
		"""
		Create some genes (and a test population)
		"""
		
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		namespaceID = self.addNamespaces([
			('symbol',      0),
			('entrez_gid',  0),
			('alternate_id', 0),
			('multi_ids', 1)
		])
		
		typeID = self.addTypes([
			('gene',),
		])
		
		# Set the zone size to 9
		self._loki.setDatabaseSetting("zone_size", "9")
		
		self.log("Creating LD Profiles ...")
		ld_ids = self.addLDProfiles([('', 'no LD adjustment', None), ('PLUS5', 'Add 5 to all boundaries', 'Adding 5')])
		self.log("OK\n")
		
		self.log("Creating Regions ...")
		b_ids = self.addTypedBiopolymers(typeID['gene'], (("G%d" % (r+1), "Gene %d" % (r+1)) for r in range(6)))
		self.log("OK\n")
		
		self.log("Creating Region Aliases ...")
		self.addBiopolymerNamespacedNames(namespaceID['symbol'], ((b_ids[r], "G%d" % (r+1)) for r in range(6)))
		self.addBiopolymerNamespacedNames(namespaceID['entrez_gid'], ((b_ids[r], "%d" % (r+1)) for r in range(6)))
		self.addBiopolymerNamespacedNames(namespaceID['alternate_id'],((b_ids[r], "R%d" % (r+1)) for r in range(6)))
		self.addBiopolymerNamespacedNames(namespaceID['multi_ids'], [(b_ids[1], "G23"), (b_ids[2], "G23")])
		self.log("OK\n")		
		
		self.log("Creating Region Boundaries ...")
		gene_bounds = dict((r+1, (20*r+5, 20*r+15)) for r in range(4))
		gene_bounds[5] = (30,50)
		gene_bounds[6] = (10, 20)
		
		self.addBiopolymerLDProfileRegions(ld_ids[''],((b_ids[r],1+(r==5),gene_bounds[r+1][0],gene_bounds[r+1][1])for r in range(6)))
		self.addBiopolymerLDProfileRegions(ld_ids['PLUS5'],((b_ids[r],1+(r==5),gene_bounds[r+1][0]+5,gene_bounds[r+1][1]+5)for r in range(6)))
		self.log("OK\n")