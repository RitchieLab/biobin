import loki_source

class Source_test_snp(loki_source.Source):
	"""
	Class to load testing dataset
	"""
	
	def download(self, options):
		"""
		Do nothing to download
		"""
		pass
		
	def update(self, options):
		"""
		Create the testing dataset
		"""
		
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# Add some roles
		self.log("Creating Roles ...")
		role_list = [
			("Exon", "Test Exon", 1, 1), 
			("Regulatory", "Test Regulatory", 0, 1), 
			("Intron", "Test Intron", 0, 0)]
		role_ids = self.addRoles(role_list)
			
		self.log("OK\n")
			
		chr_snp = {1 : [r+1 for r in range(11)], 
				   2 : [r+21 for r in range(4)]}
		
		snp_gene = { 1: [1], 2: [1], 4:[2], 5: [2, 5], 6: [5], 7: [3,5], 8: [3],
			10: [4], 11: [4], 22: [6], 23: [6]}
			
		# Add some snp merges: 1XX -> XX on chr 1 only
		self.log("Creating SNP Merge Records ...")
		self.addSNPMerges(zip((r+100 for r in chr_snp[1]), chr_snp[1]))
		self.log("OK\n")
		
		# Add some role data
		self.log("Creating SNP Roles ...")
		for snp, v in snp_gene.iteritems():
			self.addSNPEntrezRoles(((snp, g, role_ids[role_list[(snp+g)%3][0]]) for g in v))
		self.log("OK\n")
		
		# Add some SNPs
		self.log("Creating SNPs ...")
		for chr in chr_snp:
			self.addChromosomeSNPLoci(chr, ((r, r*7 - (chr-1)*20*7, r%4!=0) for r in chr_snp[chr]))
		self.log("OK\n")