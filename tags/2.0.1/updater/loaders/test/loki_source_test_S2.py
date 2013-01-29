import loki_source

class Source_test_S2 (loki_source.Source):
	
	def download(self,options):
		"""
		No downloading
		"""
		pass
		
	def update(self,options):
		"""
		Create source 2
		"""
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		path_typeid = self.addType("ontology")
		gene_typeid = self.addType("gene")
		
		relationship_id = self.addRelationship("regulates")
		
		ns_id = self.addNamespace("S2")
		symbol_ns = self.addNamespace("symbol")
		
		self.log("Creating Groups ...")
		gids = self.addTypedGroups(path_typeid, [("P3", "Pathway 3"), ("P4", "Pathway 4 (==P1)")])		
		self.addGroupNamespacedNames(ns_id, ((gids[r], "P%d" % (r+3)) for r in range(2)))
		self.log("OK\n")
		
		self.log("Creating Group Relationships ...")
		self.addGroupRelationships([(gids[0],gids[1],relationship_id)])		
		self.log("OK\n")
		
		self.log("Creating Group Membership ...")
		self.addGroupMemberTypedNamespacedNames(gene_typeid, symbol_ns, 
			[(gids[1],1,'G1'),(gids[1],2,'G2'),(gids[0],1,'G2'),(gids[0],2,'G3')])
		self.log("OK\n")