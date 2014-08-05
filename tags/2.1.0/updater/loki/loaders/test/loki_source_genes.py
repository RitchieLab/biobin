#!/usr/bin/env python

from loki import loki_source


class Source_genes(loki_source.Source):
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0 (2013-02-14)'
	#getVersionString()
	
	
	def download(self, options):
		pass
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# get or create the required metadata records
		ldprofileID = self.addLDProfiles([
			('',   'no LD adjustment',   None, None),
			('ld', 'some LD adjustment', None, None),
		])
		namespaceID = self.addNamespaces([
			('gene',       0),
			('entrez_gid', 0), # needed to resolve snp_entrez_role
			('protein',    1),
		])
		typeID = self.addTypes([
			('gene',),
		])
		
		# define genes
		self.log("adding genes to the database ...")
		listGene = [
			#(label,desc)
			('A', 'normal gene'),
			('B', 'normal gene'),
			('C', 'overlapping gene'),
			('D', 'overlapping gene'),
			('E', 'gene with 2 regions'),
			('F', 'gene with no SNPs'),
			('G', 'gene with no regions'),
			('H', 'overlapping gene'),
			('I', 'overlapping gene'),
			('P', 'gene with only nearby SNPs'),
			('Q', 'normal gene'),
			('R', 'normal gene'),
			('S', 'gene with 2 regions'),
		]
		listBID = self.addTypedBiopolymers(typeID['gene'], listGene)
		geneBID = dict(zip((g[0] for g in listGene), listBID))
		self.log(" OK: %d genes\n" % len(geneBID))
		
		# define gene aliases
		self.log("adding gene identifiers to the database ...")
		genEName = ((bid,ord(g)-64) for g,bid in geneBID.iteritems()) # A->1, B->2, ... S->19 ... Z->26
		self.addBiopolymerNamespacedNames(namespaceID['entrez_gid'], genEName)
		listGName = [
			#(biopolymer_id,name)
			# nothing has name 'Z'
			(geneBID['A'], 'A'),  (geneBID['A'], 'A2'),
			(geneBID['B'], 'B'),
			(geneBID['C'], 'C'),
			(geneBID['D'], 'D'),  (geneBID['D'], 'DE'),
			(geneBID['E'], 'E'),  (geneBID['E'], 'DE'),  (geneBID['E'], 'EF'),
			(geneBID['F'], 'F'),  (geneBID['F'], 'EF'),  (geneBID['F'], 'FG'),
			(geneBID['G'], 'G'),  (geneBID['G'], 'FG'),
			(geneBID['H'], 'H'),
			(geneBID['I'], 'I'),
			(geneBID['P'], 'P'),
			(geneBID['Q'], 'Q'),
			(geneBID['R'], 'R'),
			(geneBID['S'], 'S'),
		]
		self.addBiopolymerNamespacedNames(namespaceID['gene'], listGName)
		listPName = [
			#(biopolymer_id,name)
			(geneBID['P'],'pqr'),  (geneBID['P'],'qrp'),
			(geneBID['Q'],'pqr'),  (geneBID['Q'],'qrp'),  (geneBID['Q'],'qrs'),
			(geneBID['R'],'pqr'),  (geneBID['R'],'qrp'),  (geneBID['R'],'qrs'),
			                                              (geneBID['S'],'qrs'),
		]
		self.addBiopolymerNamespacedNames(namespaceID['protein'], listPName)
		self.log(" OK: %d identifiers\n" % (len(geneBID)+len(listGName)+len(listPName)))
		
		# TODO: name references?
		
		# define gene regions
		self.log("adding gene regions to the database ...")
		ld0 = ldprofileID['']
		ld1 = ldprofileID['ld']
		listRegion = [
			#(biopolymer_id,ldprofile_id,chr,posMin,posMax)
			(geneBID['A'], ld0, 1,  8, 22),  (geneBID['A'], ld1, 1,  6, 24), # expand both, no gain
			(geneBID['B'], ld0, 1, 28, 52),  (geneBID['B'], ld1, 1, 26, 52), # expand left, no gain
			(geneBID['C'], ld0, 1, 54, 62),  (geneBID['C'], ld1, 1, 48, 64), # expand both, gain dupe
			(geneBID['D'], ld0, 1, 58, 72),  (geneBID['D'], ld1, 1, 54, 74), # expand both, gain 1
			(geneBID['E'], ld0, 1, 78, 82),  (geneBID['E'], ld1, 1, 78, 84), # expand in, no gain
			(geneBID['E'], ld0, 1, 84, 92),  (geneBID['E'], ld1, 1, 84, 94), # expand right, no gain
			(geneBID['F'], ld0, 1, 94, 98),  (geneBID['F'], ld1, 1, 94, 99), # expand right, no gain
			# no regions for G
			(geneBID['H'], ld0, 2, 22, 42),  (geneBID['H'], ld1, 2, 22, 48), # expand to match
			(geneBID['I'], ld0, 2, 38, 48),  (geneBID['I'], ld1, 2, 22, 48), # expand to match
			(geneBID['P'], ld0, 3, 14, 18),  (geneBID['P'], ld1, 3, 16, 22), # expand both, gain 1
			(geneBID['Q'], ld0, 3, 28, 36),  (geneBID['Q'], ld1, 3, 26, 42), # expand both, gain 1 between
			(geneBID['R'], ld0, 3, 44, 52),  (geneBID['R'], ld1, 3, 38, 54), # expand both, gain 1 between
			(geneBID['S'], ld0, 3, 58, 64),  (geneBID['S'], ld1, 3, 56, 72), # expand to dupe
			(geneBID['S'], ld0, 3, 66, 72),  (geneBID['S'], ld1, 3, 56, 72), # expand to dupe
		]
		self.addBiopolymerRegions(listRegion)
		self.log(" OK: %d regions\n" % len(listRegion))
		
		# set the zone size to 7 so that a few things land right on zone edges
		self._loki.setDatabaseSetting("zone_size", "7")
	#update()
	
	
#Source_genes
