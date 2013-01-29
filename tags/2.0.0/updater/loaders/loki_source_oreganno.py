import loki_source
import itertools

class Source_oreganno(loki_source.Source):
	
	_remHost = "hgdownload.cse.ucsc.edu"
	_remPath = "/goldenPath/hg19/database/"
	
	_remFiles = ["oreganno.txt.gz","oregannoAttr.txt.gz", "oregannoLink.txt.gz"]
	
	
	@classmethod
	def getVersionString(cls):
		return '2.0a1 (2012-07-06)'
	#getVersionString()
	
	
	def download(self, options):
		"""
		Download OregAnno from UCSC
		"""
		self.downloadFilesFromFTP(self._remHost, dict(((f, self._remPath + f) for f in self._remFiles)))
		
	def update(self, options):
		"""
		Update the database with the OregAnno data from ucsc
		"""
			
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")		
		
		# Add the 'oreganno' namespace
		ns = self.addNamespace('oreganno')
		
		# Add the ensembl and entrez namespaces
		external_ns = self.addNamespaces([
			('symbol',      0),
			('entrez_gid',  0),
			('ensembl_gid', 0),
		])
		
		# Add the types of Regions
		typeids = self.addTypes([
			('regulatory_region',),
			('tfbs',),
			('gene',)
		])
		
		# Add the types of groups
		group_typeid = self.addType('regulatory_group')
		
		# Add the role for regulatory
		snp_roleid = self.addRole("regulatory", "OregAnno Regulatory Polymorphism", 1, 0)
		
		# Get the default population ID
		ldprofile_id = self.addLDProfile('', 'no LD adjustment')
		
		# build dict of gene id->oreganno id and a dict of
		# oreganno id->entrez id and oreganno id->ensembl id
		oreg_gene = {}
		oreg_tfbs = {}
		oreg_snp = {}
		link_f = self.zfile("oregannoLink.txt.gz")
		entrez_ns = external_ns['entrez_gid']
		ensembl_ns = external_ns['ensembl_gid']
		symbol_ns = external_ns['symbol']
		self.log("parsing external links ...")
		for l in link_f:
			fields = l.split()
			if fields[1] == "Gene":
				oreg_id = fields[0]
				if fields[2] == "EnsemblGene":
					gene_id = fields[3].split(',')[1]
					oreg_gene.setdefault(oreg_id,{})[ensembl_ns] = gene_id
				elif fields[2] == "EntrezGene":
					gene_id = fields[3]
					oreg_gene.setdefault(oreg_id,{})[entrez_ns] = gene_id
			elif fields[1] == "TFbs":
				oreg_id = fields[0]
				if fields[2] == "EnsemblGene":
					gene_id = fields[3].split(',')[1]
					oreg_tfbs.setdefault(oreg_id,{})[ensembl_ns] = gene_id
				elif fields[2] == "EntrezGene":
					gene_id = fields[3]
					oreg_tfbs.setdefault(oreg_id,{})[entrez_ns] = gene_id
			elif fields[1] == "ExtLink" and fields[2] == "dbSNP":
				# Just store the RS# (no leading "rs")
				oreg_snp[fields[0]] = fields[3][2:]
		
		self.log("OK\n")
		
		# Now, create a dict of oreganno id->type
		oreganno_type = {}
		self.log("parsing region attributes ...")
		attr_f = self.zfile("oregannoAttr.txt.gz")
		for l in attr_f:
			fields = l.split('\t')
			if fields[1] == "type":
				oreganno_type[fields[0]] = fields[2]
			elif fields[1] == "Gene":
				oreg_gene.setdefault(fields[0],{})[symbol_ns] = fields[2]
			elif fields[1] == "TFbs":
				oreg_tfbs.setdefault(fields[0],{})[symbol_ns] = fields[2]
				
		
		self.log("OK\n")
		
		# OK, now parse the actual regions themselves
		region_f = self.zfile("oreganno.txt.gz")
		oreganno_roles = []
		oreganno_regions = []
		oreganno_bounds = []
		oreganno_groups = {}
		oreganno_types = {}
		self.log("parsing regulatory regions ...")
		snps_unmapped = 0
		for l in region_f:
			fields = l.split()
			chrom = self._loki.chr_num.get(fields[1][3:])
			start = int(fields[2]) + 1
			stop = int(fields[3])
			oreg_id = fields[4]
			oreg_type = oreganno_type[oreg_id]
			if chrom and oreg_type == "REGULATORY POLYMORPHISM":
				entrez_id = oreg_gene.get(oreg_id,{}).get(entrez_ns)
				rsid = oreg_snp.get(oreg_id)
				if entrez_id and rsid:
					oreganno_roles.append((int(rsid), entrez_id, snp_roleid))
				else:
					snps_unmapped+=1
			elif chrom and (oreg_type == "REGULATORY REGION" or oreg_type == "TRANSCRIPTION FACTOR BINDING SITE"):
				gene_symbol = oreg_gene.get(oreg_id,{}).get(symbol_ns)
				if not gene_symbol:
					gene_symbol = oreg_tfbs.get(oreg_id,{}).get(symbol_ns)
				
				if gene_symbol:
					oreganno_groups.setdefault(gene_symbol, []).append(oreg_id)
				
				if oreg_type == "REGULATORY REGION":
					oreg_typeid = typeids['regulatory_region']
				else:
					oreg_typeid = typeids['tfbs']
				
				oreganno_types[oreg_id] = oreg_typeid
				oreganno_regions.append((oreg_typeid, oreg_id, ''))
				oreganno_bounds.append((chrom, start, stop))
				
		self.log("OK (%d regions found, %d SNPs found, %d SNPs unmapped)\n" % (len(oreganno_regions), len(oreganno_roles), snps_unmapped))
	
		self.log("Writing to database ... ")
		self.addSNPEntrezRoles(oreganno_roles)
		reg_ids = self.addBiopolymers(oreganno_regions)
		self.addBiopolymerNamespacedNames(ns, ((reg_ids[i], oreganno_regions[i][1]) for i in range(len(reg_ids))))
		bound_gen = zip(((r,) for r in reg_ids),oreganno_bounds)
		self.addBiopolymerLDProfileRegions(ldprofile_id, ((itertools.chain(*c) for c in bound_gen)))
		
		# Now, add the regulation groups
		oreg_genes = oreganno_groups.keys()
		oreg_gids = self.addTypedGroups(group_typeid, (("regulatory_%s" % k, "OregAnno Regulation of %s" % k) for k in oreg_genes))
		self.addGroupNamespacedNames(ns, zip(oreg_gids, ("regulatory_%s" % k for k in oreg_genes)))
		
		group_membership = []
		for i in range(len(oreg_gids)):
			gid = oreg_gids[i]
			gene_key = oreg_genes[i]
			gene_member = set()
			tfbs_member = {}
			member_num = 2
			for oreg_id in oreganno_groups[gene_key]:
				member_num += 1
				group_membership.append((gid, member_num, oreganno_types.get(oreg_id, 0), ns, oreg_id))
				for external_nsid, external_val in oreg_gene.get(oreg_id,{}).iteritems():
					gene_member.add((gid, 1, typeids['gene'], external_nsid, external_val))
					
				member_num += 1
				for external_nsid, external_val in oreg_tfbs.get(oreg_id,{}).iteritems():
					tfbs_member.setdefault(external_nsid,{})[external_val] = member_num
				
			group_membership.extend(gene_member)
			for ext_ns, d in tfbs_member.iteritems():
				for sym, mn in d.iteritems():
					group_membership.append((gid, mn, typeids['gene'], ext_ns, sym))
		
		self.addGroupMemberNames(group_membership)
	
		self.log("OK\n")
