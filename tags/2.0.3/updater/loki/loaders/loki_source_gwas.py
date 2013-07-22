#!/usr/bin/env python

import re
from loki import loki_source


class Source_gwas(loki_source.Source):
	
	
	##################################################
	# source interface
	
	
	@classmethod
	def getVersionString(cls):
		return '2.1 (2013-06-11)'
	#getVersionString()
	
	
	def download(self, options):
		# download the latest source files
		self.downloadFilesFromHTTP('www.genome.gov', {
			'gwascatalog.txt': '/admin/gwascatalog.txt',
		})
	#download()
	
	
	def update(self, options):
		# clear out all old data from this source
		self.log("deleting old records from the database ...")
		self.deleteAll()
		self.log(" OK\n")
		
		# process gwas cataog
		self.log("processing GWAS catalog annotations ...")
		reRS = re.compile('rs[0-9]+', re.I)
		listNone = [None]
		numInc = 0
		setGwas = set()
		with open('gwascatalog.txt','rU') as gwasFile:
			header = gwasFile.next().rstrip()
			if not header.startswith("Date Added to Catalog\tPUBMEDID\tFirst Author\tDate\tJournal\tLink\tStudy\tDisease/Trait\tInitial Sample Size\tReplication Sample Size\tRegion\tChr_id\tChr_pos\tReported Gene(s)\tMapped_gene\tUpstream_gene_id\tDownstream_gene_id\tSnp_gene_ids\tUpstream_gene_distance\tDownstream_gene_distance\tStrongest SNP-Risk Allele\tSNPs\tMerged\tSnp_id_current\tContext\tIntergenic\tRisk Allele Frequency\tp-Value\tPvalue_mlog\tp-Value (text)\tOR or beta\t95% CI (text)\t"):#Platform [SNPs passing QC]\tCNV"):
				self.log(" ERROR\n")
				self.log("unrecognized file header: %s\n" % (header,))
				return False
			for line in gwasFile:
				line = line.rstrip("\r\n")
				words = list(w.strip() for w in line.decode('latin-1').split("\t"))
				if len(words) <= 31:
					# blank line at the end is normal
					if (len(words) > 1) or words[0]:
						numInc += 1
					continue
				chm = self._loki.chr_num[words[11]] if (words[11] in self._loki.chr_num) else None
				pos = long(words[12]) if words[12] else None
				trait = words[7]
				snps = words[21] if words[20].endswith('aplotype') else words[20]
				rses = list(int(rs[2:]) for rs in reRS.findall(snps)) or listNone
				orBeta = words[30]
				allele95ci = words[31]
				riskAfreq = words[26]
				pubmedID = int(words[1]) if words[1] else None
				for rs in rses:
					setGwas.add( (rs,chm,pos,trait,snps,orBeta,allele95ci,riskAfreq,pubmedID) )
			#foreach line
		#with gwasFile
		self.log(" OK: %d catalog entries (%d incomplete)\n" % (len(setGwas),numInc))
		if setGwas:
			self.log("writing GWAS catalog annotations to the database ...")
			self.addGWASAnnotations(setGwas)
			self.log(" OK\n")
	#update()
	
#Source_gwas
