#!/usr/bin/env python

from loki import loki_source


class Source_snps(loki_source.Source):
	
	
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
		
		# define positions
		self.log("adding SNPs to the database ...")
		listSNP = [
			#(rs,chr,pos,valid)
			(11, 1, 10, 1),
			(12, 1, 20, 1),
			(13, 1, 35, 1),
			(14, 1, 35, 1),
			(15, 1, 50, 1),
			(15, 1, 55, 1),
			(16, 1, 60, 1),
			(17, 1, 70, 1),
			(18, 1, 80, 1),
			(19, 1, 90, 1),
			(21, 2, 10, 1),
			(22, 2, 20, 1),
			(23, 2, 30, 0),
			(24, 2, 40, 1),
			(25, 2, 50, 1),
			(31, 3, 10, 0),
			(32, 3, 20, 1),
			(33, 3, 30, 1),
			(34, 3, 40, 1),
			(35, 3, 50, 1),
			(36, 3, 60, 1),
			(37, 3, 70, 0),
		]
		self.addSNPLoci(listSNP)
		self.log(" OK: %d SNP positions (%d RS#s)\n" % (len(listSNP),len(set(s[0] for s in listSNP))))
		
		# define merges
		self.log("adding SNP merge records to the database ...")
		listMerge = [
			#(rsOld,rsNew)
			(9,19),
		]
		self.addSNPMerges(listMerge)
		self.log(" OK: %d merges\n" % len(listMerge))
		
		# define role codes
		self.log("adding SNP role codes to the database ...")
		listRole = [
			#(role,desc,coding,exon)
			('exon',   'exon',                1, 1),
			('utr',    'untranslated region', 0, 1),
			('intron', 'intron',              0, 0),
			('reg',    'regulatory',          1, 0),
		]
		roleID = self.addRoles(listRole)
		self.log(" OK: %d role codes\n" % len(roleID))
		
		# define SNP roles
		self.log("adding SNP roles to the database ...")
		listSNPRole = [
			#(rs,entrez_id,role_id)
			(11,0,roleID['reg']),
			(12,1,roleID['exon']),
			(13,2,roleID['utr']),
			(13,2,roleID['intron']),
			# no role for rs14 which overlaps rs13
			(15,2,roleID['reg']),
			(15,3,roleID['exon']),
			(16,3,roleID['intron']),
			(16,4,roleID['intron']),
			(17,3,roleID['reg']),
			(18,5,roleID['exon']),
			( 9,5,roleID['exon']), # rs9 merged -> rs19
			# no role for rs21
			(22,8,roleID['utr']),
			(23,8,roleID['intron']),
			(24,8,roleID['reg']),
			(24,9,roleID['exon']),
			(25,16,roleID['reg']),
			(36,19,roleID['intron']),
			(37,19,roleID['exon']),
		]
		self.addSNPEntrezRoles(listSNPRole)
		self.log(" OK: %d roles\n" % len(listSNPRole))
	#update()
	
	
#Source_snps
