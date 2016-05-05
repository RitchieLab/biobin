#!/usr/bin/env python

from loki import loki_source


class Source_paint(loki_source.Source):
	
	
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
		namespaceID = self.addNamespaces([
			('gene',  0),
			('group', 0),
		])
		relationshipID = self.addRelationships([
			('different_than',),
		])
		typeID = self.addTypes([
			('gene',),
			('group',),
		])
		
		# define groups
		self.log("adding groups to the database ...")
		listGroup = [
			#(label,description)
			('cyan',    'gene ambiguity resolved by either heuristic'),
			('magenta', 'gene ambiguity resolved only by implication heuristic'),
			('yellow',  'gene ambiguity resolved only by quality heuristic'),
			('gray',    'unresolvable gene ambiguity'),
		]
		listGID = self.addTypedGroups(typeID['group'], listGroup)
		groupGID = dict(zip((g[0] for g in listGroup), listGID))
		self.log(" OK: %d groups\n" % len(groupGID))
		
		# define group names
		self.log("adding group names to the database ...")
		listName = [
			#(group_id,name)
			(groupGID['cyan'],    'cyan'),
			(groupGID['magenta'], 'magenta'),
			(groupGID['yellow'],  'yellow'),
			(groupGID['gray'],    'gray'),
			(groupGID['gray'],    'black'),
		]
		self.addGroupNamespacedNames(namespaceID['group'], listName)
		self.log(" OK: %d names\n" % len(listName))
		
		# define group relationships
		self.log("adding group relationships to the database ...")
		listRel = [
			#(group_id,related_group_id,relationship_id)
			(groupGID['cyan'],    groupGID['magenta'], relationshipID['different_than']),
			(groupGID['magenta'], groupGID['yellow'],  relationshipID['different_than']),
			(groupGID['yellow'],  groupGID['cyan'],    relationshipID['different_than']),
		]
		self.addGroupSiblingRelationships(listRel)
		self.log(" OK: %d relationships\n" % len(listRel))
		
		# define group members
		self.log("adding group members to the database ...")
		listMember = [
			#(group_id,member,name)
			(groupGID['cyan'],    11, 'A2'),
			(groupGID['cyan'],    12, 'C'),
			(groupGID['cyan'],    13, 'D'),
			(groupGID['cyan'],    13, 'DE'),
			(groupGID['magenta'], 21, 'DE'),
			(groupGID['magenta'], 21, 'EF'),
			(groupGID['magenta'], 21, 'G'),
			(groupGID['yellow'],  31, 'EF'),
			(groupGID['yellow'],  31, 'FG'),
			(groupGID['yellow'],  31, 'G'),
			(groupGID['gray'],    41, 'F'),
			(groupGID['gray'],    41, 'FG'),
			(groupGID['gray'],    41, 'G'),
		]
		self.addGroupMemberTypedNamespacedNames(typeID['gene'], namespaceID['gene'], listMember)
		self.log(" OK: %d members (%d identifiers)\n" % (len(set(m[1] for m in listMember)),len(listMember)))
	#update()
	
	
#Source_paint
