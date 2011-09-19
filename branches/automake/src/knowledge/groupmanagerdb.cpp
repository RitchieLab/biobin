/* 
 * File:   groupmanagerdb.cpp
 * Author: torstees
 * 
 * Created on March 23, 2011, 11:22 AM
 */

#include "groupmanagerdb.h"
#include <soci-sqlite3.h>
#include "utility/locus.h"
#include "utility/strings.h"

namespace Knowledge {

uint GroupManagerDB::LoadFromDB(soci::session& sociDB, RegionManager& regions) {
	Utility::IdCollection ids;
	return LoadFromDB(sociDB, ids, regions);
}


uint GroupManagerDB::LoadFromDB(soci::session& sociDB,
	const Utility::IdCollection& pids,
	RegionManager& regions,
	Utility::StringArray& groupNames) {

	Utility::IdCollection ids(pids);		
	std::string groupConstraint = "";
	if (ids.size() > 0)
		groupConstraint = " AND group_id IN (" + Utility::Join(ids, ",") + ") ";
	if (groupNames.size() > 0)
		if (ids.size() > 0)
			groupConstraint = " AND ( group_id IN (" + Utility::Join(ids, ",") + ") OR group_name IN ('" + Utility::Join(groupNames, "', '") + "'))";
		else
			groupConstraint = " AND group_name IN ('" + Utility::Join(groupNames, "', '") + "')";

	soci::rowset<soci::row> rs = (sociDB.prepare<<"SELECT group_id from groups WHERE group_type_id="<<id<<groupConstraint);
	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); itr++) {
		soci::row const& row = *itr;
		uint groupID = row.get<int>(0);
		ids.insert(groupID);
	}

	//Basically, we can't pass an empty id set to LoadFromDB since it will assume we want
	//to load everything. So, this is likely not supposed to have anything in it, or
	//the user misspelled the names
	if (ids.size() == 0)
		return 0;
	return LoadFromDB(sociDB, ids, regions);
}

/**
 * Algorithmically, we want to pull all groups for the list of IDs
 */
uint GroupManagerDB::LoadFromDB(soci::session& sociDB, const Utility::IdCollection& ids, RegionManager& regions) {
	Utility::IdCollection geneIDs = ids;

	//This allows us to capture children of the meta group but also works for the recursive functionality
	if (ids.size() == 0)
		geneIDs.insert(id);
	
	Utility::IdCollection totalGroupsLoaded;

	while (geneIDs.size() > 0) {

		std::string idList = Utility::Join(geneIDs, ",");
		std::string idConstraint = "";
		if (ids.size() > 0) {
			idConstraint = "AND group_id IN (" + idList + ")";
		}
		geneIDs.clear();
		{
			//std::cerr<<"SELECT group_id, group_name, group_desc FROM groups WHERE group_type_id="<<id<<" "<<idConstraint<<"\n";
			soci::rowset<soci::row> rs = (sociDB.prepare<<"SELECT group_id, group_name, group_desc FROM groups WHERE group_type_id="<<id<<" "<<idConstraint);
			//std::cerr<<"ASDFASDFASDF: "<<"SELECT group_id, group_name, group_desc FROM groups WHERE group_type_id="<<id<<" "<<idConstraint<<"\n";
			for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); itr++) {
				soci::row const& row = *itr;
				uint groupID = row.get<int>(0);
				std::string gname = row.get<std::string>(1);
				std::string gdesc = row.get<std::string>(2);
				if (totalGroupsLoaded.find(groupID) == totalGroupsLoaded.end()) {
					AddGroup(groupID, gname.c_str(), gdesc.c_str());
					totalGroupsLoaded.insert(groupID);
					geneIDs.insert(groupID);
				}
			}
		}

		idList = Utility::Join(geneIDs, ",");
		//std::cerr<<"--  "<<idList<<"\n"<<"  --"<<Utility::Join(totalGroupsLoaded, ",")<<"\n";
		geneIDs.clear();
		{
			//std::cerr<<"SELECT parent_id, child_id FROM group_relationships WHERE parent_id IN ("<<idList<<")\n";
			soci::rowset<soci::row> rs = (sociDB.prepare<<"SELECT parent_id, child_id FROM group_relationships WHERE parent_id IN ("<<idList<<")");
			for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); itr++) {
				soci::row const& row = *itr;
				uint parentID	= row.get<int>(0);
				uint childID	= row.get<int>(1);
				AddAssociation(parentID, childID);
				//std::cerr<<" -- "<<parentID<<" : "<<childID<<"\n";

				if (totalGroupsLoaded.find(childID) == totalGroupsLoaded.end())
					geneIDs.insert(childID);
			}
		}
	}
	std::string idList = Utility::Join(totalGroupsLoaded, ",");
	//At this point, we should have the IDs associated with all of the groups, so we can query for the gene associations all at once
	soci::rowset<soci::row> rs = (sociDB.prepare<<"SELECT group_id, gene_id FROM group_associations WHERE group_id IN ("<<idList<<")");
	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); itr++) {
		soci::row const& row = *itr;
		uint groupID	= row.get<int>(0);
		uint geneID		= row.get<int>(1);
		AddGeneAssociation(groupID, regions(geneID), regions);
	}

	return ids.size()==0 ? totalGroupsLoaded.size() - 1 : totalGroupsLoaded.size();		// We don't want to count the Metagroup ID, so subtract one

}


}

#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

class GroupManagerDBTest : public ::testing::Test {
public:
	SnpDataset dataset;
	soci::session sociDB;
	RegionManagerDB regions;
	GroupManagerDB groups;
	std::map<std::string, uint> regionAliasToID;

	GroupManagerDBTest() : groups(1, MetaGroup::DiseaseIndependent, "test1", "a db test") {}

	virtual void SetUp() {
		remove("fake-bio2.db");

		sociDB.open(soci::sqlite3, "dbname=fake-bio2.db timeout=2500");
		sociDB<<"CREATE TABLE regions (gene_id INTEGER UNIQUE PRIMARY KEY, primary_name VARCHAR(128) UNIQUE, chrom VARCHAR(2), description TEXT)";
		sociDB<<"CREATE TABLE region_bounds (gene_id INTEGER, population_id INTEGER, start INTEGER, end INTEGER)";
		sociDB<<"CREATE TABLE region_alias_type (region_alias_type_id INTEGER UNIQUE PRIMARY KEY, region_alias_type_desc TEXT)";
		sociDB<<"CREATE TABLE region_alias (region_alias_type_id INTEGER, alias VARCHAR(64), gene_id INTEGER, gene_count INTEGER, UNIQUE(gene_id, alias, region_alias_type_id))";
		sociDB<<"CREATE TABLE group_type (group_type_id INTEGER, group_type VARCHAR(64), role_id INTEGER, download_date DATE)";
		sociDB<<"CREATE TABLE group_role (role_id INTEGER, group_role VARCHAR(64))";
		sociDB<<"CREATE TABLE groups (group_type_id INTEGER, group_id INTEGER PRIMARY KEY AUTOINCREMENT, group_name VARCHAR(32) UNIQUE, group_desc TEXT)";
		sociDB<<"CREATE TABLE group_relationships (child_id INTEGER, parent_id INTEGER, relationship INTEGER, relationship_description TEXT, UNIQUE(child_id, parent_id, relationship))";
		sociDB<<"CREATE TABLE group_associations (group_id INTEGER, gene_id INTEGER, UNIQUE(group_id, gene_id))";

		//The chromosome and positions are made up-I just want to have a little overlap, but don't feel like typing in large numbers
		//Also, we've already tested ld, so I'm just doing no ld
		sociDB<<"INSERT INTO regions VALUES (0, 'A1BG', 1, 'alpha-1-B glycoprotein')";
		sociDB<<"INSERT INTO region_bounds VALUES (0, 0, 500, 1500)";
		dataset.AddSNP(1, 1000, "rs1000", 1);
		dataset.AddSNP(1, 1100, "rs1001", 1);
		dataset.AddSNP(1, 1200, "rs1002", 1);
		dataset.AddSNP(1, 1300, "rs1003", 1);
		dataset.AddSNP(1, 1350, "rs1004", 1);

		sociDB<<"INSERT INTO regions VALUES (1, 'A2M',  1, 'alpha-2-macroglobulin')";
		sociDB<<"INSERT INTO region_bounds VALUES (1, 0, 1600, 2000)";
		dataset.AddSNP(1, 1700, "rs1010", 1);
		dataset.AddSNP(1, 1900, "rs1011", 1);
		dataset.AddSNP(1, 2000, "rs1012", 1);



		sociDB<<"INSERT INTO regions VALUES (2, 'A2MP1', 1, 'aplha-2-macroglobulin pseudogene')";
		sociDB<<"INSERT INTO region_bounds VALUES (2, 0, 1900, 2200)";
		//rs1011
		//rs1012
		dataset.AddSNP(1, 2050, "rs1013", 1);
		dataset.AddSNP(1, 2100, "rs1014", 1);
		dataset.AddSNP(1, 2150, "rs1015", 1);

		sociDB<<"INSERT INTO regions VALUES (3,'NAT1',1,'N-acetyltransferase 1 (arylamine N-acetyltransferase)')";
		sociDB<<"INSERT INTO region_bounds VALUES (3, 0, 2300, 2400)";
		dataset.AddSNP(1, 2379, "rs1016", 1);

		sociDB<<"INSERT INTO regions VALUES (4,'NAT2',1,'N-acetyltransferase 2 (arylamine N-acetyltransferase)')";
		sociDB<<"INSERT INTO region_bounds VALUES (4, 0, 2500, 2600)";
		dataset.AddSNP(1, 2599, "rs1017", 1);

		sociDB<<"INSERT INTO regions VALUES (5, 'AACCP', 1, 'arylamide acetylase pseudogene')";
		sociDB<<"INSERT INTO region_bounds VALUES (5, 0, 3000, 3500)";
		dataset.AddSNP(1, 3010, "rs1020", 1);
		dataset.AddSNP(1, 3300, "rs1021", 1);
		dataset.AddSNP(1, 3500, "rs1022", 1);

		sociDB<<"INSERT INTO regions VALUES (6, 'SERPINA3', 3, 'serpin peptidase inhibitor, clade A (alpha-1 antiproteinase, antitrypsin), member 3')";
		sociDB<<"INSERT INTO region_bounds VALUES (6, 0, 1000, 2000)";
		dataset.AddSNP(3, 1080, "rs1025", 1);
		dataset.AddSNP(3, 1200, "rs1026", 1);
		dataset.AddSNP(3, 1400, "rs1027", 1);
		dataset.AddSNP(3, 1500, "rs1028", 1);
		dataset.AddSNP(3, 1789, "rs1029", 1);
		dataset.AddSNP(3, 1999, "rs1030", 1);

		sociDB<<"INSERT INTO regions VALUES (7, 'AADAC', 3, 'arylacetamide deacetylase (esterase)')";
		sociDB<<"INSERT INTO region_bounds VALUES (7, 0, 2100, 2500)";
		dataset.AddSNP(3, 2111, "rs1031", 1);
		dataset.AddSNP(3, 2180, "rs1032", 1);
		dataset.AddSNP(3, 2200, "rs1033", 1);
		dataset.AddSNP(3, 2250, "rs1034", 1);
		dataset.AddSNP(3, 2300, "rs1035", 1);
		dataset.AddSNP(3, 2330, "rs1036", 1);
		dataset.AddSNP(3, 2380, "rs1037", 1);
		dataset.AddSNP(3, 2400, "rs1038", 1);
		dataset.AddSNP(3, 2409, "rs1039", 1);
		dataset.AddSNP(3, 2499, "rs1040", 1);


		sociDB<<"INSERT INTO regions VALUES (8, 'AAMP', 3, 'angio-associated, migratory cell protein')";
		sociDB<<"INSERT INTO region_bounds VALUES (8, 0, 2400, 2800)";
		//rs1038
		//rs1039
		//rs1040
		dataset.AddSNP(3, 2600, "rs1041", 1);
		dataset.AddSNP(3, 2750, "rs1042", 1);

		sociDB<<"INSERT INTO regions VALUES (9, 'AANAT', 3, 'angio-associated, migratory cell protein')";
		sociDB<<"INSERT INTO region_bounds VALUES (9, 0, 3000, 3500)";
		dataset.AddSNP(3, 3250, "rs1045", 1);
		dataset.AddSNP(3, 3400, "rs1046", 1);

		sociDB<<"INSERT INTO regions VALUES (10, 'AARS', 3, 'angio-associated, migratory cell protein')";
		sociDB<<"INSERT INTO region_bounds VALUES (10, 0, 4000, 5000)";
		dataset.AddSNP(3, 4400, "rs1050", 1);
		dataset.AddSNP(3, 4800, "rs1051", 1);


		regions.LoadFromDB(sociDB, 0);
		regions.AssociateSNPs(dataset);

		sociDB<<"INSERT INTO group_role VALUES (1, 'Disease Independent')";
		sociDB<<"INSERT INTO group_role VALUES (2, 'Disease Dependent')";
		sociDB<<"INSERT INTO group_role VALUES (3, 'SNP Collection')";
		sociDB<<"INSERT INTO group_role VALUES (4, 'Gene Collection')";

		sociDB<<"INSERT INTO group_type(group_type_id, group_type, role_id) VALUES(1, 'TEST', 1)";

		sociDB<<"INSERT INTO groups VALUES (1, 2, 'A', 'first root')";
		sociDB<<"INSERT INTO groups VALUES (1, 3, 'B', 'second root')";
		sociDB<<"INSERT INTO groups VALUES (1, 100, 'a1', 'child of A')";
		sociDB<<"INSERT INTO groups VALUES (1, 101, 'a2', 'child of A')";
		sociDB<<"INSERT INTO groups VALUES (1, 102, 'a3', 'child of A')";
		sociDB<<"INSERT INTO groups VALUES (1, 110, 'a1a', 'child of a1')";
		sociDB<<"INSERT INTO groups VALUES (1, 111, 'a1b', 'child of a1')";
		sociDB<<"INSERT INTO groups VALUES (1, 112, 'a1b1', 'child of a1b')";
		sociDB<<"INSERT INTO groups VALUES (1, 113, 'a1b2', 'child of a1b2')";
		sociDB<<"INSERT INTO groups VALUES (1, 114, 'a1b1a', 'child of a1b1')";
		sociDB<<"INSERT INTO groups VALUES (1, 115, 'a1b1b', 'child of a1b1')";
		sociDB<<"INSERT INTO groups VALUES (1, 116, 'a2a', 'child of a2')";
		sociDB<<"INSERT INTO groups VALUES (1, 117, 'a2a1', 'child of a2a')";
		sociDB<<"INSERT INTO groups VALUES (1, 118, 'a2a2', 'child of a2a')";
		sociDB<<"INSERT INTO groups VALUES (1, 119, 'a3a', 'child of a3')";
		sociDB<<"INSERT INTO groups VALUES (1, 120, 'a3b', 'child of a3')";
		sociDB<<"INSERT INTO groups VALUES (1, 121, 'a3c', 'child of a3')";
		sociDB<<"INSERT INTO groups VALUES (1, 122, 'a3c1', 'child of a3c')";
		sociDB<<"INSERT INTO groups VALUES (1, 123, 'a3c2', 'child of a3c')";
		sociDB<<"INSERT INTO groups VALUES (1, 124, 'a3c3', 'child of a3c')";
		sociDB<<"INSERT INTO groups VALUES (1, 125, 'a3c2a', 'child of a3c2')";
		sociDB<<"INSERT INTO groups VALUES (1, 126, 'a3c2b', 'child of a3c2')";
		sociDB<<"INSERT INTO groups VALUES (1, 200, 'b1', 'child of B')";
		sociDB<<"INSERT INTO groups VALUES (1, 201, 'b2', 'child of B')";
		sociDB<<"INSERT INTO groups VALUES (1, 202, 'b3', 'child of B')";
		sociDB<<"INSERT INTO groups VALUES (1, 203, 'b4', 'child of B')";
		sociDB<<"INSERT INTO groups VALUES (1, 204, 'b5', 'child of B')";
		sociDB<<"INSERT INTO groups VALUES (1, 205, 'b3a', 'child of B3')";
		sociDB<<"INSERT INTO groups VALUES (1, 206, 'b3b', 'child of B3')";
		sociDB<<"INSERT INTO groups VALUES (1, 207, 'b3c', 'child of B3')";
		sociDB<<"INSERT INTO groups VALUES (1, 208, 'b3d', 'child of B3')";
		sociDB<<"INSERT INTO groups VALUES (1, 209, 'b3e', 'child of B3')";
		//sociDB<<"INSERT INTO groups VALUES (1, 210, 'b3e', 'child of B3')";

		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (1, 2, 'root', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (1, 3, 'root', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (2, 100, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (2, 101, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (2, 102, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (3, 200, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (3, 201, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (3, 202, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (3, 203, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (3, 204, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (100, 110, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (101, 116, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (102, 119, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (102, 120, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (102, 121, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (110, 111, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (110, 113, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (112, 114, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (112, 115, '', '')";
		//sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (100, 110, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (116, 117, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (116, 118, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (121, 122, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (121, 123, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (121, 124, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (123, 125, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (123, 126, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (202, 205, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (202, 206, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (202, 207, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (202, 208, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (202, 209, '', '')";
		sociDB<<"INSERT INTO group_relationships(parent_id, child_id, relationship, relationship_description) VALUES (204, 300, '', '')";

		sociDB<<"INSERT INTO group_associations VALUES (110, 1)";
		sociDB<<"INSERT INTO group_associations VALUES (110, 2)";
		sociDB<<"INSERT INTO group_associations VALUES (110, 3)";
		sociDB<<"INSERT INTO group_associations VALUES (114, 0)";
		sociDB<<"INSERT INTO group_associations VALUES (114, 2)";
		sociDB<<"INSERT INTO group_associations VALUES (115, 0)";
		sociDB<<"INSERT INTO group_associations VALUES (115, 3)";
		sociDB<<"INSERT INTO group_associations VALUES (113, 2)";
		sociDB<<"INSERT INTO group_associations VALUES (113, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (113, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (119, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (119, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (120, 3)";
		sociDB<<"INSERT INTO group_associations VALUES (120, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (122, 2)";
		sociDB<<"INSERT INTO group_associations VALUES (122, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (125, 2)";
		sociDB<<"INSERT INTO group_associations VALUES (125, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (126, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (126, 7)";
		sociDB<<"INSERT INTO group_associations VALUES (124, 6)";

		sociDB<<"INSERT INTO group_associations VALUES (200, 3)";
		sociDB<<"INSERT INTO group_associations VALUES (200, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (201, 3)";
		sociDB<<"INSERT INTO group_associations VALUES (201, 10)";
		sociDB<<"INSERT INTO group_associations VALUES (205, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (205, 6)";
		sociDB<<"INSERT INTO group_associations VALUES (206, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (206, 7)";
		sociDB<<"INSERT INTO group_associations VALUES (206, 8)";
		sociDB<<"INSERT INTO group_associations VALUES (207, 6)";
		sociDB<<"INSERT INTO group_associations VALUES (207, 8)";
		sociDB<<"INSERT INTO group_associations VALUES (208, 5)";
		sociDB<<"INSERT INTO group_associations VALUES (208, 7)";
		sociDB<<"INSERT INTO group_associations VALUES (208, 9)";
		sociDB<<"INSERT INTO group_associations VALUES (209, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (209, 6)";
		sociDB<<"INSERT INTO group_associations VALUES (209, 7)";
		sociDB<<"INSERT INTO group_associations VALUES (203, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (203, 10)";
		sociDB<<"INSERT INTO group_associations VALUES (204, 3)";
		sociDB<<"INSERT INTO group_associations VALUES (204, 4)";
		sociDB<<"INSERT INTO group_associations VALUES (204, 10)";


		//std::map<uint, uint> geneIdToIndex;
		//regions.GenerateLookupTable(geneIdToIndex);
		groups.LoadFromDB(sociDB, regions);

	}

	virtual void TearDown() {
		remove("fake-bio2.db");
	}
};


TEST_F(GroupManagerDBTest, VerifyLoading) {
	EXPECT_EQ(32, groups.Size());
	EXPECT_EQ(1, groups.id);
	EXPECT_EQ("test1", groups.name);
	EXPECT_EQ("a db test", groups.description);

	EXPECT_EQ("A", groups[0].name);
	EXPECT_EQ("first root", groups[0].desc);
	EXPECT_EQ(3, groups[1].id);
	EXPECT_EQ("a1", groups(100).name);

	EXPECT_EQ(119, groups(119).id);
	EXPECT_EQ("a3a", groups(119).name);

	EXPECT_EQ(3, groups(2).groups.size());
	EXPECT_EQ(1, groups(100).groups.size());

	EXPECT_EQ(119, groups(119).id);
	EXPECT_EQ(3, groups(102).groups.size());
	EXPECT_EQ(3, groups(121).groups.size());
	EXPECT_EQ(5, groups(3).groups.size());


	EXPECT_EQ(3, groups(110).regions.size());
	EXPECT_EQ(3, groups(113).regions.size());
	EXPECT_EQ(119, groups(119).id);
	EXPECT_EQ(2, groups(119).regions.size());
	EXPECT_EQ(2, groups(125).regions.size());

	GroupManager mgr2(groups);
	EXPECT_EQ(1, mgr2.id);
	EXPECT_EQ("test1", mgr2.name);
	EXPECT_EQ("a db test", mgr2.description);

	EXPECT_EQ(3, mgr2(121).groups.size());
	EXPECT_EQ(5, mgr2(3).groups.size());


	EXPECT_EQ(3, mgr2(110).regions.size());
	EXPECT_EQ(3, mgr2(113).regions.size());
}
#endif
