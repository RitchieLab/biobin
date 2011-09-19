/* 
 * File:   groupmanager.cpp
 * Author: torstees
 * 
 * Created on March 11, 2011, 12:13 PM
 */

#include "groupmanager.h"
#include <soci-sqlite3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "utility/exception.h"

namespace Knowledge {

uint GroupManager::maxGeneCount = 30;

void GroupManager::GenerateGeneGeneModels(GeneGeneModelArchive& geneArchive, 
	 std::map<uint, Utility::IdCollection>& regionLookup,
	 RegionManager& regions,
	 uint idx,
	 std::ostream& os,
	 Utility::IdCollection& visited) {

	Group &group = groups[idx];
	if (visited.find(idx) == visited.end()) {
		visited.insert(idx);
		if (regionLookup[idx].size() > maxGeneCount) {
			Utility::IdCollection::iterator itr = group.groups.begin();
			Utility::IdCollection::iterator end = group.groups.end();

			while (itr != end) {
				GenerateGeneGeneModels(geneArchive, regionLookup, regions, *itr++, os, visited);
			}
		} else {
			// Let's start creating gene gene models
			geneArchive.AddRegions(regionLookup[idx], regions);
		}
	}
}

void GroupManager::GenerateGeneGeneModels(GeneGeneModelArchive& geneArchive, RegionManager& regions, std::ostream& os) {

	//If group Type is one of the dataset types, we don't want to generate models 
	if (groupType == MetaGroup::DiseaseIndependent || groupType == MetaGroup::DiseaseDependent) {
		std::map<uint, Utility::IdCollection > regionLookup;
		BuildRegionCollections(regionLookup);

		Utility::IdCollection::iterator itr = root.begin();
		Utility::IdCollection::iterator end = root.end();
		Utility::IdCollection visited;
		while (itr != end)
			GenerateGeneGeneModels(geneArchive, regionLookup, regions, *itr++, os, visited);
	}
}

/* 
void GroupManager::CollectLeafIndexes(uint snpThreshold, 
			RegionManager& regions, 
			Utility::IdCollection& groupIDs,
		std::map<uint, Utility::IdCollection>& regionLookup) {
	BuildRegionCollections(regionLookup);

	Utility::IdCollection::iterator itr = root.begin();
	Utility::IdCollection::iterator end = root.end();
	Utility::IdCollection visited;
	while (itr != end)
		CollectLeafIndexes(snpThreshold, itr++, regions, groupIDs, regionLookup, visited);
}

void GroupManager::CollectLeafIndexes(uint snpThreshold, 
			uint groupIdx, 
			RegionManager& regions, 
			Utility::IdCollection& groupIdxs, 
			std::map<uint, Utility::IdCollection >& regionLookup, 
			Utility::IdCollection& visited) {
	
	if (visited.find(groupIdx) == visited.end()) {
		visited.insert(groupIdx);
		Utility::IdCollection& geneIDs = regionLookup[groupIdx];
		uint snpCount = 0;
		
		Utility::IdCollection snps;
		Utility::IdCollection geneIdx = geneIDs.begin();
		Utility::IdCollection geneEnd = geneIDs.end();
		
		while (geneIdx != geneEnd) {
			Region& r = regions[*geneIdx++];
			snps.insert(r.snps.begin(), r.snps.end());
		}
		Group& group = groups[groupIdx];
		// If we are still too big and we have children...
		if (snps.size() > snpThreshold && group.groups.size() > 0) {
			Utility::IdCollection groupItr = group.groups.begin();
			Utility::IdCollection groupEnd = group.groups.end();
			
			while (groupItr != groupEnd) 
				CollectLeafIndexes(snpThreshold, *groupItr++, regions, groupIdxs, regionLookup, visited);
			
		} else {
			groupIdxs.insert(groupIdx);
		}
		
	}	
}
*/
Utility::IdCollection GroupManager::Root() {
	return root;
}

void GroupManager::BuildRegionCollections(uint idx, std::map<uint, Utility::IdCollection >& regions) {
	// Make sure we haven't already visited this one to avoid infinite loops
	if (regions.find(idx) == regions.end()) {
		Utility::IdCollection childRegions;
		regions[idx] = childRegions;

		Utility::IdCollection::iterator itr = groups[idx].groups.begin();
		Utility::IdCollection::iterator end = groups[idx].groups.end();

		while (itr != end) {
			BuildRegionCollections(*itr, regions);
			Utility::IdCollection &localRegions = regions[*itr];
			childRegions.insert(localRegions.begin(), localRegions.end());
			itr++;
		}

		childRegions.insert(groups[idx].regions.begin(), groups[idx].regions.end());

		regions[idx] = childRegions;
	}
}

void GroupManager::BuildRegionCollections(std::map<uint, Utility::IdCollection >& regions) {
	uint groupCount = groups.size();

	for (uint i=0; i<groupCount; i++) {
		BuildRegionCollections(i, regions);
	}
}

void GroupManager::AddGroup(uint id, const char *name, const char *desc) {
	if (id != this->id) {
		Group g(id, name, desc);
		groupLookup[id]	= groups.size();
		idxByName[name]	= groups.size();
		groups.push_back(g);
	} 
	// Let's ignore the ones that are the same as their group type. Those are just
	// there to keep things lined up properly
}

void GroupManager::AddAssociation(uint groupID, uint childID) {
	std::map<uint, uint>::iterator missing = groupLookup.end();
	if (groupID == id) {
		if (groupLookup.find(childID) != missing)
			root.insert(groupLookup[childID]);
	} else {
		if (groupLookup.find(groupID) != missing && groupLookup.find(childID) != missing) {
			int pIndex = groupLookup[groupID];
			int cIndex = groupLookup[childID];

			if (pIndex > -1 && cIndex > -1)
				groups[pIndex].groups.insert(cIndex);
		}
	}
}

void GroupManager::AddGeneAssociation(uint groupID, uint regionIndex, RegionManager& regions) {
	std::map<uint, uint>::iterator group = groupLookup.find(groupID);

	// Make sure that the group isn't one that we've decided to drop for some reason
	if (regionIndex < (uint)-1 && group != groupLookup.end() && group->second >= 0) {
		groups[group->second].regions.insert(regionIndex);
		regions.AddMetaID(id, groupType, regionIndex);
	} /* This is for debugging-I think it's related to purging pre-existing DB and am looking into it-EST july 2011
	else {
		std::cerr<<"The following id was said to be associated with a gene, but we couldn't find it:  "<<groupID<<"\n";
	} */
}

/* Stripping this out because it needs to be rethought...the children indexes must be updated, but efficiently
void GroupManager::Squeeze() {
	std::vector<Group> newGroups;
	std::vector<Group>::iterator itr = groups.begin();
	std::vector<Group>::iterator end = groups.end();
	groupLookup.clear();


	while (itr != end) {
		if (itr->groups.size() > 0 || itr->regions.size() > 0) {
			std::cerr<<"  --> "<<itr->id<<"\n";
			groupLookup[itr->id] = newGroups.size();
			newGroups.push_back(*itr);
		}
		else
			std::cerr<<"!!Squeezing out: "<<itr->id<<"\n";
		itr++;
	}

	groups = newGroups;
}*/

Group& GroupManager::operator[](uint idx) {
	return groups[idx];
}

Group& GroupManager::operator()(uint id) {
	return groups[groupLookup[id]];
}

Group& GroupManager::operator[](std::string& name) {
	std::map<std::string, uint>::iterator itrIdx = idxByName.find(name);
	if (itrIdx == idxByName.end())
		throw Utility::Exception::MissingKey("idxByName", name.c_str(), idxByName.size());
	return groups[itrIdx->second];
}
/**
 * File Format:
 * [Meta-Group Name] Description (can have many words)
 * GROUP [name] Description (can have many words)
 * gene1
 * gene2
 * ...
 * GROUP name Description
 * gene1
 * gene2
 *
 * Genes are on separate lines
 */
/**
 * Loads the contents of a group manager text file
 * @param filename
 * @param regionAliasLookup
 */
std::map<std::string, uint> GroupManager::LoadArchive(const char *filename, 
	 RegionManager& regions,
	 uint id,
	 Utility::StringArray& unmatchedAliases) {
/*	std::ifstream file(filename);

	if (!file.good()) {
		std::cerr<<"Unable to open file, "<<filename<<". Disease dependent data was not loaded.\n";
		exit(1);
	}
*/


	std::string contents = Utility::LoadContents(filename);
	Utility::StringArray strings = Utility::StringSplit(contents.c_str(), "\nGROUP");
	Utility::StringArray identity = Utility::Split(strings[0].c_str(), " ");
	name = identity[0];
	groupType = MetaGroup::ConvertType(identity[1].c_str());
	if (groupType == (uint)-1)
		throw Utility::Exception::General((std::string("Unknown Group type: ") + identity[1]).c_str());
	
	description = Utility::Join(identity, " ", 2, -1);

	uint count = strings.size();

	std::map<std::string, uint> nameToId;		///< Used to look up the index based on the group name
	nameToId[name] = id;
	Utility::StringArray childrenDesc;			///< Captured during parsing-will be parsed at very end

	for (uint i=1; i<count; i++) {
		Utility::StringArray data = Utility::Split(strings[i].c_str(), "\n");
		identity = Utility::Split(data[0].c_str(), " ");
		std::string gname = identity[0];
		std::string gdesc = Utility::Join(identity, " ", 1, -1);
		Group group(id++, gname, gdesc);
		nameToId[gname] = group.id;
		//Add the genes to the group
		uint size = data.size();
		for (uint n=1; n<size; n++) {
			std::string alias = data[n];
			Utility::StringArray words = Utility::Split(alias.c_str(), " ");
			if (words.size() > 0) {
				if (words[0] == "CHILDREN") {
					childrenDesc.push_back(alias);
				}
				uint id = regions(alias);
				if (id < (uint)-1) {
					group.regions.insert(id);
					regions.AddMetaID(group.id, groupType, id);
				}
				else
					unmatchedAliases.push_back(alias);
			}
		}
		std::cerr<<"User Defined Group ("<<group.name<<" : Type "<<MetaGroup::ConvertType(groupType)<<" ) : "<<group.regions.size()<<"\n";
		groupLookup[group.id] = groups.size();
		groups.push_back(group);
	}

	//Iterate over childrenDesc and set up any group relationships
	Utility::StringArray::iterator itr = childrenDesc.begin();
	Utility::StringArray::iterator end = childrenDesc.end();

	while (itr != end) {
		Utility::StringArray words = Utility::Split(itr->c_str(), " ");
		if (words.size() > 2) {
			uint targetID = nameToId[words[1]];
			for (uint i=2; i<words.size(); i++) {
				AddAssociation(targetID, nameToId[words[i]]);
			}
		} else {
			std::cerr<<"Incomplete CHILDREN definition. Each line should contain the target name followed by 1 or more child groups.\n";
			std::cerr<<"Ignoring line: "<<*itr<<"\n";
		}
		itr++;
	}



	return nameToId;
}


void GroupManager::ListGroupAssociation(uint idx, std::ostream& os, RegionManager& regions, Utility::IdCollection& visited, Knowledge::SnpDataset& snps, uint tabCount) {
	Group& g = groups[idx];
	//std::cerr<<tabCount<<"--\n";
	Utility::IdCollection regionsLocal = regionCollection[groupLookup[g.id]];
	Utility::IdCollection::iterator rlItr = regionsLocal.begin();
	Utility::IdCollection::iterator rlEnd = regionsLocal.end();
	Utility::IdCollection snpList;
	while (rlItr != rlEnd) {
		snpList.insert(regions[*rlItr].snps.begin(), regions[*rlItr].snps.end());
		rlItr++;
	}
	os<<std::string(tabCount, '\t')<<g.name<<":"<<g.id<<" ("<<g.groups.size()<<","<<regionCollection[groupLookup[g.id]].size()<<","<<snpList.size()<<"):\n";

	if (visited.find(g.id) == visited.end()) {
		visited.insert(g.id);

		Utility::IdCollection::iterator itr = g.regions.begin();
		Utility::IdCollection::iterator end = g.regions.end();

		while (itr != end) {
			regions[*itr++].ListGroupAssociations(os, tabCount+1, snps);
		}

		itr = g.groups.begin();
		end = g.groups.end();
		while (itr != end)
			ListGroupAssociation(*itr++, os, regions, visited, snps, tabCount+1);
	}
}

void GroupManager::ListAssociations(std::ostream& os, RegionManager& regions, Knowledge::SnpDataset& snps){
	os<<"Associations "<<name<<":"<<id<<":\n";

	Utility::IdCollection::iterator itr = root.begin();
	Utility::IdCollection::iterator end = root.end();

	if (regionCollection.size() == 0)
		BuildRegionCollections(regionCollection);
	Utility::IdCollection visited;

	while (itr != end)
		ListGroupAssociation(*itr++, os, regions, visited, snps, 1);


}

}


#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

class GroupManagerTest : public ::testing::Test {
public:
	SnpDataset dataset;
	soci::session sociDB;
	Knowledge::RegionManagerDB regions;
	GroupManager mgr;
	std::map<std::string, uint> regionAliasToID;

	GroupManagerTest() : mgr(1, MetaGroup::DiseaseIndependent, "test1", "a simple test") {}

	virtual void SetUp() {
		remove("fake-bio.db");
		sociDB.open(soci::sqlite3, "dbname=fake-bio.db timeout=2500");
		sociDB<<"CREATE TABLE regions (gene_id INTEGER UNIQUE PRIMARY KEY, primary_name VARCHAR(128) UNIQUE, chrom VARCHAR(2), description TEXT)";
		sociDB<<"CREATE TABLE region_bounds (gene_id INTEGER, population_id INTEGER, start INTEGER, end INTEGER)";
		sociDB<<"CREATE TABLE region_alias_type (region_alias_type_id INTEGER UNIQUE PRIMARY KEY, region_alias_type_desc TEXT)";
		sociDB<<"CREATE TABLE region_alias (region_alias_type_id INTEGER, alias VARCHAR(64), gene_id INTEGER, gene_count INTEGER, UNIQUE(gene_id, alias, region_alias_type_id))";

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
		mgr.AddGroup(2, "A", "first root");						// (0, 1, 2, 3, 4, 5, 6, 7)
		mgr.AddGroup(100, "a1", "child of 1");					// (0, 1, 2, 3, 4, 5)
		mgr.AddGroup(110, "a1a", "child of a1");				// 1, 2, 3
		mgr.AddGroup(111, "a1b", "child of a1");				// (0, 2, 3, 4, 5)
		mgr.AddGroup(112, "a1b1", "child of a1b");			// (0, 2, 3)
		mgr.AddGroup(114, "a1b1a", "child of a1b1");			// 0, 2
		mgr.AddGroup(115, "a1b1b", "child of a1b1");			// 0, 3
		mgr.AddGroup(113, "a1b2", "child of a1b");			// 2, 4, 5
		mgr.AddGroup(101, "a2", "another child of 1");		// (2, 5, 6, 7)
		mgr.AddGroup(116, "a2a", "child of a2");				// 2 (5, 6, 7)
		mgr.AddGroup(117, "a2a1", "child of a2a");			// 5, 6
		mgr.AddGroup(118, "a2a2", "child of a2a");			// 2, 6, 7
		mgr.AddGroup(102, "a3", "another child of 1");		// (2, 3, 4, 5, 6, 7)
		mgr.AddGroup(119, "a3a", "child of a3");				//	4, 5
		mgr.AddGroup(120, "a3b", "child of a3");				// 3, 4
		mgr.AddGroup(121, "a3c", "child of a3");				// (2, 4, 5)
		mgr.AddGroup(122, "a3c1", "child of a3c");			// 2, 4
		mgr.AddGroup(123, "a3c2", "child of a3c");			// (2, 5, 7)
		mgr.AddGroup(125, "a3c2a", "child of a3c2");			// 2, 5
		mgr.AddGroup(126, "a3c2b", "child of a3c2");			// 5, 7
		mgr.AddGroup(124, "a3c3", "child of a3c");			// 6
		mgr.AddGroup(3, "B", "second root");					// (3, 4, 5, 6, 7, 8, 9, 10)
		mgr.AddGroup(200, "b1", "child of B");					// 3, 4
		mgr.AddGroup(201, "b2", "child of B");					// 3, 10
		mgr.AddGroup(202, "b3", "child of B");					// (4, 5, 6, 7, 8, 9)
		mgr.AddGroup(205, "b3a", "child of B3");				// 5, 6
		mgr.AddGroup(206, "b3b", "child of B3");				// 5, 7, 8
		mgr.AddGroup(207, "b3c", "child of b3");				// 6, 8
		mgr.AddGroup(208, "b3d", "child of b3");				// 5, 7, 9
		mgr.AddGroup(209, "b3e", "child of b3");				// 4, 6, 7
		mgr.AddGroup(203, "b4", "child of B");					// 4, 10
		mgr.AddGroup(204, "b5", "child of b");					// 3, 4, 10
		//mgr.AddGroup(300, "asdf", "This should get deleted during squeeze");


		mgr.AddAssociation(1, 2);
		mgr.AddAssociation(1, 3);
		mgr.AddAssociation(2, 100);
		mgr.AddAssociation(100, 110);
		mgr.AddAssociation(110, 111);
		mgr.AddAssociation(110, 112);

		mgr.AddAssociation(112, 114);
		mgr.AddAssociation(112, 115);
		mgr.AddAssociation(100, 113);
		mgr.AddAssociation(2, 101);
		mgr.AddAssociation(101, 116);
		mgr.AddAssociation(116, 117);
		mgr.AddAssociation(116, 118);
		mgr.AddAssociation(2, 102);
		mgr.AddAssociation(102, 119);
		mgr.AddAssociation(102, 120);
		mgr.AddAssociation(102, 121);
		mgr.AddAssociation(121, 122);
		mgr.AddAssociation(121, 123);
		mgr.AddAssociation(123, 125);
		mgr.AddAssociation(123, 126);
		mgr.AddAssociation(121, 124);
		mgr.AddAssociation(3, 200);
		mgr.AddAssociation(3, 201);
		mgr.AddAssociation(3, 202);
		mgr.AddAssociation(202, 205);
		mgr.AddAssociation(202, 206);
		mgr.AddAssociation(202, 207);
		mgr.AddAssociation(202, 208);
		mgr.AddAssociation(202, 209);
		mgr.AddAssociation(3, 203);
		mgr.AddAssociation(3, 204);
		mgr.AddAssociation(204, 300);


		mgr.AddGeneAssociation(110, 1, regions);
		mgr.AddGeneAssociation(110, 2, regions);
		mgr.AddGeneAssociation(110, 3, regions);
		mgr.AddGeneAssociation(114, 0, regions);
		mgr.AddGeneAssociation(114, 2, regions);
		mgr.AddGeneAssociation(115, 0, regions);
		mgr.AddGeneAssociation(115, 3, regions);
		mgr.AddGeneAssociation(113, 2, regions);
		mgr.AddGeneAssociation(113, 4, regions);
		mgr.AddGeneAssociation(113, 5, regions);
		mgr.AddGeneAssociation(119, 4, regions);
		mgr.AddGeneAssociation(119, 5, regions);
		mgr.AddGeneAssociation(120, 3, regions);
		mgr.AddGeneAssociation(120, 4, regions);
		mgr.AddGeneAssociation(122, 2, regions);
		mgr.AddGeneAssociation(122, 4, regions);
		mgr.AddGeneAssociation(125, 2, regions);
		mgr.AddGeneAssociation(125, 5, regions);
		mgr.AddGeneAssociation(126, 5, regions);
		mgr.AddGeneAssociation(126, 7, regions);
		mgr.AddGeneAssociation(124, 6, regions);

		mgr.AddGeneAssociation(200, 3, regions);
		mgr.AddGeneAssociation(200, 4, regions);
		mgr.AddGeneAssociation(201, 3, regions);
		mgr.AddGeneAssociation(201, 10, regions);
		mgr.AddGeneAssociation(205, 5, regions);
		mgr.AddGeneAssociation(205, 6, regions);
		mgr.AddGeneAssociation(206, 5, regions);
		mgr.AddGeneAssociation(206, 7, regions);
		mgr.AddGeneAssociation(206, 8, regions);
		mgr.AddGeneAssociation(207, 6, regions);
		mgr.AddGeneAssociation(207, 8, regions);
		mgr.AddGeneAssociation(208, 5, regions);
		mgr.AddGeneAssociation(208, 7, regions);
		mgr.AddGeneAssociation(208, 9, regions);
		mgr.AddGeneAssociation(209, 4, regions);
		mgr.AddGeneAssociation(209, 6, regions);
		mgr.AddGeneAssociation(209, 7, regions);
		mgr.AddGeneAssociation(203, 4, regions);
		mgr.AddGeneAssociation(203, 10, regions);
		mgr.AddGeneAssociation(204, 3, regions);
		mgr.AddGeneAssociation(204, 4, regions);
		mgr.AddGeneAssociation(204, 10, regions);

		//This is to be used by the load archive test
		regionAliasToID["A1BG"]		= 0;
		regionAliasToID["A2M"]		= 1;
		regionAliasToID["A2MP1"]	= 2;
		regionAliasToID["NAT1"]		= 3;
		regionAliasToID["NAT2"]		= 4;
		regionAliasToID["AACCP"]	= 5;
		regionAliasToID["SERPINA3"]	= 6;
		regionAliasToID["AADAC"]	= 7;
		regionAliasToID["AAMP"]		= 8;
		regionAliasToID["AANAT"]	= 9;
		regionAliasToID["AARS"]		= 10;

		std::ofstream file("test-associations.txt");
		mgr.ListAssociations(file, regions, dataset);
		
	}

	virtual void TearDown() {
		remove("fake-bio.db");
		remove("test-associations.txt");
	}
};



TEST_F(GroupManagerTest, ArchiveLoading) {
	char file_contents[] = "test1 a simple test\n"
	"GROUP A first root\n"
	"GROUP a1 child of 1\n"
	"GROUP a1a child of a1\nA2M\nA2MP1\nNAT1\n"
	"GROUP a1b child of a1\n"
	"GROUP a1b1 child of a1b\n"
	"GROUP a1b1a child of a1b\nA1BG\nA2MP1\n"
	"GROUP a1b1b child of a1b\nA1BG\nNAT1\n"
	"GROUP a1b2 child of a1\nA2MP1\nNAT2\nAACCP\n"
	"GROUP a2 another child of 1\n"
	"GROUP a2a child of a2\nA2MP1\n"
	"GROUP a2a1 child of a2a\nAACCP\nSERPINA3\n"
	"GROUP a2a2 child of a2a\nA2MP1\nSERPINA3\nAADAC\n"
	"GROUP a3 another child of 1\n"
	"GROUP a3a child of a3\nNAT2\nAACCP\n"
	"GROUP a3b child of a3\nNAT1\nNAT2\n"
	"GROUP a3c child of a3\n"
	"GROUP a3c1 child of a3c\nA2MP1\nNAT2\n"
	"GROUP a3c2 child of a3c\n"
	"GROUP a3c2a child of a3c2nA2MP1\nAACCP\n"
	"GROUP a3c2b child of a3c2\nAACCP\nAADAC\n"
	"GROUP a3c3 child of a3c\nSEPINA3\n"
	"GROUP B second root\n"
	"GROUP b1 child of b\nNAT1\nNAT2\n"
	"GROUP b2 child of b\nNAT1\nAARS\n"
	"GROUP b3 child of b\n"
	"GROUP b3a child of b3\nAACCP\nSERPINA3\n"
	"GROUP b3b child of b3\nAACCP\nAADAC\nAANAT\n"
	"GROUP b3c child of b3\nSERPINA3\nAAMP\n"
	"GROUP b3d child of b3\nAACCP\nAADAC\nAANAT\n"
	"GROUP b3e child of b3\nNAT2\nSERPINA3\nAADAC\n"
	"GROUP b4 child of b\nNAT2\nAARS\n"
	"GROUP b5 child of b\nNAT1\nNAT2\nAARS\n"
	"CHILDREN A a1 a2 a3\n"
	"CHILDREN a1 a1a a1b\n"
	"CHILDREN a1a a1b a1b1\n"
	"CHILDREN a1b1 a1b1a a1b1b\n"
	"CHILDREN a2 a2a\n"
	"CHILDREN a2a a2a1 a2a2\n"
	"CHILDREN a3 a3a a3b a3c\n"
	"CHILDREN a3c a3c1 a3c2 a3c3\n"
	"CHILDREN a3c2 a3c2a a3c2b\n"
	"CHILDREN B b1 b2 b3 b4 b5\n"
	"CHILDREN b3 b3a b3b b3c b3d b3e\n";

	std::ofstream ofile("group-archive.test");
	ofile<<file_contents;
	ofile.close();

	GroupManager archived(0, MetaGroup::DiseaseIndependent, "test1", "a simple test");
	Utility::StringArray unmatched;
	std::string filename = "group-archive.test";
	std::map<std::string, uint> idLookup = archived.LoadArchive(filename.c_str(), regions, 1, unmatched);

	GeneGeneModelArchive archive;
	GroupManager::maxGeneCount = 2;
	Utility::IdCollection visited;

	std::map<uint, Utility::IdCollection> regioncollections;

	archived.BuildRegionCollections(regioncollections);
	uint id = idLookup["B"];
	uint index = archived.GetIndex(id);
	archived.GenerateGeneGeneModels(archive, regioncollections, regions, index, std::cerr, visited);

	std::stringstream ss;
	archive.WriteToArchive(ss, false);
	EXPECT_EQ("3	4	1\n3	10	1\n4	10	1\n5	6	1\n6	8	1\n", ss.str());
	remove("group-archive.test");
}


TEST_F(GroupManagerTest, GeneGeneModelGenerationGcOfTwo) {
	/*
		mgr.AddGroup(3, "B", "second root");					// (3, 4, 5, 6, 7, 8, 9, 10)
		mgr.AddGroup(200, "b1", "child of B");					// 3, 4
		mgr.AddGroup(201, "b2", "child of B");					// 3, 10
		mgr.AddGroup(202, "b3", "child of B");					// (4, 5, 6, 7, 8, 9)
		mgr.AddGroup(205, "b3a", "child of B3");				// 5, 6
		mgr.AddGroup(206, "b3b", "child of B3");				// 5, 7, 8
		mgr.AddGroup(207, "b3c", "child of b3");				// 6, 8
		mgr.AddGroup(208, "b3d", "child of b3");				// 5, 7, 9
		mgr.AddGroup(209, "b3e", "child of b3");				// 4, 6, 7
		mgr.AddGroup(203, "b4", "child of B");					// 4, 10
		mgr.AddGroup(204, "b5", "child of b");					// 3, 4, 10
	 */

	GeneGeneModelArchive archive;
	GroupManager::maxGeneCount = 2;
	Utility::IdCollection visited;
	std::map<uint, Utility::IdCollection> regioncollections;

	mgr.BuildRegionCollections(regioncollections);
	uint index = mgr.GetIndex(3);
	mgr.GenerateGeneGeneModels(archive, regioncollections, regions, index, std::cerr, visited);

	std::stringstream ss;
	archive.WriteToArchive(ss, false);
	EXPECT_EQ("3	4	1\n3	10	1\n4	10	1\n5	6	1\n6	8	1\n", ss.str());

}

TEST_F(GroupManagerTest, GeneGeneModelGenerationComplexGeneCountLarge) {
	/*
		mgr.AddGroup(3, "B", "second root");					// (3, 4, 5, 6, 7, 8, 9, 10)
		mgr.AddGroup(200, "b1", "child of B");					// 3, 4
		mgr.AddGroup(201, "b2", "child of B");					// 3, 10
		mgr.AddGroup(202, "b3", "child of B");					// (4, 5, 6, 7, 8, 9)
		mgr.AddGroup(205, "b3a", "child of B3");				// 5, 6
		mgr.AddGroup(206, "b3b", "child of B3");				// 5, 7, 8
		mgr.AddGroup(207, "b3c", "child of b3");				// 6, 8
		mgr.AddGroup(208, "b3d", "child of b3");				// 5, 7, 9
		mgr.AddGroup(209, "b3e", "child of b3");				// 4, 6, 7
		mgr.AddGroup(203, "b4", "child of B");					// 4, 10
		mgr.AddGroup(204, "b5", "child of b");					// 3, 4, 10
	 */

	GeneGeneModelArchive archive;
	GroupManager::maxGeneCount = 5;
	Utility::IdCollection visited;
	std::map<uint, Utility::IdCollection> regioncollections;

	mgr.BuildRegionCollections(regioncollections);
	uint index = mgr.GetIndex(3);
	mgr.GenerateGeneGeneModels(archive, regioncollections, regions, index, std::cerr, visited);

	std::stringstream ss;
	archive.WriteToArchive(ss, false);
	EXPECT_EQ("3	4	1\n3	10	1\n4	6	1\n4	7	1\n4	10	1\n5	6	1\n5	7	1\n5	8	1\n5	9	1\n6	7	1\n6	8	1\n7	8	1\n7	9	1\n", ss.str());

}

TEST_F(GroupManagerTest, GeneGeneModelGenerationComplexGeneCount) {
	/*
		mgr.AddGroup(202, "b3", "child of B");					// (4, 5, 6, 7, 8, 9)
		mgr.AddGroup(205, "b3a", "child of B3");				// 5, 6
		mgr.AddGroup(206, "b3b", "child of B3");				// 5, 7, 8
		mgr.AddGroup(207, "b3c", "child of b3");				// 6, 8
		mgr.AddGroup(208, "b3d", "child of b3");				// 5, 7, 9
		mgr.AddGroup(209, "b3e", "child of b3");				// 4, 6, 7
		mgr.AddGroup(203, "b4", "child of B");					// 4, 10
		mgr.AddGroup(204, "b5", "child of b");					// 3, 4, 10
	 */

	GeneGeneModelArchive archive;
	GroupManager::maxGeneCount = 5;
	Utility::IdCollection visited;
	std::map<uint, Utility::IdCollection> regioncollections;

	mgr.BuildRegionCollections(regioncollections);
	uint index = mgr.GetIndex(202);
	mgr.GenerateGeneGeneModels(archive, regioncollections, regions, index, std::cerr, visited);

	std::stringstream ss;
	archive.WriteToArchive(ss, false);
	EXPECT_EQ("4	6	1\n4	7	1\n5	6	1\n5	7	1\n5	8	1\n5	9	1\n6	7	1\n6	8	1\n7	8	1\n7	9	1\n", ss.str());

}

TEST_F(GroupManagerTest, GeneGeneModelGenerationSingleNode) {
	GeneGeneModelArchive archive;
	GroupManager::maxGeneCount = 5;
	Utility::IdCollection visited;
	//First, we'll just produce the lookup ourselves, and call the child function
	std::map<uint, Utility::IdCollection> regioncollections;
	//		mgr.AddGroup(206, "b3b", "child of B3");				// 5, 7, 8

	mgr.BuildRegionCollections(regioncollections);
	uint index = mgr.GetIndex(206);
	mgr.GenerateGeneGeneModels(archive, regioncollections, regions, index, std::cerr, visited);
	EXPECT_EQ(3, archive.Size());
	std::stringstream ss;
	archive.WriteToArchive(ss, false);
	EXPECT_EQ("5	7	1\n5	8	1\n7	8	1\n", ss.str());



}

TEST_F(GroupManagerTest, RegionCollections) {
	std::map<uint, Utility::IdCollection> regioncollections;
	//std::map<uint, Utility::IdCollection>::iterator itr;
	EXPECT_EQ("6", Utility::Join(mgr(124).regions));
	uint index = mgr.GetIndex(124);

	mgr.BuildRegionCollections(index, regioncollections);
	EXPECT_EQ(1, regioncollections.size());
	EXPECT_EQ(6, *(regioncollections[index].begin()));

	regioncollections.clear();
	index		= mgr.GetIndex(126);
	mgr.BuildRegionCollections(index, regioncollections);
	Utility::IdCollection regionids = regioncollections[index];
	EXPECT_EQ(2, regionids.size());
	EXPECT_EQ("5,7", Utility::Join<std::set<uint> >(regionids, ","));

	index = mgr.GetIndex(123);			// This is a small tree
	regioncollections.clear();
	mgr.BuildRegionCollections(index, regioncollections);
	EXPECT_EQ(3, regioncollections.size());
	regionids = regioncollections[mgr.GetIndex(126)];
	//Double check that it found the same stuff we found last time
	EXPECT_EQ(2, regionids.size());
	EXPECT_EQ("5,7", Utility::Join<std::set<uint> >(regionids, ","));

	regionids = regioncollections[mgr.GetIndex(125)];
	EXPECT_EQ(2, regionids.size());
	EXPECT_EQ("2,5", Utility::Join<std::set<uint> >(regionids, ","));

	regionids = regioncollections[index];
	EXPECT_EQ(3, regionids.size());
	EXPECT_EQ(2, *(regionids.begin()));
	EXPECT_EQ("2,5,7", Utility::Join<std::set<uint> >(regionids, ","));

	//mgr.AddGroup(2, "B", "second root");					// (3, 4, 5, 6, 7, 8, 9, 10)
	regioncollections.clear();
	index = mgr.GetIndex(3);
	mgr.BuildRegionCollections(index, regioncollections);
	EXPECT_EQ(11, regioncollections.size());

	//		mgr.AddGroup(206, "b3b", "child of B3");				// 5, 7, 8
	regionids = regioncollections[mgr.GetIndex(206)];
	EXPECT_EQ(3, regionids.size());
	EXPECT_EQ("5,7,8", Utility::Join(regionids, ","));

	//		mgr.AddGroup(207, "b3c", "child of b3");				// 6, 8
	regionids = regioncollections[mgr.GetIndex(207)];
	EXPECT_EQ(2, regionids.size());
	EXPECT_EQ("6,8", Utility::Join(regionids, ","));

	//		mgr.AddGroup(208, "b3d", "child of b3");				// 5, 7, 9
	regionids = regioncollections[mgr.GetIndex(208)];
	EXPECT_EQ(3, regionids.size());
	EXPECT_EQ("5,7,9", Utility::Join(regionids, ","));

	//mgr.AddGroup(202, "b3", "child of B");					// (4, 5, 6, 7, 8, 9)
	regionids = regioncollections[mgr.GetIndex(202)];
	EXPECT_EQ(6, regionids.size());
	EXPECT_EQ("4,5,6,7,8,9", Utility::Join(regionids, ","));

	regionids = regioncollections[index];
	EXPECT_EQ(8, regionids.size());
	EXPECT_EQ("3,4,5,6,7,8,9,10", Utility::Join<std::set<uint> >(regionids, ","));

	
	//Test the root function

	regioncollections.clear();
	mgr.BuildRegionCollections(regioncollections);
	EXPECT_EQ(32, regioncollections.size());
	//		mgr.AddGroup(1, "A", "first root");						// (0, 1, 2, 3, 4, 5, 6, 7)
	index = mgr.GetIndex(1);
	regionids = regioncollections[index];
	EXPECT_EQ(8, regionids.size());
	EXPECT_EQ("0,1,2,3,4,5,6,7", Utility::Join(regionids, ","));

	index = mgr.GetIndex(3);
	regionids = regioncollections[index];
	EXPECT_EQ(8, regionids.size());
	EXPECT_EQ("3,4,5,6,7,8,9,10", Utility::Join<std::set<uint> >(regionids, ","));

}

TEST_F(GroupManagerTest, ListAssociations) {
	std::stringstream ss;
	Utility::IdCollection visited;
	ss.str("");
	mgr.ListGroupAssociation(5, ss, regions, visited, dataset, 1);
	EXPECT_EQ("\ta1b1a:114 (0,2,10):\n\t\tA1BG (RS1000 RS1001 RS1002 RS1003 RS1004 )\n\t\tA2MP1 (RS1011 RS1012 RS1013 RS1014 RS1015 )\n", ss.str());

	ss.str("");
	visited.clear();
	mgr.ListGroupAssociation(8, ss, regions, visited, dataset, 0);
	EXPECT_EQ("a2:101 (1,0,0):\n\ta2a:116 (2,0,0):\n\t\ta2a1:117 (0,0,0):\n\t\ta2a2:118 (0,0,0):\n", ss.str());

	ss.str("");
	mgr.ListGroupAssociation(4, ss, regions, visited, dataset, 0);
	EXPECT_EQ("a1b1:112 (2,3,11):\n\ta1b1a:114 (0,2,10):\n\t\tA1BG (RS1000 RS1001 RS1002 RS1003 RS1004 )\n\t\tA2MP1 (RS1011 RS1012 RS1013 RS1014 RS1015 )\n\ta1b1b:115 (0,2,6):\n\t\tA1BG (RS1000 RS1001 RS1002 RS1003 RS1004 )\n\t\tNAT1 (RS1016 )\n", ss.str());

	ss.str("");
	mgr.ListGroupAssociation(30, ss, regions, visited, dataset, 0);
	EXPECT_EQ("b4:203 (0,2,3):\n\tNAT2 (RS1017 )\n\tAARS (RS1050 RS1051 )\n", ss.str());
}

TEST_F(GroupManagerTest, Basics) {
	EXPECT_EQ(32, mgr.Size());
	EXPECT_EQ(1, mgr.id);
	EXPECT_EQ("test1", mgr.name);
	EXPECT_EQ("a simple test", mgr.description);



	EXPECT_EQ("A", mgr[0].name);
	EXPECT_EQ("first root", mgr[0].desc);
	EXPECT_EQ(100, mgr[1].id);
	EXPECT_EQ("a1", mgr(100).name);

	EXPECT_EQ(119, mgr(119).id);
	EXPECT_EQ("a3a", mgr(119).name);

	EXPECT_EQ(3, mgr(1).groups.size());
	EXPECT_EQ(2, mgr(100).groups.size());

	EXPECT_EQ(119, mgr(119).id);
	EXPECT_EQ(3, mgr(102).groups.size());
	EXPECT_EQ(3, mgr(121).groups.size());
	EXPECT_EQ(5, mgr(3).groups.size());


	EXPECT_EQ(3, mgr(110).regions.size());
	EXPECT_EQ(3, mgr(113).regions.size());
	EXPECT_EQ(119, mgr(119).id);
	EXPECT_EQ(2, mgr(119).regions.size());
	EXPECT_EQ(2, mgr(125).regions.size());

	GroupManager mgr2(mgr);
	EXPECT_EQ(1, mgr2.id);
	EXPECT_EQ("test1", mgr2.name);
	EXPECT_EQ("a simple test", mgr2.description);

	EXPECT_EQ(3, mgr2(121).groups.size());
	EXPECT_EQ(5, mgr2(3).groups.size());


	EXPECT_EQ(3, mgr2(110).regions.size());
	EXPECT_EQ(3, mgr2(113).regions.size());

}

#endif

