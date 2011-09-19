/* 
 * File:   groupmanager.h
 * Author: torstees
 *
 * This is responsible for aggregating groups under a single "metagroup".
 * This class will be responsible for deciding who 
 *
 * Created on March 11, 2011, 12:13 PM
 */

#ifndef GROUPMANAGER_H
#define	GROUPMANAGER_H


#include "regionmanager.h"
#include "snpdataset.h"
#include "group.h"
#include "genegenemodelarchive.h"

namespace Knowledge {



class GroupManager {
public:
	GroupManager(uint id, MetaGroup::Type groupType, const char *name = "??", 
			const char *desc = "") : id(id), name(name), description(desc),
			groupType(groupType) {}
	GroupManager(const GroupManager& orig);
	virtual ~GroupManager() {}

	/**
	 * Generate gene gene models based on maxGenecount
    * @param geneArchive
    * @param regions
    * @param maxGeneCount
    * @param os
    */
	void GenerateGeneGeneModels(GeneGeneModelArchive& geneArchive,
		RegionManager& regions, std::ostream& os);
	void GenerateGeneGeneModels(GeneGeneModelArchive& geneArchive, 
		std::map<uint, Utility::IdCollection>& regionLookup,
		RegionManager& regions,
		uint idx, 
		std::ostream& os,
		Utility::IdCollection& visited);


	/**
	 * Traverses tree populating groupIDs with indexes for any leaf nodes
	 * @param snpCountThreshold The minimum number of snps present before stepping down to the next level
    * @param regions
    * @param snps
    * @param groupIDs
    * /
	void CollectLeafIndexes(uint snpCountThreshold, 
		RegionManager& regions, 
		Utility::IdCollection& groupIDs,
		std::map<uint, Utility::IdCollection>& regionLookup);

	/ * *
	 * Recursive function that actually populates the groupIDs structure with leaf indexes
    * @param snpCountThresh
    * @param regions
    * @param snps
    * @param groupIDs
    * @param regionLookup
    * @param visited
    * /
	void CollectLeafIndexes(uint snpCountThresh, 
		uint index, 
		RegionManager& regions, 
		Utility::IdCollection& groupIDs, 
		std::map<uint, Utility::IdCollection>& regionLookup, 
		Utility::IdCollection& visited);
	*/
	void AddGroup(uint id, const char *name, const char *desc = "");
	
	/**
	 * Associate a group with another group
    * @param groupID
    * @param childGroupID
    */
	void AddAssociation(uint groupID, uint childGroupID);
	
	/**
	 * Associate a SNP gene (via index) with a group (via ID)
    * @param groupID
    * @param regionIndex
    * @param regions
    */
	void AddGeneAssociation(uint groupID, uint regionIndex, RegionManager& regions);
	//void Squeeze();																/// Purge all empty (unassociated) groups (do only after adding groups and genes)

	
	/**
	 * Lookup via index
    * @param idx
    * @return 
    */
	Group& operator[](uint idx);
	
	/**
	 * Lookup via ID
    * @param id
    * @return 
    */
	Group& operator()(uint id);
	
	/**
	 * Lookup via name
    * @param name
    * @return 
    */
	Group& operator[](std::string& name);

	uint Size();

	uint GetIndex(uint id);

	/**
	 * Load data from file
    * @param filename
    * @param regions
    * @param id
    * @param unmatchedAliases
    * @return 
    */
	std::map<std::string, uint> LoadArchive(const char *filename,
		RegionManager& regions, uint id,
		Utility::StringArray& unmatchedAliases);
	
	/**
	 * Builds up the collection of region IDs that are associated to a given group
	 * These structures are important because the system can simply query
	 * the map for all IDs associated with a given group-and not have to step
	 * into children to perform the collection. It also allows tree traversal
	 * to quickly determine if the local node is sufficient for model generation
	 * or if the we should step into the children.
    * @param regions
    */
	void BuildRegionCollections(std::map<uint, Utility::IdCollection >& regions);
	void BuildRegionCollections(uint idx, std::map<uint, Utility::IdCollection >& regions);

	/**
	 * Generate a tree hierarchy representing the contents of the group
    * @param os
    */
	void ListAssociations(std::ostream& os, RegionManager& regions,
		Knowledge::SnpDataset& snps);
	
	/**
	 * Writes a group relationship report to the stream
    * @param idx
    * @param os
    * @param regions
    * @param visited		-- Used to short circuit any loops
    * @param snps			-- snp lookup in case we are going to write out contained RS 
    * @param tabCount	-- For pretty printing
    */
	void ListGroupAssociation(uint idx, std::ostream& os, RegionManager& regions,
		Utility::IdCollection& visited, Knowledge::SnpDataset& snps, uint tabCount);
	void PrintLookup() {
		std::map<uint, uint>::iterator itr = groupLookup.begin();
		std::map<uint, uint>::iterator end = groupLookup.end();
		uint i=0;

		while (itr != end) {
			std::cerr<<++i<<" this["<<itr->first<<"] = \t"<<itr->second<<"\n";
			itr++;
		}
	}
	
	Utility::IdCollection Root();												///< Return the root entries

	static uint maxGeneCount;													///< Max gene count in order to determine whether models are to be generated or not

	uint id;																			///< DB ID associated with the group type
	std::string name;																///< Name associated in the database
	std::string description;													///< Description from the database
protected:
	std::vector<Group> groups;													///< vector of group indexes
	Utility::IdCollection root;												///< This represents the root level of the tree
	std::map<uint, uint> groupLookup;										///< Converts ID to current index. -1 indicates the group is no longer present
	MetaGroup::Type groupType;													///< Disease independent, Disease dependent, etc...
	std::map<uint, Utility::IdCollection > regionCollection;			///< Calculated on the fly
	std::map<std::string, uint> idxByName;									///< Lookup for index
};


inline
GroupManager::GroupManager(const GroupManager& other) : id(other.id), name(other.name), description(other.description), groupType(other.groupType) {
	groups			= other.groups;
	root				= other.root;
	groupLookup		= other.groupLookup;
	regionCollection = other.regionCollection;
	idxByName		= other.idxByName;
}

inline
uint GroupManager::GetIndex(uint id) {
	return groupLookup[id];
}

inline
uint GroupManager::Size() {
	return groups.size();
}
}

#endif	/* GROUPMANAGER_H */

