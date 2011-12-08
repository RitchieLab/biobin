/*
 * GroupCollection.h
 *
 *  Created on: Nov 14, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_GROUPCOLLECTION_H
#define KNOWLEDGE_GROUPCOLLECTION_H

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "Group.h"

using boost::unordered_map;
using boost::unordered_set;
using std::string;
using std::vector;

namespace Knowledge{

class RegionCollection;

class GroupCollection{

public:
	/**
	 * Creates a collection of groups from a single source, identified by the
	 * id and name.
	 *
	 * @param id The database ID for the group collection (category)
	 * @param name The name of the group collection (category)
	 */
	GroupCollection(uint id, const string& name="");
	/**
	 * Destroys the collection.  This must go through the map and destroy all
	 * of the created groups.
	 */
	virtual ~GroupCollection();

	/**
	 * Adds a single group to the collection and adds it to the internal map.
	 *
	 * @param id The database ID of the group (must be unique)
	 * @param name The name of the group
	 * @param desc A description of the group's functionality
	 */
	void addGroup(uint id, const string& name, const string& desc="");
	/**
	 * Adds a parent/child relationship between two groups.  This will also add
	 * the reltationship to the group objects themselves.
	 *
	 * @param parent_id The Database ID of the parent group
	 * @param child_id The Database ID of the child group
	 */
	void addRelationship(uint parent_id, uint child_id);
	/**
	 * Adds an association between a group and a region.  This will also add the
	 * relationships to the objects themselves.  Note that the region must be
	 * able to be found in the given RegionCollection object.
	 *
	 * @param group_id The Database ID of the containing group
	 * @param region_id The database ID of the region
	 * @param regions The collection of regions where the given ID can be found
	 */
	void addAssociation(uint group_id, uint region_id, RegionCollection& regions);

	/**
	 * Determines if a given group is the special "invalid" group
	 *
	 * @param other The group to check for validity
	 */
	bool isValid(const Group& other) const;

	/**
	 * Accessing a given group by ID.  NOTE: unlike the standard operator[] for
	 * maps, this will not create an object if one is not found.  Instead, it
	 * will return the special _group_not_found Group.
	 *
	 * @param id The databse ID of the group to be returned
	 */
	Group& operator[](uint id);
	/**
	 * Accessing a given group by ID.  Again, will not create an object if one
	 * is not found, but will return the special _group_not_found Group.
	 *
	 * @param id Te database ID of the group to be returned
	 */
	const Group& operator[](uint id) const;

	/**
	 * Load the given ids into memory.  The RegionCollection is simply a pass
	 * through to the addAssociation method.  NOTE: this function must be
	 * defined in an inherited class.
	 */
	virtual uint Load(RegionCollection& regions, const vector<string>& groupNames,
			const unordered_set<uint>& ids) = 0;

	virtual uint Load(RegionCollection& regions, const vector<string>& groupNames);
	/**
	 * Load all IDs into memory.  Calls Load(regions, ids) with an empty set of
	 * ids
	 */
	virtual uint Load(RegionCollection& regions);

	uint LoadArchive(RegionCollection& regions, const string& filename,
			vector<string>& unmatched_aliases);

	virtual uint getMaxGroup() = 0;

	uint size() { return _group_map.size(); }

protected:
	// The ID of the type (source) of groups
	uint _id;
	// The name of the set of groups
	string _name;
	// The maximum group number (for loading from archive)
	uint _max_group;
	// Mapping of ids -> Groups
	unordered_map<uint, Group*> _group_map;
	// mapping of parent->child relationships
	unordered_map<uint, unordered_set<uint> > _group_relationships;
	// mapping of group->region relationships
	unordered_map<uint, unordered_set<Region*> > _group_associations;

private:
	Group _group_not_found;

	uint initGroupFromArchive(const vector<string>& split_line);
};

}


#endif /* KNOWLEDGE_GROUPCOLLECTION_H */
