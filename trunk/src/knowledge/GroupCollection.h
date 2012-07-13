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
class Information;

/*!
 * \brief A collection of groups, all from a single source.
 * This class represents a collection of Group objects from a single source in
 * the database.
 */
class GroupCollection{

public:
	/*!
	 * Creates a collection of groups from a single source, identified by the
	 * id and name.
	 *
	 * \param id The database ID for the group collection (category)
	 * \param name The name of the group collection (category)
	 */
	GroupCollection(RegionCollection& reg);
	/*!
	 * Destroys the collection.  This must go through the map and destroy all
	 * of the created groups.
	 */
	virtual ~GroupCollection();

	/*!
	 * Determines if a given group is the special "invalid" group
	 *
	 * \param other The group to check for validity
	 *
	 * \return A boolean that is false <==> other.ID() == _group_not_found.ID()
	 */
	bool isValid(const Group& other) const;

	/*!
	 * Accessing a given group by ID.  NOTE: unlike the standard operator[] for
	 * maps, this will not create an object if one is not found.  Instead, it
	 * will return the special _group_not_found Group.
	 *
	 * \param id The databse ID of the group to be returned
	 *
	 * \return The group with the given ID, or the special "group not found"
	 */
	Group& operator[](uint id);
	/*!
	 * Accessing a given group by ID.  Again, will not create an object if one
	 * is not found, but will return the special _group_not_found Group.
	 *
	 * \param id The database ID of the group to be returned
	 *
	 * \return The group with the given ID, or the special "group not found"
	 */
	const Group& operator[](uint id) const;

	/*!
	 * \brief Loads the groups into memory by name and ID.
	 * Load the given ids into memory.  The RegionCollection is simply a pass
	 * through to the addAssociation method.  NOTE: this function must be
	 * defined in an inherited class.
	 *
	 * \param groupNames A list of names of groups to load
	 * \param ids A list of IDs of groups to load
	 *
	 * \return 0 if successful, anything else, otherwise.
	 */
	virtual void Load(const vector<string>& groupNames,
			const unordered_set<uint>& ids) = 0;

	/*!
	 * \brief Loads groups into memory by name.
	 * Loads groups with the given names into memory,  Again, the
	 * RegionCollection is a pass through to addAssociation, and we assume that
	 * all regions from the database have been loaded.
	 *
	 * \param groupNames The names of the groups to load.
	 *
	 * \return 0 if successful, anything else, otherwise.
	 */
	void Load(const vector<string>& groupNames);
	/*!
	 * \brief Load all IDs into memory.
	 * Load all IDs into memory.  Calls Load(regions, ids) with an empty set of
	 * ids.
	 *
	 * \param regions
	 *
	 * \return 0 if successful, anything else, otherwise.
	 */
	void Load();

	/*!
	 * \brief Loads Groups from a file.
	 * Loads a collection of groups from a file with specific format.  Will
	 * number groups starting at the getMaxGroup() group ID and auto-increment.
	 * This method is typically used to load custom groups into memory.
	 *
	 * \param filename The name of the file containing the group definition
	 * \param[out] unmatched_aliases A list of the aliases of Regions that were
	 * not found in the given RegionCollection
	 *
	 * \return 0 if successful, anything else, otherwise.
	 */
	void LoadArchive(const string& filename,
			vector<string>& unmatched_aliases);

	/*!
	 * Returns the number of Groups in the collection.
	 *
	 * \return The size of the collection.
	 */
	uint size() { return _group_map.size(); }

	/*!
	 * Adds a parent/child relationship between two groups.  This will also add
	 * the reltationship to the group objects themselves.
	 *
	 * \param parent_id The Database ID of the parent group
	 * \param child_id The Database ID of the child group
	 */
	void addRelationship(uint parent_id, uint child_id);
	/*!
	 * Adds an association between a group and a region.  This will also add the
	 * relationships to the objects themselves.  Note that the region must be
	 * able to be found in the given RegionCollection object.
	 *
	 * \param group_id The Database ID of the containing group
	 * \param region_id The database ID of the region
	 * \param regions The collection of regions where the given ID can be found
	 */
	void addAssociation(uint group_id, uint region_id);

	//! A configuration value of group names to include
	static vector<string> c_group_names;
	//! A configuration value of group IDs to include
	static unordered_set<uint> c_id_list;

protected:
	/*!
	 * Adds a single group to the collection and adds it to the internal map.
	 *
	 * \param id The database ID of the group (must be unique)
	 * \param name The name of the group
	 * \param desc A description of the group's functionality
	 */
	Group* AddGroup(uint id, const string& name, const string& desc="");

	/*!
	 * \brief Returns the maximum ID contained within the database.
	 * Returns the maximum ID of any group within the data store.  If any
	 * custom groups have been added, it will return the largest ID of the
	 * groups that are added.
	 *
	 * NOTE:  This method is NOT thread-safe, and if two parallel data stores
	 * are used, this could cause problems.
	 *
	 * \return The maximum ID within the database (including custom-loaded groups)
	 */
	virtual uint getMaxGroup() = 0;

	//! The maximum group number (for loading from archive)
	uint _max_group;
	//! Mapping of ids -> Groups
	unordered_map<uint, Group*> _group_map;
	//! mapping of parent->child relationships
	unordered_map<uint, unordered_set<uint> > _group_relationships;
	//! mapping of group->region relationships
	unordered_map<uint, unordered_set<Region*> > _group_associations;

	static const vector<string> child_types;
	RegionCollection& _regions;

	Information* _info;

private:
	static Group _group_not_found;

	/*!
	 * A helper function that initializes a group from an archive file.
	 */
	uint initGroupFromArchive(const string& src_name, const vector<string>& split_line);
};

}


#endif /* KNOWLEDGE_GROUPCOLLECTION_H */
