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
#include <boost/program_options.hpp>

#include "Group.h"

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
	enum Ambiguity_ENUM { STRICT, RESOLVABLE, PERMISSIVE };

	class AmbiguityModel{
	public:
		AmbiguityModel() : _data(RESOLVABLE){}
		AmbiguityModel(const Ambiguity_ENUM& d):_data(d){}
		operator const char*() const{
			switch(_data){
			case GroupCollection::STRICT:
				return "strict";
			case GroupCollection::RESOLVABLE:
				return "resolvable";
			case GroupCollection::PERMISSIVE:
				return "permissive";
			default:
				return "unknown";
			}
		}
		operator int() const{return _data;}

	private:
		Ambiguity_ENUM _data;
	};


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
	Group& operator[](unsigned int id);
	/*!
	 * Accessing a given group by ID.  Again, will not create an object if one
	 * is not found, but will return the special _group_not_found Group.
	 *
	 * \param id The database ID of the group to be returned
	 *
	 * \return The group with the given ID, or the special "group not found"
	 */
	const Group& operator[](unsigned int id) const;

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
	virtual void Load(const std::vector<std::string>& groupNames,
			const boost::unordered_set<unsigned int>& ids) = 0;

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
	void Load(const std::vector<std::string>& groupNames);
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
	void LoadArchive(const std::string& filename,
			std::vector<std::string>& unmatched_aliases);

	/*!
	 * Returns the number of Groups in the collection.
	 *
	 * \return The size of the collection.
	 */
	unsigned int size() { return _group_map.size(); }

	/*!
	 * Adds a parent/child relationship between two groups.  This will also add
	 * the reltationship to the group objects themselves.
	 *
	 * \param parent_id The Database ID of the parent group
	 * \param child_id The Database ID of the child group
	 */
	void addRelationship(unsigned int parent_id, unsigned int child_id);
	/*!
	 * Adds an association between a group and a region.  This will also add the
	 * relationships to the objects themselves.  Note that the region must be
	 * able to be found in the given RegionCollection object.
	 *
	 * \param group_id The Database ID of the containing group
	 * \param region_id The database ID of the region
	 * \param regions The collection of regions where the given ID can be found
	 */
	void addAssociation(unsigned int group_id, unsigned int region_id);

	//! A configuration value of group names to include
	static std::vector<std::string> c_group_names;
	//! A configuration value of group IDs to include
	static boost::unordered_set<unsigned int> c_id_list;
	//! The ambiguity setting to use
	static AmbiguityModel c_ambiguity;

protected:
	/*!
	 * Adds a single group to the collection and adds it to the internal map.
	 *
	 * \param id The database ID of the group (must be unique)
	 * \param name The name of the group
	 * \param desc A description of the group's functionality
	 */
	Group* AddGroup(unsigned int id, const std::string& name, const std::string& desc="");

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
	virtual unsigned int getMaxGroup() = 0;

	//! The maximum group number (for loading from archive)
	unsigned int _max_group;
	//! Mapping of ids -> Groups
	boost::unordered_map<unsigned int, Group*> _group_map;
	//! mapping of parent->child relationships
	boost::unordered_map<unsigned int, boost::unordered_set<unsigned int> > _group_relationships;
	//! mapping of group->region relationships
	boost::unordered_map<unsigned int, boost::unordered_set<Region*> > _group_associations;

	RegionCollection& _regions;

	Information* _info;

private:
	static Group _group_not_found;

	/*!
	 * A helper function that initializes a group from an archive file.
	 */
	unsigned int initGroupFromArchive(const std::string& src_name, const std::vector<std::string>& split_line);
};


//namespace std{
std::istream& operator>>(std::istream& in, GroupCollection::AmbiguityModel& model_out);
std::ostream& operator<<(std::ostream& o, const GroupCollection::AmbiguityModel& m);
//}


}
//void validate(boost::any& v, const std::vector<std::string>& values, Knowledge::GroupCollection::AmbiguityModel* m, int);


#endif /* KNOWLEDGE_GROUPCOLLECTION_H */
