/*
 * GroupCollection.h
 *
 *  Created on: Nov 14, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_GROUPCOLLECTION_H
#define KNOWLEDGE_GROUPCOLLECTION_H

#include <string>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "group.h"

using boost::unordered_map;
using boost::unordered_set;
using std::string;

class RegionCollection;

namespace Knowledge{

class GroupCollection{

public:
	GroupCollection(uint id);
	virtual ~GroupCollection();

	void addGroup(uint id, const char *name, const char *desc="");
	void addRelationship(uint parent_id, uint child_id);
	void addAssociation(uint group_id, uint region_id, RegionCollection& regions);

	Group& operator[](uint id);
	Group& operator[](const string& name);
	const Group& operator[](uint id) const;
	const Group& operator[](const string& name) const;

protected:
	uint _id;
	string _name;
	unordered_map<uint, Group> _group_map;
	unordered_map<uint, unordered_set<uint> > _group_relationships;
	unordered_map<uint, Region&> _group_associations;
	unordered_map<string, uint> _group_alias;

private:
	Group _group_not_found;
};

}


#endif /* KNOWLEDGE_GROUPCOLLECTION_H */
