/*
 * Group.h
 *
 *  Created on: Nov 18, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_GROUP_H
#define KNOWLEDGE_GROUP_H

#include <string>
#include <set>

using std::string;
using std::set;


namespace Knowledge{

class Region;

class Group{

public:
	typedef set<Region*>::const_iterator const_region_iterator;
	typedef set<Group*>::const_iterator const_group_iterator;

	Group(uint id, const string& name, const string& desc = "");

	template <class T_iter>
	void addChildren(T_iter& begin, const T_iter& end);

	template <class T_iter>
	void addParents(T_iter& begin, const T_iter& end);

	template <class T_iter>
	void addRegions(T_iter& begin, const T_iter& end);

	void addParent(Group* parent) {_parents.insert(parent);}
	void addChild(Group* child) {_children.insert(child);}
	void addRegion(Region* region) {_regions.insert(region);}

	const_region_iterator regionBegin() const {return _regions.begin();}
	const_region_iterator regionEnd() const {return _regions.end();}
	const_group_iterator parentBegin() const {return _parents.begin();}
	const_group_iterator parentEnd() const {return _parents.end();}
	const_group_iterator childBegin() const {return _children.begin();}
	const_group_iterator childEnd() const {return _children.end();}

	uint getID() const {return _id;}
	const string& getName() const {return _name;}
	const string& getDescription() const {return _description;}


private:
	// Prohibit copying and assignment
	Group(const Group& other);
	Group& operator=(const Group& other);

	uint _id;
	string _name;
	string _description;

	set<Region*> _regions;
	set<Group*> _children;
	set<Group*> _parents;

};

}


#endif /* KNOWLEDGE_GROUP_H */
