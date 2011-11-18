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

class Region;

namespace Knowledge{

class Group{

public:
	Group(uint id, const string& name, string desc = "");

	template <class T_iter>
	void addChildren(T_iter& begin, const T_iter& end);

	template <class T_iter>
	void addParents(T_iter& begin, const T_iter& end);

	template <class T_iter>
	void addRegions(T_iter& begin, const T_iter& end);

private:
	uint _id;
	string _name;
	string _description;

	set<Region*> _regions;
	set<Group*> _children;
	set<Group*> _parents;

};

}


#endif /* KNOWLEDGE_GROUP_H */
