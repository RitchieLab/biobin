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

#include "Region.h"

using std::string;
using std::set;


namespace Knowledge{

/*!
 * \brief A representation of a group of Regions.
 * This class represents a group of Region classes.  Regions within the group
 * are directly related; however, Groups can have relationships among themselves
 * in a hierarchical nature.  To maintain this nature, a Group can have both
 * parents and children.  For ease of implementation and traversal, Groups
 * maintain both forward links to children and backwards links to parents.
 * While the terms "child" and "parent" are used, the actual structure of the
 * GroupCollection may not be hierarchical or even acyclic, but the relationships
 * are directed.
 */
class Group{

public:
	/*!
	 * A typedef to hide the implementation of the storage of Regions
	 */
	typedef set<Region*>::const_iterator const_region_iterator;
	/*!
	 * A typedef to hide the implementation of the storage of Groups
	 */
	typedef set<Group*>::const_iterator const_group_iterator;

	/*!
	 * \brief Constructs a group.
	 * Creates a Group object that contains a numeric ID, which is unique, as
	 * well as a name, which should also be unique.  Additionally, a description
	 * can be given, but this is currently not used anywhere.
	 *
	 * \param id Unique ID of the group.
	 * \param name The canonical name of the group.
	 * \param desc A description of the group.
	 */
	Group(uint id, const string& name, const string& desc = "");

	/*!
	 * \brief Adds children groups.
	 * Adds a set of children groups through the use of an iterator.  The
	 * iterator must be incrementable, and when derefenced return a Group*.
	 *
	 * \tparam begin The start of the collection of Groups
	 * \tparam end The end of the collection of Groups to add
	 */
	template <class T_iter>
	void addChildren(T_iter& begin, const T_iter& end);
	/*!
	 * \brief Adds a single child Group.
	 * Adds a Group as a child of this Group.
	 *
	 * \param child The Group to add as a child.
	 */
	void addChild(Group& child) {_children.insert(&child); child.addParent(*this);}

	/*!
	 * \brief Adds a collection of Regions.
	 * Associates a collection of Region objects with this group using
	 * iterators.  The iterators must be incrementable and return Region*
	 * objects when dereferenced.
	 *
	 * \tparam begin The beginning iterator of the collection of Regions
	 * \tparam end The ending iterator of the collection of Regions
	 */
	template <class T_iter>
	void addRegions(T_iter& begin, const T_iter& end);
	/*!
	 * \brief Adds a single Region
	 * Associates a single Region object with this Group.
	 *
	 * \param region The Region to associate with this group.
	 */
	void addRegion(Region& region) {_regions.insert(&region);}

	/*!
	 * Gets the starting iterator of the collection of Regions associated with
	 * this Group.
	 *
	 * \return The starting iterator of all the Regions.
	 */
	const_region_iterator regionBegin() const {return _regions.begin();}
	/*!
	 * Gets the ending iterator of the collection of Regions associated with
	 * this Group.
	 *
	 * \return The ending iterator of all the Regions.
	 */
	const_region_iterator regionEnd() const {return _regions.end();}
	/*!
	 * Gets the beginning iterator for all parents of the current Group.
	 *
	 * \return The beginning iterator of all parents of the group.
	 */
	const_group_iterator parentBegin() const {return _parents.begin();}
	/*!
	 * Gets the ending iterator for all parents of the current Group.
	 *
	 * \return The ending iterator of parents of the group.
	 */
	const_group_iterator parentEnd() const {return _parents.end();}
	/*!
	 * Gets the starting iterator of all children of the current Group.
	 *
	 * \return The starting iterator of children of this group.
	 */
	const_group_iterator childBegin() const {return _children.begin();}
	/*!
	 * Gets the ending iterator of all children of the Group.
	 *
	 * \return The ending iterator of children of the group.
	 */
	const_group_iterator childEnd() const {return _children.end();}

	/*!
	 * Returns the unique ID of the group.
	 *
	 * \return The ID of the group.
	 */
	uint getID() const {return _id;}

	/*!
	 * Returns the canonical name of the group.
	 *
	 * \return The name of the group.
	 */
	const string& getName() const {return _name;}
	/*!
	 * Returns the description of the group (useful for display).
	 *
	 * \return The description of the group.
	 */
	const string& getDescription() const {return _description;}

	/*!
	 * \brief Comparison operator for STL groups.
	 * This
	 *
	 * \param other Other Group to compare this to
	 *
	 * \return A boolean comparing the group IDs
	 */
	bool operator<(const Group& other) const {return _id < other._id;}


private:
	// Prohibit copying and assignment
	Group(const Group& other);
	Group& operator=(const Group& other);

	/*!
	 * Adds parents by use of iterators.  Probably will not be used, because
	 * we now only add children, which is then responsible for adding the
	 * backward links, so addParent is the likelier candidate
	 */
	template <class T_iter>
	void addParents(T_iter& begin, const T_iter& end);

	/*!
	 * Adds a Group as a parent of the current Group.  Should only be called
	 * within this class.
	 */
	void addParent(Group& parent) {_parents.insert(&parent);}

	uint _id;
	string _name;
	string _description;

	set<Region*> _regions;
	set<Group*> _children;
	set<Group*> _parents;

};

template <class T_iter>
void Group::addChildren(T_iter& begin, const T_iter& end){
	while(begin != end){
		_children.insert(*begin);
		(*begin)->addParent(*this);
		++begin;
	}
}

template <class T_iter>
void Group::addParents(T_iter& begin, const T_iter& end){
	while(begin != end){
		_parents.insert(*begin);
		++begin;
	}
}

template <class T_iter>
void Group::addRegions(T_iter& begin, const T_iter& end){
	while(begin != end){
		_regions.insert(*begin);
		++begin;
	}
}

}

#endif /* KNOWLEDGE_GROUP_H */
