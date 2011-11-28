/*
 * Group.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#include "Group.h"

#include "Region.h"

namespace Knowledge{

Group::Group(uint id, const string& name, const string& desc) :
		_id(id), _name(name), _description(desc) {}

template <class T_iter>
void Group::addChildren(T_iter& begin, const T_iter& end){
	while(begin != end){
		_children.insert(*begin);
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



