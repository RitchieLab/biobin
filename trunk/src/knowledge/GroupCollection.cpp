/*
 * GroupCollection.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#include "GroupCollection.h"
#include "RegionCollection.h"

namespace Knowledge{

GroupCollection::GroupCollection(uint id, const string& name) :
	_id(id), _name(name), _group_not_found(-1,"Not Found") {}

GroupCollection::~GroupCollection(){
	unordered_map<uint, Group*>::iterator itr = _group_map.begin();
	unordered_map<uint, Group*>::const_iterator end = _group_map.end();
	while(itr != end){
		delete itr->second;
		++itr;
	}
}

void GroupCollection::addGroup(uint id, const string& name, const string& desc){
	Group* new_group = new Group(id, name, desc);
	_group_map.insert(std::make_pair(id, new_group));
}

void GroupCollection::addRelationship(uint parent_id, uint child_id){
	unordered_map<uint, Group*>::iterator parent_itr = _group_map.find(parent_id);
	unordered_map<uint, Group*>::iterator child_itr = _group_map.find(child_id);
	unordered_map<uint, Group*>::const_iterator end = _group_map.end();
	if (parent_itr != end && child_itr != end){
		_group_relationships[parent_id].insert(child_id);
		parent_itr->second->addChild(child_itr->second);
		child_itr->second->addParent(parent_itr->second);
	}
}

void GroupCollection::addAssociation(uint group_id, uint region_id,
		RegionCollection& regions){
	unordered_map<uint, Group*>::iterator itr = _group_map.find(group_id);
	if (itr != _group_map.end()){
		Region& child = regions[region_id];
		if (regions.isValid(child)){
			itr->second->addRegion(&child);
			child.addGroup(_id, *(itr->second));
			_group_associations[group_id].insert(&child);
		}
	}
}

Group& GroupCollection::operator[](uint id){
	unordered_map<uint, Group*>::iterator itr = _group_map.find(id);
	if (itr != _group_map.end()){
		return *(itr->second);
	}
	return _group_not_found;
}

const Group& GroupCollection::operator[](uint id) const {
	unordered_map<uint, Group*>::const_iterator itr = _group_map.find(id);
	if (itr != _group_map.end()){
		return *(itr->second);
	}
	return _group_not_found;
}

bool GroupCollection::isValid(const Group& other) const{
	return other.getID() != _group_not_found.getID();
}

uint GroupCollection::Load(RegionCollection& regions){
	unordered_set<uint> empty_set;
	return Load(regions, empty_set);
}

}




