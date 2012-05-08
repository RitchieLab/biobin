/*
 * GroupCollection.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#include "GroupCollection.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "RegionCollection.h"

using std::ifstream;
using std::getline;
using std::stringstream;

using boost::algorithm::split;
using boost::algorithm::join;
using boost::algorithm::is_any_of;
using boost::to_upper;

namespace Knowledge{

vector<string> GroupCollection::c_group_names;
unordered_set<uint> GroupCollection::c_id_list;
uint GroupCollection::_max_group = 0;

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
		parent_itr->second->addChild(*(child_itr->second));
		//child_itr->second->addParent(*(parent_itr->second));
	}
}

void GroupCollection::addAssociation(uint group_id, uint region_id,
		RegionCollection& regions){
	unordered_map<uint, Group*>::iterator itr = _group_map.find(group_id);
	if (itr != _group_map.end()){
		Region& child = regions[region_id];
		if (regions.isValid(child)){
			itr->second->addRegion(child);
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

uint GroupCollection::Load(RegionCollection& regions,
		const vector<string>& group_names){
	unordered_set<uint> empty_set;
	return Load(regions, group_names, empty_set);
}

uint GroupCollection::Load(RegionCollection& regions){
	vector<string> empty_set;
	return Load(regions, empty_set);
}

uint GroupCollection::LoadArchive(RegionCollection& regions,
		const string& filename, vector<string>& unmatched_aliases){

	// Open the file
	ifstream data_file(filename.c_str());
	if (!data_file.is_open()){
		std::cerr<<"WARNING: cannot find " << filename <<", ignoring.";
		return 1;
	}

	// Read the definition of the meta-group
	string group_def;
	getline(data_file, group_def);
	vector<string> result;
	split(result, group_def, is_any_of(" \n\t"), boost::token_compress_on);
	_name = result[0];
	// For now, we drop the description

	// Start reading groups
	string group_init;
	getline(data_file, group_init);
	split(result, group_init, is_any_of(" \n\t"), boost::token_compress_on);
	to_upper(result[0]);
	if (!(result[0] == "GROUP")){
		std::cerr<<"WARNING: Invalid formatting in archive file.  "
				<< filename << " will be ignored.";
		return 2;
	}

	uint num_groups = 0;

	uint curr_group = initGroupFromArchive(result);
	if (curr_group != (uint)-1){
		++num_groups;
	}
	string curr_line;
	while(data_file.good()){
		getline(data_file, curr_line);
		split(result, curr_line, is_any_of(" \n\t"), boost::token_compress_on);
		to_upper(result[0]);
		if (result[0] == "GROUP"){
			curr_group = initGroupFromArchive(result);
			if (curr_group != (uint)-1){
				++num_groups;
			}
		}else if (curr_group != (uint)-1){
			vector<string>::const_iterator a_itr = result.begin();
			vector<string>::const_iterator a_end = result.end();
			while(a_itr != a_end){
				RegionCollection::const_region_iterator r_itr =
						regions.aliasBegin(*a_itr);
				RegionCollection::const_region_iterator r_end =
						regions.aliasEnd(*a_itr);
				if (r_itr == r_end){
					unmatched_aliases.push_back(*a_itr);
				}
				while(r_itr != r_end){
					addAssociation(curr_group, (*r_itr)->getID(), regions);
					++r_itr;
				}
				++a_itr;
			}
		}
	}

	return 0;
}

uint GroupCollection::initGroupFromArchive(const vector<string>& split_line){
	if (split_line.size() < 2){
		std::cerr << "Invalid formatting: no group name given, ignoring.";
		return -1;
	}
	// I know there's at least 2 elements to the vector now, so this is safe.
	vector<string>::const_iterator itr = ++(++(split_line.begin()));
	vector<string>::const_iterator end = split_line.end();
	stringstream ss;
	if (itr != end){
		ss << *itr;
	}
	while (++itr != end){
		ss << " " << *itr;
	}
	addGroup(++_max_group, split_line[1], ss.str());
	return _max_group;
}

}




