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
#include "Information.h"

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
Group GroupCollection::_group_not_found(-1, "Not found", "");

GroupCollection::GroupCollection(RegionCollection& reg) : _regions(reg) {}

GroupCollection::~GroupCollection(){
	unordered_map<uint, Group*>::iterator itr = _group_map.begin();
	unordered_map<uint, Group*>::const_iterator end = _group_map.end();
	while(itr != end){
		delete itr->second;
		++itr;
	}
}

Group* GroupCollection::AddGroup(uint id, const string& name, const string& desc){
	if (_group_map.find(id) != _group_map.end()){
		return _group_map[id];
	}
	else{
		Group* new_group = new Group(id, name, desc);
		_group_map.insert(std::make_pair(id, new_group));
		return new_group;
	}
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

void GroupCollection::addAssociation(uint group_id, uint region_id){
	unordered_map<uint, Group*>::iterator itr = _group_map.find(group_id);
	if (itr != _group_map.end()){
		Region& child = _regions[region_id];
		if (_regions.isValid(child)){
			itr->second->addRegion(child);
			child.addGroup(*(itr->second));
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

void GroupCollection::Load(const vector<string>& group_names){
	unordered_set<uint> empty_set;
	Load(group_names, empty_set);
}

void GroupCollection::Load(){
	vector<string> empty_set;
	Load(empty_set);
}

void GroupCollection::LoadArchive(const string& filename, vector<string>& unmatched_aliases){

	// Open the file
	ifstream data_file(filename.c_str());
	if (!data_file.is_open()){
		std::cerr<<"WARNING: cannot find " << filename <<", ignoring.";
		return;
	}

	// Read the definition of the meta-group
	string group_def;
	getline(data_file, group_def);
	vector<string> result;
	split(result, group_def, is_any_of(" \n\t"), boost::token_compress_on);
	string src_name = result[0];
	// For now, we drop the description

	// Start reading groups
	string group_init;
	getline(data_file, group_init);
	split(result, group_init, is_any_of(" \n\t"), boost::token_compress_on);
	to_upper(result[0]);
	if (!(result[0] == "GROUP")){
		std::cerr<<"WARNING: Invalid formatting in archive file.  "
				<< filename << " will be ignored.";
		return;
	}

	uint num_groups = 0;

	uint curr_group = initGroupFromArchive(src_name, result);
	if (curr_group != (uint)-1){
		++num_groups;
	}
	string curr_line;
	while(data_file.good()){
		getline(data_file, curr_line);
		split(result, curr_line, is_any_of(" \n\t"), boost::token_compress_on);
		to_upper(result[0]);
		if (result[0] == "GROUP"){
			curr_group = initGroupFromArchive(src_name, result);
			if (curr_group != (uint)-1){
				++num_groups;
			}
		}else if (curr_group != (uint)-1){
			vector<string>::const_iterator a_itr = result.begin();
			vector<string>::const_iterator a_end = result.end();
			while(a_itr != a_end){
				RegionCollection::const_region_iterator r_itr =
						_regions.aliasBegin(*a_itr);
				RegionCollection::const_region_iterator r_end =
						_regions.aliasEnd(*a_itr);
				if (r_itr == r_end){
					unmatched_aliases.push_back(*a_itr);
				}
				while(r_itr != r_end){
					addAssociation(curr_group, (*r_itr)->getID());
					++r_itr;
				}
				++a_itr;
			}
		}
	}
}

uint GroupCollection::initGroupFromArchive(const string& src_name, const vector<string>& split_line){
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
	AddGroup(++_max_group, src_name + ":" + split_line[1], ss.str());
	return _max_group;
}

}




