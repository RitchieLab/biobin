/*
 * RegionCollection.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: jrw32
 */

#include "RegionCollection.h"

#include <boost/algorithm/string.hpp>

#include "Locus.h"
#include "Region.h"
#include "Information.h"

using boost::algorithm::split;
using boost::algorithm::is_any_of;

using boost::unordered_map;
using boost::unordered_set;
using std::string;
using std::vector;
using std::set;

namespace Knowledge{

// The configured population to use here
string RegionCollection::pop_str = "NO-LD";
// The amount of gene boundary expansion (if using no-ld)
int RegionCollection::gene_expansion = 0;

vector<string> RegionCollection::c_region_files;
vector<string> RegionCollection::c_region_filter;

RegionCollection::~RegionCollection(){
	// Go through the map of id->Region* and delete all of the regions
	unordered_map<uint, Region*>::iterator itr = _region_map.begin();
	unordered_map<uint, Region*>::const_iterator end = _region_map.end();

	while (itr != end){
		delete itr->second;
		++itr;
	}

	delete _dataset;
}

Region& RegionCollection::operator[](const uint idx){
	unordered_map<uint,Region*>::iterator element = _region_map.find(idx);
	if (element != _region_map.end()){
		return *(element->second);
	}
	return region_not_found;
}

const Region& RegionCollection::operator[](const uint idx) const{
	unordered_map<uint,Region*>::const_iterator element = _region_map.find(idx);
	if (element != _region_map.end()){
		return *(element->second);
	}
	return region_not_found;
}

RegionCollection::const_region_iterator RegionCollection::aliasBegin(const string& alias) const{
	unordered_map<string,set<Region*> >:: const_iterator idx_itr = _alias_map.find(alias);
	if (idx_itr != _alias_map.end()){
		return idx_itr->second.begin();
	}
	return empty_region_set.begin();
}

RegionCollection::const_region_iterator RegionCollection::aliasEnd(const string& alias) const{
	unordered_map<string,set<Region*> >:: const_iterator idx_itr = _alias_map.find(alias);
	if (idx_itr != _alias_map.end()){
		return idx_itr->second.end();
	}
	return empty_region_set.end();
}

RegionCollection::const_region_iterator RegionCollection::locusBegin(const Locus* loc) const{
	unordered_map<const Locus*, set<Region*> >::const_iterator itr = _locus_map.find(loc);
	if(itr != _locus_map.end()){
		return (*itr).second.begin();
	}
	return empty_region_set.begin();
}

RegionCollection::const_region_iterator RegionCollection::locusEnd(const Locus* loc) const{
	unordered_map<const Locus*, set<Region*> >::const_iterator itr = _locus_map.find(loc);
	if(itr != _locus_map.end()){
		return (*itr).second.end();
	}
	return empty_region_set.end();
}

void RegionCollection::insertRegion(Region& region){
	if(_region_map.find(region.getID()) == _region_map.end()){
		_region_map.insert(std::make_pair(region.getID(), &region));
	}
	Region::const_alias_iterator alias_itr = region.aliasBegin();

	_alias_map[region.getName()].insert(&region);
	while(alias_itr != region.aliasEnd()){
		_alias_map[*alias_itr++].insert(&region);
	}
}

Knowledge::Region* RegionCollection::AddRegion(const string& name, uint id, short chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const string& aliases) {

	if (_region_map.find(id) == _region_map.end()) {
		// Insert the region into the map
		Region& new_region = *(new Region(name, id, chrom, effStart, effStop,
				trueStart, trueStop));

		// WARNING: inserting a region with an identical ID will result in a memory leak!
		_region_map.insert(std::make_pair(id, &new_region));
		new_region.addAliases(aliases);

		// Add all aliases, including the canonical name
		_alias_map[name].insert(&new_region);
		vector<string> aliasList;
		split(aliasList, aliases, is_any_of(","));
		vector<string>::const_iterator itr = aliasList.begin();
		vector<string>::const_iterator end = aliasList.end();
		while (itr != end) {
			_alias_map[*itr++].insert(&new_region);
		}

		return &new_region;
	} else {
		return (*_region_map.find(id)).second;
	}

}

Knowledge::Region* RegionCollection::AddRegion(const string& name, uint id, short chrom, uint start, uint stop, const string& aliases) {
	return AddRegion(name, id, chrom, start, stop, start, stop, aliases);
}

uint RegionCollection::Load(){
	unordered_set<uint> empty_set;
	return this->Load(empty_set);
}

uint RegionCollection::Load(const unordered_set<uint>& ids){
	vector<string> empty_list;
	return this->Load(ids,empty_list);
}

uint RegionCollection::Load(const vector<string>& alias_list){
	unordered_set<uint> empty_set;
	return this->Load(empty_set, alias_list);
}

bool RegionCollection::isValid(const Region& other) const{
	return (region_not_found.getID() != other.getID());
}


};


