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

namespace Knowledge{

// The configured population to use here
string RegionCollection::pop_str = "NO-LD";
// The amount of gene boundary expansion (if using no-ld)
int RegionCollection::gene_expansion = 0;

RegionCollection::~RegionCollection(){
	// Go through the map of id->Region* and delete all of the regions
	unordered_map<uint, Region*>::iterator itr = _region_map.begin();
	unordered_map<uint, Region*>::const_iterator end = _region_map.end();

	while (itr != end){
		delete itr->second;
		++itr;
	}

	delete _dataset;
	delete _info;
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

RegionCollection::const_region_iterator RegionCollection::positionBegin(short chrom, uint pos) const{
	unordered_map<short,interval_map<uint, set<Region*> > >::const_iterator idx_itr =
			_region_bounds.find(chrom);
	if (idx_itr != _region_bounds.end()){
		interval_map<uint, set<Region*> >::const_iterator interval_itr =
				(*idx_itr).second.find(pos);
		if (interval_itr != (*idx_itr).second.end()){
			return (*interval_itr).second.begin();
		}
	}
	return empty_region_set.begin();
}

RegionCollection::const_region_iterator RegionCollection::positionEnd(short chrom, uint pos) const{
	unordered_map<short,interval_map<uint, set<Region*> > >::const_iterator idx_itr =
			_region_bounds.find(chrom);
	if (idx_itr != _region_bounds.end()){
		interval_map<uint, set<Region*> >::const_iterator interval_itr =
				(*idx_itr).second.find(pos);
		if (interval_itr != (*idx_itr).second.end()){
			return (*interval_itr).second.end();
		}
	}
	return empty_region_set.end();
}

Knowledge::Region* RegionCollection::AddRegion(const string& name, uint id, short chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const string& aliases) {

	if (_region_map.find(id) == _region_map.end()) {
		// Insert the region into the map
		Region& new_region = *(new Region(name, id, chrom, effStart, effStop,
				trueStart, trueStop));

		// WARNING: inserting a region with an identical ID will result in a memory leak!
		_region_map.insert(std::make_pair(id, &new_region));
		new_region.addAliases(aliases);

		set<Region*> new_set;
		new_set.insert(&new_region);

		// Make it available to the interval map
		_region_bounds[chrom].add(
				std::make_pair(interval<uint>::closed(effStart, effStop),
						new_set));

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


