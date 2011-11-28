/*
 * RegionCollection.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: jrw32
 */

#include "RegionCollection.h"

#include "Locus.h"
#include "Region.h"

#include "utility/strings.h"


namespace Knowledge{

RegionCollection::~RegionCollection(){
	// Go through the map of id->Region* and delete all of the regions
	unordered_map<uint, Region*>::iterator itr = _region_map.begin();
	unordered_map<uint, Region*>::const_iterator end = _region_map.end();

	while (itr != end){
		delete itr->second;
		++itr;
	}
}

Region& RegionCollection::operator[](const uint idx){
	unordered_map<uint,Region*>::iterator element = _region_map.find(idx);
	if (element != _region_map.end()){
		return *(element->second);
	}
	return region_not_found;
}

/*
const set<Region&>& RegionCollection::operator[](const string& alias){
	unordered_map<string,uint>::iterator idx_itr = _alias_map.find(alias);
	if (idx_itr != _alias_map.end()){
		unordered_map<uint,Region>::iterator element = _region_map.find(idx_itr->second);
		if (element != _region_map.end()){
			return element->second;
		}
	}
	return region_not_found;
}
*/


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

/**
 * Erases all regions which are not associated with any SNPs in the dataset
 * This should be called only by AssociateSNPs
 */
void RegionCollection::Squeeze(){
	unordered_map<uint, Region*>::iterator end = _region_map.end();
	for (unordered_map<uint, Region*>::iterator itr = _region_map.begin(); itr!=end;itr++){
		if (itr->second->locusCount() == 0){
			// Erase all of the aliases associated with this region
			for (Region::const_alias_iterator alias_itr=itr->second->aliasBegin(); alias_itr != itr->second->aliasEnd(); alias_itr++){
				_alias_map[*alias_itr].erase(itr->second);
			}
			// And get rid of the region
			_region_map.erase(itr);
			delete itr->second;
		}
	}
}

void RegionCollection::AddRegion(const string& name, uint id, short chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const string& aliases) {
	assert(chrom > 0);
	// Insert the region into the map

	Region& new_region = *(new Region(name, id, chrom, effStart, effStop, trueStart, trueStop));

	// WARNING: inserting a region with an identical ID will result in a memory leak!
	_region_map.insert(std::make_pair(id, &new_region));
	new_region.addAliases(aliases);

	set<Region*> new_set;
	new_set.insert(&new_region);

	// Make it available to the interval map
	_region_bounds[chrom].add(
			std::make_pair(interval<uint>::closed(effStart, effStop), new_set));

	// Add all aliases, including the canonical name
	_alias_map[name].insert(&new_region);
	Utility::StringArray aliasList = Utility::Split(aliases.c_str(), ",");
	Utility::StringArray::const_iterator itr = aliasList.begin();
	Utility::StringArray::const_iterator end = aliasList.end();
	while (itr != end){
		_alias_map[*itr++].insert(&new_region);
	}

}

void RegionCollection::AddRegion(const string& name, uint id, short chrom, uint start, uint stop, const string& aliases) {
	AddRegion(name, id, chrom, start, stop, start, stop, aliases);
}

uint RegionCollection::Load(){
	return this->Load(0);
}

uint RegionCollection::Load(const uint pop_id){
	unordered_set<uint> empty_set;
	return this->Load(pop_id,empty_set);
}

uint RegionCollection::Load(const uint pop_id, const unordered_set<uint>& ids){
	vector<string> empty_list;
	return this->Load(pop_id,ids,empty_list);
}

uint RegionCollection::Load(const uint pop_id, const vector<string>& alias_list){
	unordered_set<uint> empty_set;
	return this->Load(pop_id, empty_set, alias_list);
}

bool RegionCollection::isValid(const Region& other) const{
	return (region_not_found.getID() == other.getID());
}

template <class T_iter>
void RegionCollection::associateLoci(T_iter& begin, const T_iter& end){
	while (begin != end){
		interval_map<uint, set<Region*> >& chrom_map =
				_region_bounds[begin->getChrom()];
		interval_map<uint, set<Region*> >::const_iterator region_itr =
				chrom_map.find(begin->getPos());
		if (region_itr != chrom_map.end()){
			set<Region*>::iterator set_itr=region_itr->second.begin();
			set<Region*>::const_iterator set_end=region_itr->second.end();
			while (set_itr != set_end){
				(*set_itr)->addLocus(*begin);
				++set_itr;
			}
		}
		++begin;
	}
}
};


