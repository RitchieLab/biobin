/*
 * RegionCollection.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: jrw32
 */

#include "RegionCollection.h"


namespace Knowledge{

Region& RegionCollection::operator[](const uint idx){
	unordered_map<uint,Region>::iterator element = region_map.find(idx);
	if (element != region_map.end()){
		return element->second;
	}
	return region_not_found;
}

Region& RegionCollection::operator[](const string& alias){
	unordered_map<string,uint>::iterator idx_itr = alias_map.find(alias);
	if (idx_itr != alias_map.end()){
		unordered_map<uint,Region>::iterator element = region_map.find(idx_itr->second);
		if (element != region_map.end()){
			return element->second;
		}
	}
	return region_not_found;
}


const Region& RegionCollection::operator[](const uint idx) const{
	unordered_map<uint,Region>::const_iterator element = region_map.find(idx);
	if (element != region_map.end()){
		return element->second;
	}
	return region_not_found;
}

const Region& RegionCollection::operator[](const string& alias) const{
	unordered_map<string,uint>:: const_iterator idx_itr = alias_map.find(alias);
	if (idx_itr != alias_map.end()){
		unordered_map<uint,Region>::const_iterator element = region_map.find(idx_itr->second);
		if (element != region_map.end()){
			return element->second;
		}
	}
	return region_not_found;
}

/**
 * Erases all regions which are not associated with any SNPs in the dataset
 * This should be called only by AssociateSNPs
 */
void RegionCollection::Squeeze(){
	unordered_map<uint, Region>::iterator end = region_map.end();
	for (unordered_map<uint, Region>::iterator itr = region_map.begin(); itr!=end;itr++){
		if (itr->second.SnpCount() == 0){
			// Erase all of the aliases associated with this region
			for (Utility::StringArray::iterator alias_itr=itr->second.aliases.begin(); alias_itr != itr->second.aliases.begin(); alias_itr++){
				alias_map.erase(*alias_itr);
			}
			// And get rid of the region
			region_map.erase(itr);
		}
	}
}

void RegionCollection::AddRegion(const char *name, uint id, char chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const char* aliases) {
	assert(chrom > 0);
	// Insert the region into the map
	region_map.insert(std::make_pair(id, Region(name, id, chrom, effStart, effStop, trueStart, trueStop)));
	region_map[id].AddAliases(aliases);

	// Add all aliases, including the canonical name
	alias_map[name] = id;
	Utility::StringArray aliasList = Utility::Split(aliases, ",");
	Utility::StringArray::const_iterator itr = aliasList.begin();
	Utility::StringArray::const_iterator end = aliasList.end();
	while (itr != end){
		alias_map[*itr++] = id;
	}

}

void RegionCollection::AddRegion(const char *name, uint id, char chrom, uint start, uint stop, const char *aliases) {
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

bool RegionCollection::isValid(const Region& other){
	return (region_not_found.id == other.id);
}
};


