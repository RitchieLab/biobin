/*
 * Region.cpp
 *
 *  Created on: Nov 18, 2011
 *      Author: jrw32
 */

#include "Region.h"
#include "Locus.h"
#include "Group.h"

#include <sstream>

using std::stringstream;
using std::deque;
using std::string;
using std::set;
using std::map;
using std::vector;


namespace Knowledge{

Region::Region(const string& name, uint id) : _name(name), _id(id) {}

Region::Region(const string& name, uint id, short chrom, uint start, uint end) :
		_name(name), _id(id){
	addDefaultBoundary(chrom, start, end);
}

Region::Region(const string& name, uint id, short chrom,
			uint eff_start, uint eff_end, uint true_start, uint true_end) :
			_name(name), _id(id) {
	addPopulationBoundary(chrom, eff_start, eff_end);
	addDefaultBoundary(chrom, true_start, true_end);
}

void Region::addLocus(const Locus& locus){
	_locus_set.insert(&locus);
}

void Region::addGroup(Group& container){
	_group_set.insert(&container);
}

string Region::getAliasString(const string& sep) const{
	deque<string>::const_iterator itr = _aliases.begin();

	stringstream ss;

	if (itr != _aliases.end()){
		ss<<*itr;
	}

	while (++itr != _aliases.end()){
		ss<<sep<<*itr;
	}

	return ss.str();
}

bool Region::containsLocus(const Locus& other) const{
	return (_locus_set.count(&other) != 0);
}

void Region::addAliases(const string& aliases, const string& sep){

	int init_pos = 0;
	int end_pos = aliases.find(sep);
	do {
		_aliases.push_back(aliases.substr(init_pos, end_pos-init_pos));
		init_pos = end_pos + sep.size();
		end_pos = aliases.find(sep, init_pos);
	} while((int) end_pos != (int) string::npos);
}

bool Region::operator<(const Region& other) const{
	// order by chromosome, then start, then end position, then by id
	return (_def_bounds == other._def_bounds) ?
			_id < other._id : _def_bounds < other._def_bounds;
}
}

