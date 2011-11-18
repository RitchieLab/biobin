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

namespace Knowledge{

Region::Region(const string& name, uint id, short chrom, uint start, uint end) :
		_name(name), _chrom(chrom), _id(id), _true_start(start), _true_end(end),
		_eff_start(start), _eff_end(end) {}

Region::Region(const string& name, uint id, short chrom,
			uint eff_start, uint eff_end, uint true_start, uint true_end) :
			_name(name), _chrom(chrom), _id(id), _true_start(true_start),
			_true_end(true_end), _eff_start(eff_start), _eff_end(eff_end){}

void Region::addLocus(const Locus& locus){
	_locus_map[locus.getID()] = &locus;
}

template <class T_iter>
void Region::addLoci(T_iter& begin, const T_iter& end){
	while (begin != end){
		addLocus(*(*begin));
		++begin;
	}
}

string Region::getAliasString(const string& sep) const{
	set<string>::const_iterator end = _aliases.end();
	set<string>::const_iterator itr = _aliases.begin();

	stringstream ss;

	if (itr != end){
		ss<<*itr;
	}

	while (++itr != end){
		ss<<sep<<*itr;
	}

	return ss.str();
}

bool Region::containsLocus(const Locus& other) const{
	return (_locus_map.find(other.getID()) != _locus_map.end());
}

void Region::addAliases(const string& aliases, const string& sep){

	int init_pos = 0;
	int end_pos = aliases.find(sep);
	do {
		_aliases.insert(aliases.substr(init_pos, end_pos-init_pos));
		init_pos = end_pos + sep.size();
		end_pos = aliases.find(sep, init_pos);
	} while(end_pos != string::npos);
}


}

