/*
 * Allele.cpp
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#include "Allele.h"

using std::string;
using std::deque;
using std::map;

using std::ostream;

namespace Knowledge{

deque<string> Allele::s_string_pool;
map<const string, unsigned short> Allele::s_string_map;


Allele::Allele(const string& data, unsigned short pos):	_pos(pos) {
	map<const string, unsigned short>::const_iterator itr = s_string_map.find(data);
	if(itr != s_string_map.end()){
		_data_idx = (*itr).second;
	} else {
		// This probably deserves a comment:
		// If we are here, we have determined that the string we have been passed
		// is unique, so we need to insert it into the static deque of strings,
		// then, insert the string, index pair into the map so we can find
		// the string when constructing later Alleles.

		// Note that this line is essentially the equivalent to:
		_data_idx = s_string_pool.size();
		s_string_map.insert(std::make_pair(data,s_string_pool.size()));
		s_string_pool.push_back(data);
		//_data_idx = (*s_string_map.insert(s_string_map.end(), std::make_pair(
		//		*s_string_pool.insert(s_string_pool.end(), data),
		//		s_string_pool.size()))).second;
	}

}

bool Allele::operator<(const Allele& other) const{
	return (_data_idx < other._data_idx);
}

void Allele::print(ostream& o, const string& sep) const{
	o << s_string_pool[_data_idx];
}

}

ostream& operator<<(ostream& o, const Knowledge::Allele& l){
	l.print(o);
	return o;
}

