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
map<const string*, unsigned short, Allele::str_cmp> Allele::s_string_map;


Allele::Allele(const string& data, float freq, unsigned short pos):	_freq(freq), _pos(pos) {
	map<const string*, unsigned short, Allele::str_cmp>::const_iterator itr = s_string_map.find(&data);
	if(itr != s_string_map.end()){
		_data_idx = (*itr).second;
	} else {
		_data_idx = (*s_string_map.insert(s_string_map.end(), std::make_pair(&(*s_string_pool.insert(s_string_pool.end(),data)),s_string_pool.size()))).second;
	}

}

// I want to oder them from largest to smallest alele frequency
bool Allele::operator<(const Allele& other) const{
	return (_freq==other._freq ? s_string_pool[_data_idx] < s_string_pool[other._data_idx] : _freq > other._freq);
}

/*
bool Allele::operator>(const Allele& other) const{
	return (other < *this);
}
bool Allele::operator==(const Allele& other) const{
	return (!(*this < other) && !(other < *this));
}
*/

void Allele::print(ostream& o, const string& sep) const{
	o << s_string_pool[_data_idx] << sep << _freq;
}

}

ostream& operator<<(ostream& o, const Knowledge::Allele& l){
	l.print(o);
	return o;
}

