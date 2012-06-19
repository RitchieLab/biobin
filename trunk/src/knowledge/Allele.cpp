/*
 * Allele.cpp
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#include "Allele.h"

namespace Knowledge{

set<string> Allele::s_string_pool;

// I want to oder them from largest to smallest alele frequency
bool Allele::operator<(const Allele& other) const{
	return (_freq==other._freq ? *_data < *(other._data) : _freq > other._freq);
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
	o << *_data << sep << _freq;
}

}

ostream& operator<<(ostream& o, const Knowledge::Allele& l){
	l.print(o);
	return o;
}

