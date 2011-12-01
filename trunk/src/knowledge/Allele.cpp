/*
 * Allele.cpp
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#include "Allele.h"

namespace Knowledge{

bool Allele::operator<(const Allele& other) const{
	return (_freq==other._freq ? _data < other._data : _freq < other._freq);
}
bool Allele::operator>(const Allele& other) const{
	return !(*this < other || *this == other);
}
bool Allele::operator==(const Allele& other) const{
	return (_freq==other._freq && _data == other._data);
}

void Allele::print(ostream& o, const string& sep) const{
	o << _data << sep << _freq;
}

}



