/*
 * Locus.cpp
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#include "Locus.h"
#include <algorithm>
#include <sstream>
#include <boost/algorithm/string.hpp>

using std::new_handler;
using std::set_new_handler;
using std::sort;
using std::find;
using std::iter_swap;
using std::string;
using std::set;
using std::vector;
using std::pair;
using std::ostream;
using boost::pool;

namespace Knowledge{

string __vinit[] = {"1","2","3","4","5","6","7","8","9","10",
		"11","12","13","14","15","16","17","18","19","20","21","22","X",
		"Y","XY","MT"};

// Make sure the # below matches the # of elements in the array above!!
const vector<string> Locus::_chrom_list(__vinit, __vinit + (sizeof(__vinit) / sizeof(__vinit[0])));

const string Locus::invalid_chrom("");

pool<> Locus::s_locus_pool(sizeof(Locus));

void* Locus::operator new(size_t size) {
	if (size != sizeof(Locus)) {
		return ::operator new(size);
	}
	void * ret_val = 0;
	while (!ret_val) {
		ret_val = s_locus_pool.malloc();

		// Do rituals for out-of-memory conditions here!
		if (!ret_val) {
			new_handler gh = set_new_handler(0);
			set_new_handler(gh);

			if (gh) {
				(*gh)();
			} else {
				throw std::bad_alloc();
			}
		}
	}

	return ret_val;
}

void Locus::operator delete(void* deadObj, size_t size) {
	if (deadObj == 0) {
		return;
	}
	if (size != sizeof(Locus)) {
		::operator delete(deadObj);
		return;
	}

	s_locus_pool.free(deadObj);
}

Locus::Locus(short chrom, uint pos, const string& id, const string& ref):
		_chrom(chrom), _pos(pos), _id(id){
	if (id.size() == 0 || _id == "."){
		createID(ref);
	}

}

Locus::Locus(const string& chrom_str, uint pos, const string& id, const string& ref):
		_chrom(getChrom(chrom_str)), _pos(pos), _id(id){
	if (_id.size() == 0 || _id == "."){
		createID(ref);
	}
}

bool Locus::operator <(const Locus& other) const{
	return _chrom == other._chrom ?
			_pos < other._pos :
			_chrom < other._chrom;
}

unsigned int Locus::distance(const Locus& other) const{
	return _chrom == other._chrom ?	abs(_pos - other._pos) : -1;
}

const string& Locus::getChromStr(unsigned short chrom){
	if (chrom < 0 || (uint) chrom > _chrom_list.size()){
		return invalid_chrom;
	}else{
		return _chrom_list[chrom-1];
	}
}

unsigned short Locus::getChrom(const string& chrom_str){
	string eval_str = boost::to_upper_copy(chrom_str);

	// remove the 'CHR' at the beginning, should it exist
	if (eval_str.substr(0,3) == "CHR"){
		eval_str.erase(0,3);
	}

	// Catch a couple of special cases here...
	if (eval_str == "X|Y"){
		eval_str = "XY";
	}else if (eval_str == "M"){
		eval_str = "MT";
	}

	//Now, try to find in our vector, and give the position
	vector<string>::const_iterator chr_pos = find(_chrom_list.begin(),
			_chrom_list.end(), eval_str);

	// If this is true, we did not find an exact match, return -1
	if (chr_pos == _chrom_list.end()){
		return UNKNOWN_CHROM;
	}else{
		return chr_pos - _chrom_list.begin() + 1;
	}

}

void Locus::createID(const string& ref){
	std::stringstream ss;
	ss << "chr" << getChromStr(_chrom) << "-" << _pos;
	if(!ref.empty()){
		ss << "-" << ref;
	}
	if(!ref.empty()){
		ss << "-" << ref;
	}
	_id = ss.str();
}

void Locus::print(ostream& o, const string& sep) const{
	o << getChromStr() << sep << _pos << sep << _id;
}
/*
void Locus::printAlleles(ostream& o, const string& sep) const{
	vector<Allele>::const_iterator itr = _alleles.begin();
	if (itr != _alleles.end()){
		o << *itr;
	}
	while(++itr != _alleles.end()){
		o << sep << *itr;
	}
}
*/
}

ostream& operator<<(ostream& o, const Knowledge::Locus& l){
	l.print(o);
	return o;
}




