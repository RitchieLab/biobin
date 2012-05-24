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

namespace Knowledge{

string __vinit[] = {"1","2","3","4","5","6","7","8","9","10",
		"11","12","13","14","15","16","17","18","19","20","21","22","X",
		"Y","XY","MT"};

// Make sure the # below matches the # of elements in the array above!!
const vector<string> Locus::_chrom_list(__vinit, __vinit + (sizeof(__vinit) / sizeof(__vinit[0])));

const string Locus::invalid_chrom("");

Locus::Locus(short chrom, uint pos, bool rare, const string& id):
		_chrom(chrom), _pos(pos), _id(id), _is_rare(rare){
	if (id.size() == 0 || _id == "."){
		createID();
	}

}

Locus::Locus(const string& chrom_str, uint pos, bool rare, const string& id):
		_pos(pos), _id(id), _is_rare(rare){
	_chrom = getChrom(chrom_str);
	if (_id.size() == 0 || _id == "."){
		createID();
	}
}

void Locus::addAllele(const string& allele, float freq){
	_alleles.insert(Allele(allele, freq, _alleles.size()));
}

float Locus::majorAlleleFreq() const{
	return _alleles.begin()->getFreq();
}

float Locus::minorAlleleFreq() const{
	if (_alleles.size() > 1){
		return (++_alleles.begin())->getFreq();
	}else{
		return 1-majorAlleleFreq();
	}
}

bool Locus::operator <(const Locus& other) const{
	return other._chrom==_chrom ? other._pos > _pos : other._chrom > _chrom;
}

bool Locus::isMinor(const string& allele) const{
	return !(allele == _alleles.begin()->getData());
}

uint Locus::distance(const Locus& other) const{
	return other._chrom == _chrom ? abs(other._pos - _pos) : -1;
}

short Locus::encodeGenotype(uint a1, uint a2) const{
	if (a1 == (uint)-1 || a2 == (uint)-1){
		return -1;
	}
	return a1 * _alleles.size() + a2;
}

pair<uint, uint> Locus::decodeGenotype(short encoded_type) const{
	pair<uint, uint> to_return = std::make_pair(-1,-1);
	if (encoded_type != -1){
		to_return.first = encoded_type / _alleles.size();
		to_return.second = encoded_type % _alleles.size();
	}
	return to_return;
}


const string& Locus::getChromStr(short chrom){
	if (chrom < 0 || (uint) chrom > _chrom_list.size()){
		return invalid_chrom;
	}else{
		return _chrom_list[chrom-1];
	}
}

short Locus::getChrom(const string& chrom_str){
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
		return -1;
	}else{
		return chr_pos - _chrom_list.begin() + 1;
	}

}

void Locus::createID(){
	std::stringstream ss;
	ss << "chr" << getChromStr(_chrom) << "-" <<_pos;
	_id = ss.str();
}

void Locus::print(ostream& o, const string& sep, bool printAlleles) const{
	o << getChromStr() << sep << _pos << sep << _id;
	if (printAlleles){
		this->printAlleles(o,sep);
	}
}

void Locus::printAlleles(ostream& o, const string& sep) const{
	set<Allele>::const_iterator itr = _alleles.begin();
	set<Allele>::const_iterator end = _alleles.end();
	if (itr != end){
		o << *itr;
	}
	while(++itr != end){
		o << sep << *itr;
	}
}

}

ostream& operator<<(ostream& o, const Knowledge::Locus& l){
	l.print(o);
	return o;
}




