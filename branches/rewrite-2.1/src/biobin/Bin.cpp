/*
 * Bin.cpp
 *
 *  Created on: Dec 6, 2011
 *      Author: jrw32
 */

#include "Bin.h"

#include <sstream>

#include "PopulationManager.h"
#include "binmanager.h"

using std::stringstream;
using std::set;
using std::list;

namespace BioBin{

Bin::Bin(const PopulationManager& pop_mgr, Knowledge::Group* grp) :
		_is_group(true), _is_intergenic(false), _cached(false), _chrom(-1),
		_name(grp->getName()), _pop_mgr(pop_mgr) {
	_member.group = grp;
}

Bin::Bin(const PopulationManager& pop_mgr, Knowledge::Region* reg) :
		_is_group(false), _is_intergenic(false), _cached(false), _chrom(reg->getChrom()),
		_name(reg->getName()), _pop_mgr(pop_mgr){
	_member.region = reg;
}

Bin::Bin(const PopulationManager& pop_mgr, short chrom, int bin) :
		_is_group(false), _is_intergenic(true),	_cached(false), _chrom(chrom),
		_pop_mgr(pop_mgr) {
	stringstream ss;
	ss << "chr" << Knowledge::Locus::getChromStr(chrom) << ":"
			<< bin*BinManager::IntergenicBinStep << "K-"
			<< bin*BinManager::IntergenicBinStep + BinManager::IntergenicBinWidth
			<< "K";
	_name = ss.str();
	_member.bin_no = bin;
}

Bin::Bin(const Bin& other) : _member(other._member),
		_is_group(other._is_group), _is_intergenic(other._is_intergenic),
		_cached(false), _chrom(other._chrom), _name(other._name),
		_extra_data(other._extra_data), _pop_mgr(other._pop_mgr) {}

bool Bin::operator<(const Bin& other) const{
	bool ret_val = false;
	if(_is_group){
		if(other._is_group){
			ret_val = _member.group == other._member.group ?
					_name < other._name :
					*_member.group < *other._member.group;
		}else{
			ret_val = true;
		}
	}else if (_is_intergenic){
		if(other._is_intergenic){
			ret_val = (_chrom == other._chrom ?
					(_member.bin_no == other._member.bin_no ?
							_name < other._name :
							_member.bin_no < other._member.bin_no) :
					_chrom < other._chrom);
		}else{
			ret_val = false;
		}
	}else{
		// Must be a region
		if (other._is_group){
			ret_val = false;
		}else if (other._is_intergenic){
			ret_val = true;
		}else{
			ret_val = (_member.region == other._member.region) ?
					_name < other._name :
					*_member.region < *other._member.region;
		}
	}

	return ret_val;

}

unsigned int Bin::getSize() const{
	int ret_val = 0;
	if(!_cached){
		list<Knowledge::Locus*>::const_iterator itr = _variants.begin();
		while(itr != _variants.end()){
			ret_val += _pop_mgr.genotypeContribution(**itr);
			++itr;
		}
		_size_cache = ret_val;
		_cached = true;
	}else{
		ret_val = _size_cache;
	}
	return ret_val;
}

} // namespace BioBin



