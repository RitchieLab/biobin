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
using BioBin::Utility::Phenotype;

namespace BioBin{

Bin::Bin(const PopulationManager& pop_mgr, Knowledge::Group* grp, const Utility::Phenotype& pheno) :
		_is_group(true), _is_intergenic(false), _cached_case(false), _cached_control(false),
		_chrom(-1), _name(grp->getName()), _pop_mgr(pop_mgr), _pheno(pheno) {
	_member.group = grp;
}

Bin::Bin(const PopulationManager& pop_mgr, Knowledge::Region* reg, const Utility::Phenotype& pheno) :
		_is_group(false), _is_intergenic(false), _cached_case(false), _cached_control(false),
		_chrom(reg->getChrom()), _name(reg->getName()), _pop_mgr(pop_mgr), _pheno(pheno){
	_member.region = reg;
}

Bin::Bin(const PopulationManager& pop_mgr, short chrom, int bin, const Utility::Phenotype& pheno) :
		_is_group(false), _is_intergenic(true),	_cached_case(false), _cached_control(false),
		_chrom(chrom), _pop_mgr(pop_mgr), _pheno(pheno) {
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
		_cached_case(false), _cached_control(false), _chrom(other._chrom), _name(other._name),
		_extra_data(other._extra_data), _pop_mgr(other._pop_mgr), _pheno(other._pheno) {}

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

unsigned int Bin::getCaseSize() const {
	if (!_cached_case) {
		set<Knowledge::Locus*>::const_iterator itr = _variants.begin();
		int case_t = 0;
		while (itr != _variants.end()) {
			case_t += _pop_mgr.genotypeContribution(**itr,
					&(_pheno.getStatus().second));
			++itr;
		}
		_size_case_cache = case_t;
		_cached_case = true;
	}
	return _size_case_cache;
}

unsigned int Bin::getControlSize() const {
	if (!_cached_control) {
		set<Knowledge::Locus*>::const_iterator itr = _variants.begin();
		int control_t = 0;
		while (itr != _variants.end()) {
			control_t += _pop_mgr.genotypeContribution(**itr,
					&(_pheno.getStatus().first));
			++itr;
		}
		_size_control_cache = control_t;
		_cached_control = true;
	}
	return _size_control_cache;
}

} // namespace BioBin



