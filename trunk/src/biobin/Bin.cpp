/*
 * Bin.cpp
 *
 *  Created on: Dec 6, 2011
 *      Author: jrw32
 */

#include "Bin.h"

#include <sstream>

#include "knowledge/Group.h"
#include "knowledge/Region.h"
#include "knowledge/Locus.h"

using std::stringstream;

namespace BioBin{

Bin::Bin(Knowledge::Group* grp) : _is_group(true), _is_intergenic(false),
		_chrom(-1), _name(grp->getName()) {
	_member.group = grp;
}

Bin::Bin(Knowledge::Region* reg) : _is_group(false), _is_intergenic(false),
		_chrom(reg->getChrom()), _name(reg->getName()){
	_member.region = reg;
}

Bin::Bin(short chrom, int bin) :_is_group(false), _is_intergenic(true),
		_chrom(chrom) {
	stringstream ss;
	ss << "chr" << Knowledge::Locus::getChromStr(chrom) << "-" << bin;
	_name = ss.str();
	_member.bin_no = bin;
}

bool Bin::operator<(const Bin& other) const{
	bool ret_val = false;
	if(_is_group){
		if(other._is_group){
			ret_val = _member.group < other._member.group;
		}else{
			ret_val = true;
		}
	}else if (_is_intergenic){
		if(other._is_intergenic){
			ret_val = (_chrom == other._chrom ?
					_member.bin_no < other._member.bin_no :
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
			ret_val = _member.region < other._member.region;
		}
	}

	return ret_val;

}

} // namespace BioBin



