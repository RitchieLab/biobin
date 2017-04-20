/* 
 * File:   converter.cpp
 * Author: torstees
 * 
 * Created on May 6, 2011, 1:25 PM
 */

#include <iostream>

#include "Converter.h"
#include "Chain.h"

using std::pair;
using std::set;
using std::map;
using std::string;
using std::make_pair;

using Knowledge::Locus;

namespace Knowledge {
namespace Liftover {

const float Converter::MIN_MAPPING_FRACTION = 0.95;

const pair <short, pair<int,int> > Converter::FAILED_REGION =
		make_pair(-1,make_pair(0,0));

Converter::Converter(const string& orig) :
		_origBuild(orig){}

Converter::~Converter() {
	map<short, set<Chain*> >::iterator itr = _chains.begin();

	while (itr != _chains.end()){
		set<Chain*>::iterator s_itr = (*itr).second.begin();
		set<Chain*>::iterator s_end = (*itr).second.end();

		while(s_itr != s_end){
			delete *s_itr;
			++s_itr;
		}

		(*itr).second.clear();

		++itr;
	}
	_chains.clear();
}

pair<short, pair<int, int> > Converter::convertRegion(short chrom, int start, int end) const {
	pair<short, pair<int, int> > ret_val = FAILED_REGION;

	map<short, set<Chain*> >::const_iterator map_itr = _chains.find(chrom);
	if (map_itr != _chains.end()){
		set<Chain*>::const_iterator itr = (*map_itr).second.begin();
		while(ret_val == FAILED_REGION && itr != (*map_itr).second.end()){
			if((*itr)->overlaps(start, end)){
				pair<int, int> new_reg = (*itr)->convertRegion(start, end, MIN_MAPPING_FRACTION);

				if(new_reg.first != new_reg.second){
					ret_val = make_pair((*itr)->getNewChrom(), new_reg);
				}
			}
			++itr;
		}
	}

	return ret_val;

}
Locus* Converter::convertLocus(const Locus& old_loc) const {
	pair<short, pair<int, int> > new_region = convertRegion(old_loc.getChrom(),
			old_loc.getPos(), old_loc.getPos() + 1);

	if (new_region == FAILED_REGION) {
		return 0;
	} else {
		Locus* converted = new Locus(new_region.first, new_region.second.first,
				old_loc.getID());

		return converted;

	}
}

} // namespace Lifover
} // namespace Knowledge


