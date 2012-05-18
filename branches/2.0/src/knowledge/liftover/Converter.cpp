/* 
 * File:   converter.cpp
 * Author: torstees
 * 
 * Created on May 6, 2011, 1:25 PM
 */

#include "Converter.h"
#include "Chain.h"

using std::make_pair;

namespace Knowledge {
namespace Liftover {

const float Converter::MIN_MAPPING_FRACTION = 0.95;

const pair <short, pair<int,int> > Converter::FAILED_REGION =
		make_pair(-1,make_pair(-1,-1));

Converter::Converter(const string& orig) :
		_origBuild(orig){}

Converter::~Converter() {
	map<short, set<Chain*> >::iterator itr = _chains.begin();
	map<short, set<Chain*> >::iterator end = _chains.end();

	while (itr != end){
		set<Chain*>::iterator s_itr = (*itr).second.begin();
		set<Chain*>::iterator s_end = (*itr).second.end();

		while(s_itr != s_end){
			delete *s_itr;
			++s_itr;
		}

		++itr;
	}
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

} // namespace Lifover
} // namespace Knowledge


