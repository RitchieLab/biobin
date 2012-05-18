/*
 * Chain.cpp
 *
 *  Created on: May 18, 2012
 *      Author: jrw32
 */

#include "Chain.h"
#include <algorithm>

using std::max;
using std::min;

namespace Knowledge{

namespace Liftover{

pair<int, int> Chain::convertRegion(int start, int end, float minMappingFrac) const{
	pair<int, int> ret_val;
	if (!overlaps(start, end)){
		ret_val = make_pair(NOT_INTERSECTING,NOT_INTERSECTING);
	}else{
		Segment first_seg;
		Segment end_seg;
		// Get all segments that intersect this region
		set<Segment>::const_iterator seg_itr = _data.lower_bound(start);
		bool is_first = true;
		int total_size = 0;
		while(seg_itr != _data.end() && (*seg_itr).getOldEnd() <= end){
			if((*seg_itr).getOldEnd() >= start){
				if(is_first){
					first_seg = *seg_itr;
					is_first = false;
				}
				end_seg = *seg_itr;
				total_size += (*seg_itr).size();
			}
			++seg_itr;
		}

		if (is_first){
			ret_val = make_pair(DELETED,DELETED);
		}else{
			total_size -= end_seg.size();
			int front_diff = max(0, min(first_seg.size(), start-first_seg.getOldStart()));
			int end_diff = max(0, min(end_seg.size(), end-end_seg.getOldStart()));

			total_size += end_diff - front_diff;

			int new_s = _is_fwd ? (first_seg.getNewStart() + front_diff) : (end_seg.getNewStart() - end_diff);
			int new_e = _is_fwd ? (end_seg.getNewStart() + end_diff) : (first_seg.getNewStart() - front_diff);

			if((end - start) / static_cast<float>(total_size) < minMappingFrac){
				ret_val = make_pair(PARTIALLY_DELETED,PARTIALLY_DELETED);
			}else{
				ret_val = make_pair(new_s, new_e);
			}
		}
	}

	return ret_val;
}

void Chain::addSegment(int old_s, int old_e, int new_s){
	_data.insert(_data.end(),Segment(old_s,old_e,new_s));
}

}
}
