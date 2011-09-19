/* 
 * File:   regioncontainer.h
 * Author: torstees
 *
 * Created on June 15, 2011, 10:33 AM
 */

#ifndef REGIONCONTAINER_H
#define	REGIONCONTAINER_H

#include "utility/types.h"
#include <boost/icl/split_interval_map.hpp>

namespace Knowledge {
class RegionContainer {
public:
	struct Region {
		Region() : lBound(0), rBound(0), index(0) { }
		Region(uint lBound, uint rBound, int index) : lBound(lBound), rBound(rBound), index(index) { }
		Region(const Region& other) : lBound(other.lBound), rBound(other.rBound), index(other.index) {}

		bool operator<(const Region& other) const;
		bool operator==(const Region& other) const;

		uint lBound;					///< Lower bound of the region (offset from beginning of chromosome)
		uint rBound;					///< Upper bound of the region (offset from beginning of chromosome)

		//This isn't used to evaluate equality
		uint index;						///< Used to refer to the index within the region array 
	};

	struct RegionSet {
		RegionSet()  {}
		RegionSet(uint lBound, uint rBound, int index);
		RegionSet(const RegionSet& other);

		RegionSet& operator+=(const RegionSet& other);
		bool operator==(const RegionSet& other) const;

		/**
		 I think it makes better sense to simply store the index here, 
		 since the region itself has the other information...no need
		 to duplicate those endpoints
		 */
		std::set<Region> regions;
	};
	RegionContainer() {}
	RegionContainer(uint lBound, uint rBound, int index) {
		AddSegment(lBound, rBound, index);
	}
	RegionContainer(const RegionContainer& orig) : intervals(orig.intervals) {}
	virtual ~RegionContainer() {}
	
	void AddSegment(uint begin, uint end, uint id);

	typedef boost::icl::discrete_interval<uint> Interval;
	typedef boost::icl::split_interval_map<uint, RegionSet> IntervalMap;

	bool IsCoveredBySegment(uint point);
	
	/**
	 * Returns the ChromosomeRegions that overlap with the point
    * @param point Offset from beginning of chromosome
    * @param regions
    * @return 
    */
	bool GetRegionCoverage(uint point, std::set<Region>& regions);

protected:
	IntervalMap	intervals;
};


inline
bool RegionContainer::Region::operator<(const Region& other) const {
	if (lBound == other.lBound) 
		return rBound < other.rBound;
	return lBound < other.lBound;
}

inline
bool RegionContainer::Region::operator==(const Region& other) const {
	if (lBound == other.lBound)
		return rBound == other.rBound;
	return lBound == other.lBound;
}

inline
RegionContainer::RegionSet::RegionSet(uint lBound, uint rBound, int index) {
	regions.insert(RegionContainer::Region(lBound, rBound, index));
}

inline
RegionContainer::RegionSet::RegionSet(const RegionSet& other) : regions(other.regions) { }

inline
bool RegionContainer::RegionSet::operator==(const RegionSet& other) const {
	return regions == other.regions;
}

inline
RegionContainer::RegionSet& RegionContainer::RegionSet::operator+=(const RegionSet& other)  {
	regions.insert(other.regions.begin(), other.regions.end());
	return *this;
}

inline
void RegionContainer::AddSegment(uint begin, uint end, uint id) {
	intervals += std::make_pair(RegionContainer::Interval(begin, end, boost::icl::interval_bounds::closed()), RegionContainer::RegionSet(begin, end, id));
}

inline
bool RegionContainer::IsCoveredBySegment(uint point) {
	return boost::icl::contains(intervals, point);
}

/**
 * Returns the ChromosomeRegions that overlap with the point
 * @param point Offset from beginning of chromosome
 * @param regions
 * @return 
 */
inline
bool RegionContainer::GetRegionCoverage(uint point, std::set<RegionContainer::Region>& regions) {
	RegionContainer::IntervalMap::iterator itr = intervals.lower_bound(Interval(point, point, boost::icl::interval_bounds::closed()));
	RegionContainer::IntervalMap::iterator end = intervals.upper_bound(Interval(point, point, boost::icl::interval_bounds::closed()));

	while (itr != end) {
		regions.insert(itr->second.regions.begin(), itr->second.regions.end());
		itr++;
	}
	return regions.size() > 0;
}

}

#endif	/* REGIONCONTAINER_H */

