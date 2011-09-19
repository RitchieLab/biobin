/* 
 * File:   chromosome.h
 * Author: torstees
 *
 * This is basically a way to quickly identify the SNP indexes within a range
 * It is a wrapper around a map and some basic functions, and nothing more
 *
 * Created on March 3, 2011, 10:13 AM
 */

#ifndef CHROMOSOME_H
#define	CHROMOSOME_H

#include "utility/locus.h"
#include "def.h"
#include <map>
#include <assert.h>
//#include <boost/icl/interval_set.hpp>
#include "regioncontainer.h"

namespace Knowledge {


class Chromosome {
public:	
	Chromosome() : offset(0) { }
	Chromosome(const Chromosome& orig) : offset(orig.offset), snps(orig.snps), segments(orig.segments) { }
	virtual ~Chromosome() {}

	uint offset;										///< In case we need to downcast positions that are genome oriented

	void Clear();
	/**
	 Collections SNP Indexes bound by the left and right boundaries, returning number matched

	 * @note The set is not cleared before hand, so it just adds SNPs to the set
	 */
	void GetSNPs(uint leftbound, uint rightbound, Utility::IdCollection& subset) const;

	void AddSNP(uint pos, uint idx);

	uint Size();

	std::vector<uint> GetSnpIndexes();


	/************************************************************************
	 * Region overlap tree - Client can build regions on a chromosome and check
	 * to see if a point is present in at least one of the regions. I switched
	 * to use the regionContainer which allows me to return region indexes. 
	 * 
	 * NOTE: This is a late comer in terms of functionality, so there are probably
	 * a number of places where we are doing something fairly brute force, that 
	 * can be switched to use this instead. I am not going to attempt to refactor
	 * anything to use this due to time constraints...but don't think it isn't 
	 * being used because it shouldn't be used....It probably is just code that
	 * came first. 
    */
	void AddSegment(uint begin, uint end, uint id);
	bool IsCoveredBySegment(uint point);

	/**
	 * Adds region indexes to the collection which contain the position, pos
    * @param pos
    * @param regionIndexes
    */
	bool GetSegmentCoverage(uint pos, Utility::IdCollection& regionIndexes);
private:
	
	//typedef boost::icl::discrete_interval<int> Segment;
	//typedef boost::icl::interval_set<int> SegmentContainer;
	typedef std::multimap<uint, uint> SnpLookup;
	SnpLookup snps;						///< position => snp_index
	RegionContainer segments;
};


inline
uint Chromosome::Size() {
	return snps.size();
}

inline
bool Chromosome::GetSegmentCoverage(uint pos, Utility::IdCollection& regionIndexes) {
	std::set<RegionContainer::Region> regions;
	segments.GetRegionCoverage(pos, regions);
	std::set<RegionContainer::Region>::iterator itr = regions.begin();
	std::set<RegionContainer::Region>::iterator end = regions.end();
	
	while (itr != end) {
		regionIndexes.insert(itr->index);
		itr++;
	}
	return regions.size() > 0;
}

inline
std::vector<uint> Chromosome::GetSnpIndexes() {
	Chromosome::SnpLookup::iterator itr = snps.begin();
	SnpLookup::iterator end = snps.end();
	std::vector<uint> snpData;
	while (itr != end) 
		snpData.push_back(itr++->second);
	return snpData;
}

inline
void Chromosome::AddSegment(uint begin, uint end, uint id) {
	segments.AddSegment(begin, end, id);
	//segments += Segment(begin, end, boost::icl::interval_bounds::closed());
}

inline
bool Chromosome::IsCoveredBySegment(uint point) {
	return segments.IsCoveredBySegment(point);
	//return boost::icl::contains(segments, point);
}


inline
void Chromosome::GetSNPs(uint left, uint right, Utility::IdCollection& subset) const {
	SnpLookup::const_iterator itr = snps.lower_bound(left);
	SnpLookup::const_iterator end = snps.upper_bound(right);
	//uint snpCount = snps.size();
	while (itr != end) {
		subset.insert(itr++->second);
	}
}

inline
void Chromosome::AddSNP(uint pos, uint idx) {
	//If it dies here, then we have redundant SNP positions
//	if (snps.find(pos) != snps.end())
//		std::cerr<<"Duplicate SNP found: "<<pos<<" "<<idx<<"\n";
	snps.insert(std::pair<uint, uint>(pos, idx));
	//snps[pos] = idx;
}

inline
void Chromosome::Clear() {
	snps.clear();
	offset = 0;
}


}


#endif	/* CHROMOSOME_H */

