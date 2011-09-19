/* 
 * File:   region.h
 * Author: torstees
 *
 * Region:
 * Represents a single subchain from UCSC's chain file and is ultimately
 * responsible for converting points from one build to the next using
 * the subchain's offset details.
 *
 * This class doesn't contain all information necessary for conversion, however.
 * It is assumed to be used in conjuction with the chain object, since it
 * doesn't store details about strand.
 *
 *
 * RegionSet:
 * This helper class is used to aggregate the various non-contiguous regions
 * mapped on the chain. It is just a way to keep up with all of the subchains
 * from a single query in one place.
 *
 * Created on April 28, 2011, 3:17 PM
 */

#ifndef LREGION_H
#define	LREGION_H
#include <set>
#include <string>

namespace LiftOver {

struct BuildConversion {
	BuildConversion(const char *chrom, int start, int stop);
	BuildConversion();

	void Align();

	bool operator<(const BuildConversion& other) const;


	long score;					///< The conversion's score
	std::string lChrom;		///< Local Chromosome
	int lStart;					///< Local Start
	int lStop;					///< Local Stop
	std::string rChrom;		///< Remote Chromosome
	int rStart;					///< Remote Start
	int rStop;					///< Remote Stop
};

typedef std::set<BuildConversion> ConversionSet;

struct Region {
	Region();
	Region(int lStart, int rStart, int length, int num);
	Region(const Region& other);

	int Estimate(int pos, bool direction) const;
	int GetLocalStart(int start, bool direction) const;
	int GetLocalStop(int stop, bool direction) const;
	
	bool operator<(const Region& other) const;
	bool operator==(const Region& other) const {
		if (localStart == other.localStart)
			return lineNumber==other.lineNumber;
		return localStart==other.localStart;
	}


	/**
	 * Offsets are pretty straightforward, except possibly the remote
	 * offset, which can tricky. These are predetermined by the
	 * object instantiating the regions, and in the case of inverted
	 * regions (direction is false) the remoteStart descends with each
	 * sub-region loaded.
	 *
	 * As stated elsewhere, direction is stored in the chain and propagates
	 * downward during estimation. 
	 */
	int localStart;			///< Precalculated offset (required for conversion)
	int remoteStart;			///< Precalculated offset
	int length;					///< Subregion size

	int lineNumber;			///< For debugging
};




struct RegionSet {
	RegionSet();
	RegionSet(int lStart, int rStart, int length, int num);
	RegionSet(const RegionSet& other);

	RegionSet& operator+=(const RegionSet& other);
	bool operator==(const RegionSet& other) const;
	bool EstimateConversion(bool direction, int start, int stop, BuildConversion& c);

	std::set<Region> regions;
};

}
#endif	/* REGION_H */

