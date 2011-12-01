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

#ifndef KNOWLEDGE_LIFTOVER_REGIONCONVERTER_H
#define	KNOWLEDGE_LIFTOVER_REGIONCONVERTER_H

#include <set>
#include <string>

using std::set;
using std::string;

namespace Knowledge {
namespace Liftover {

class BuildConversion;

class RegionConverter {
public:
	RegionConverter(int lStart, int rStart, int length, int num);

	int estimate(int pos, bool direction) const;
	int getLocalStart(int start) const;
	int getLocalStop(int stop) const;
	
	bool operator<(const RegionConverter& other) const;
	bool operator==(const RegionConverter& other) const {
		return (_localStart == other._localStart) ?
				_lineNumber==other._lineNumber:
				_localStart==other._localStart;
	}

	template <class T_iter>
	static bool estimateConversion(T_iter itr, T_iter end,
			bool direction, int start, int stop, BuildConversion& c_out);

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
private:
	int _localStart;			///< Precalculated offset (required for conversion)
	int _remoteStart;			///< Precalculated offset
	int _length;					///< Subregion size

	int _lineNumber;			///< For debugging
};

}
}
#endif	/* REGION_H */

