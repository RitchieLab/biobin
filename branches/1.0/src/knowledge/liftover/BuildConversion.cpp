/*
 * BuildConversion.cpp
 *
 *  Created on: Nov 30, 2011
 *      Author: jrw32
 */

#include "BuildConversion.h"

namespace Knowledge{
namespace Liftover{
BuildConversion::BuildConversion(const string& chrom, int start, int stop)
		: _score(0), _lChrom(chrom), _lStart(start), _lStop(stop),
		_rChrom(""), _rStart(0), _rStop(0) {}

void BuildConversion::align() {
	if (_rStart > _rStop) {
		int t = _rStart;
		_rStart = _rStop;
		_rStop = t;
	}
}

// Loser is Lower--so, in a set, the one you want is the rbegin
bool BuildConversion::operator<(const BuildConversion& other) const {
	// order by score, length, other start
	if (_score == other._score)  {
		int length = _rStop-_rStart;
		int oLength = other._rStop-other._rStart;

		if (length < oLength)
			return _rStart < other._rStart;
		return length < oLength;
	}
	return _score < other._score;
}

}
}

