/* 
 * File:   region.cpp
 * Author: torstees
 * 
 * Created on April 28, 2011, 3:17 PM
 */

#include "RegionConverter.h"

#include "BuildConversion.h"

namespace Knowledge{
namespace Liftover {

RegionConverter::RegionConverter(int lStart, int rStart, int length, int num)
			: _localStart(lStart), _remoteStart(rStart),
			_length(length), _lineNumber(num) { }

int RegionConverter::estimate(int pos, bool direction) const {
	//Truncate any points occuring before this segment
	if (pos < _localStart) {
		return _remoteStart;
	}
	if (pos-_localStart > _length)
		return _remoteStart + (direction?1:-1)*_length;
	return _remoteStart + (direction?1:-1)*(pos-_localStart);
}

int RegionConverter::getLocalStart(int start) const {
	return (start > _localStart) ? start : _localStart;
}

int RegionConverter::getLocalStop(int stop) const {
	return (stop < _localStart + _length) ? stop :_localStart + _length;
}

bool RegionConverter::operator<(const RegionConverter& other) const  {
	return _localStart < other._localStart;
}



}
}

#ifdef TEST_APP

#include <gtest/gtest.h>
using namespace LiftOver;

TEST(LoRegionTest, RegionTest) {
	Region r1;
	Region r2(100, 300, 400, 1);
	Region r3(r2);
	Region r4(550, 700, 150, 2);
	r1 = r2;
	EXPECT_EQ(true, r1==r2);
	EXPECT_EQ(100, r1.localStart);
	EXPECT_EQ(300, r1.remoteStart);
	EXPECT_EQ(true, r2==r3);
	EXPECT_EQ(false, r2<r3);
	EXPECT_EQ(true, r3<r4);
	EXPECT_EQ(false, r4<r3);
	EXPECT_EQ(false, r3==r4);
	EXPECT_EQ(false, r2==r4);

	//We need to test that it does the right calculation
	//in both directions as well as truncates properly
	EXPECT_EQ(350, r2.Estimate(150, true));
	EXPECT_EQ(700, r2.Estimate(550, true));
	EXPECT_EQ(650, r4.Estimate(600, false));
	EXPECT_EQ(550, r4.Estimate(750, false));
}

TEST(LoRegionTest, RegionSetTestPositive) {
	RegionSet r1(200, 300, 200, 1);

	BuildConversion b("chr1", 150, 900);
	EXPECT_EQ(true, r1.EstimateConversion(true, 150, 900, b));
	EXPECT_EQ(300, b.rStart);
	EXPECT_EQ(500, b.rStop);

	RegionSet r2(450, 600, 150, 2);
	r1+=r2;
	EXPECT_EQ(2, r1.regions.size());
	EXPECT_EQ(true, r1.EstimateConversion(true, 150, 900, b));
	EXPECT_EQ(300, b.rStart);
	EXPECT_EQ(750, b.rStop);

	RegionSet r3(800, 1000, 200, 3);
	r1+=r3;
	EXPECT_EQ(true, r1.EstimateConversion(true, 150, 900, b));
	EXPECT_EQ(300, b.rStart);
	EXPECT_EQ(1100, b.rStop);
}

TEST(LoRegionTest, RegionSetTestNegative) {
	RegionSet r1(200, 1200, 200, 1);
	BuildConversion b("chr1", 150, 900);
	EXPECT_EQ(true, r1.EstimateConversion(false, 150, 900, b));
	EXPECT_EQ(1000, b.rStart);
	EXPECT_EQ(1200, b.rStop);

	RegionSet r2(450, 750, 150, 2);
	r1+=r2;
	EXPECT_EQ(true, r1.EstimateConversion(false, 150, 900, b));
	EXPECT_EQ(600, b.rStart);
	EXPECT_EQ(1200, b.rStop);

	RegionSet r3(800, 500, 200, 3);
	r1+=r3;
	EXPECT_EQ(true, r1.EstimateConversion(false, 150, 900, b));
	EXPECT_EQ(400, b.rStart);
	EXPECT_EQ(1200, b.rStop);

	//Verify that we truncate the end points correctly
	EXPECT_EQ(true, r1.EstimateConversion(false, 50, 1200, b));
	EXPECT_EQ(300, b.rStart);
	EXPECT_EQ(1200, b.rStop);
}


#endif

