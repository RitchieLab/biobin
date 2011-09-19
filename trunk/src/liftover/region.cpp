/* 
 * File:   region.cpp
 * Author: torstees
 * 
 * Created on April 28, 2011, 3:17 PM
 */

#include "region.h"
namespace LiftOver {

RegionSet::RegionSet() {}

RegionSet::RegionSet(int lStart, int rStart, int length, int num) {
	regions.insert(Region(lStart, rStart, length, num));
}

RegionSet::RegionSet(const RegionSet& other) {
	regions = other.regions;
}

bool RegionSet::operator==(const RegionSet& other) const {
	return regions == other.regions;
}

RegionSet& RegionSet::operator+=(const RegionSet& other) {
	regions.insert(other.regions.begin(), other.regions.end());
	return *this;
}

bool RegionSet::EstimateConversion(bool direction, int start, int stop, BuildConversion& c) {
//void RegionSet::EstimateConversion(const char *lchrom, const char *rchrom,
//		bool direction, int start, int stop,
//		std::set<BuildConversion>& conversionOptions) {

	if (regions.size() > 0) {
		//BuildConversion c(lchrom, start, stop);
		//c.rChrom = rchrom;
		std::set<Region>::iterator r = regions.begin();
		c.lStart = r->GetLocalStart(start, direction);
		c.rStart = r->Estimate(start, direction);

		if (regions.size() == 1) {
			c.rStop = r->Estimate(stop, direction);
			c.lStop	= r->GetLocalStop(stop, direction);
		}
		else {
			std::set<Region>::reverse_iterator s = regions.rbegin();
			c.rStop = s->Estimate(stop, direction);
			c.lStop	= s->GetLocalStop(stop, direction);
		}
		//conversionOptions.insert(c);
		if (c.rStop < c.rStart) {
			int t = c.rStop;
			c.rStop = c.rStart;
			c.rStart = t;
		}
		return true;
	}
	return false;
}

BuildConversion::BuildConversion(const char *chrom, int start, int stop)
		: score(0), lChrom(chrom), lStart(start), lStop(stop),
		rChrom(""), rStart(0), rStop(0) {}

BuildConversion::BuildConversion() : score(0),
		lChrom(""), lStart(0), lStop(0),
		rChrom(""), rStart(0), rStop(0) {}

void BuildConversion::Align() {
	if (rStart > rStop) {
		int t = rStart;
		rStart = rStop;
		rStop = t;
	}
}

// Loser is Lower--so, in a set, the one you want is the rbegin
bool BuildConversion::operator<(const BuildConversion& other) const {
	// order by score, length, other start
	if (score == other.score)  {
		int length = rStop-rStart;
		int oLength = other.rStop-other.rStart;

		if (length < oLength)
			return rStart < other.rStart;
		return length < oLength;
	}
	return score < other.score;
}

Region::Region() : localStart(0), remoteStart(0), length(0), lineNumber(0) {}

Region::Region(int lStart, int rStart, int length, int num)
			: localStart(lStart), remoteStart(rStart),
			length(length), lineNumber(num) { }

Region::Region(const Region& other)  : localStart(other.localStart),
			remoteStart(other.remoteStart), length(other.length),
			lineNumber(other.lineNumber) {}

int Region::Estimate(int pos, bool direction) const {
	//Truncate any points occuring before this segment
	if (pos < localStart) {
		return remoteStart;
	}
	if (pos-localStart > length)
		return remoteStart + (direction?1:-1)*length;
	return remoteStart + (direction?1:-1)*(pos-localStart);
}

int Region::GetLocalStart(int start, bool direction) const {
	if (start > localStart)
		return start;
	return localStart;
}

int Region::GetLocalStop(int stop, bool direction) const {
	if (stop < localStart + length)
		return stop;
	return localStart + length;
}

bool Region::operator<(const Region& other) const  {
	return localStart < other.localStart;
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

