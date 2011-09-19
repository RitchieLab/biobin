/* 
 * File:   regioncontainer.cpp
 * Author: torstees
 * 
 * Created on June 15, 2011, 10:33 AM
 */

#include "regioncontainer.h"




#ifdef TEST_APP

#include <gtest/gtest.h>
using namespace Knowledge;

TEST(RegionContainerTest, BasicFunctionality) {
	RegionContainer regions;
	regions.AddSegment(5,20,1);
	regions.AddSegment(25,30, 2);
	regions.AddSegment(15,20, 3);
	regions.AddSegment(35,40,4);
	
	EXPECT_TRUE(regions.IsCoveredBySegment(5));
	EXPECT_TRUE(regions.IsCoveredBySegment(10));
	EXPECT_TRUE(regions.IsCoveredBySegment(20));
	EXPECT_TRUE(regions.IsCoveredBySegment(36));
	EXPECT_FALSE(regions.IsCoveredBySegment(2));
	EXPECT_FALSE(regions.IsCoveredBySegment(45));
	
	std::set<RegionContainer::Region> covered;
	EXPECT_TRUE(regions.GetRegionCoverage(10,covered));
	EXPECT_EQ(1, covered.size());				// only entry in the set
	EXPECT_EQ(1, covered.begin()->index);
	covered.clear();
	
	EXPECT_TRUE(regions.GetRegionCoverage(16, covered));
	EXPECT_EQ(2, covered.size());
	EXPECT_EQ(3, (++covered.begin())->index);		// 2nd entry in the set
	
	covered.clear();
	EXPECT_TRUE(regions.GetRegionCoverage(25, covered));
	EXPECT_EQ(1, covered.size());				// only entry in the set
	EXPECT_EQ(2, covered.begin()->index);

	covered.clear();
	EXPECT_FALSE(regions.GetRegionCoverage(24, covered));
	
	
}

#endif
