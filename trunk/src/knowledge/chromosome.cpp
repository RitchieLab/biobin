/* 
 * File:   chromosome.cpp
 * Author: torstees
 * 
 * Created on March 3, 2011, 10:13 AM
 */

#include "chromosome.h"


#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

TEST(ChromosomeTest, RangedCollection) {
	Chromosome chrom;
	chrom.AddSNP(1,0);
	chrom.AddSNP(5,1);
	chrom.AddSNP(10,2);
	chrom.AddSNP(14,3);
	chrom.AddSNP(15,4);

	EXPECT_EQ(5, chrom.Size());
	Utility::IdCollection snps;
	chrom.GetSNPs(0,5, snps);
	EXPECT_EQ(2, snps.size());
	Utility::IdCollection::iterator itr=  snps.begin();
	EXPECT_EQ(0, *itr++);
	EXPECT_EQ(1, *itr);

	snps.clear();
	chrom.GetSNPs(5, 10, snps);
	EXPECT_EQ(2, snps.size());
	itr = snps.begin();
	EXPECT_EQ(1, *itr++);
	EXPECT_EQ(2, *itr);

	snps.clear();
	chrom.GetSNPs(10, 15, snps);
	EXPECT_EQ(3, snps.size());
	itr = snps.begin();
	EXPECT_EQ(2, *itr++);
	EXPECT_EQ(3, *itr++);
	EXPECT_EQ(4, *itr);
}



#include <boost/icl/interval_set.hpp>
using namespace boost;
using namespace boost::icl;

// I'm just picking this to host a quick test to learn to use the ICL
// This might go into the chromosome eventually, but for now, it's 
// just a playground for future reference
TEST(ChromosomeTest, IclTest) {
	typedef discrete_interval<int> Interval;
	interval_set<int> intervals;
	intervals += Interval(3, 5, interval_bounds::closed());
	intervals += Interval(4, 8, interval_bounds::closed());
	intervals += Interval(10, 12, interval_bounds::closed());

	EXPECT_TRUE(boost::icl::contains(intervals, (int)4));
	EXPECT_TRUE(boost::icl::contains(intervals, 5));
	EXPECT_TRUE(boost::icl::contains(intervals, 7));
	EXPECT_TRUE(boost::icl::contains(intervals, Interval(4,4)));
	EXPECT_TRUE(boost::icl::contains(intervals, 3));
	EXPECT_TRUE(boost::icl::contains(intervals, 8));
	EXPECT_FALSE(boost::icl::contains(intervals, 9));
	EXPECT_TRUE(boost::icl::contains(intervals, 12));
	EXPECT_FALSE(boost::icl::contains(intervals, 13));


/*
	discrete_interval<int> interval1 = construct<discrete_interval<int> >(3, 7, interval_bounds::closed());
	EXPECT_TRUE(contains(interval1, 4));
	EXPECT_FALSE(contains(interval1, 3));
*/
}

#endif //TEST_APP
