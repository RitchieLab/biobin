/* 
 * File:   chain.cpp
 * Author: torstees
 *
 * Represents a complete chain as described in the liftover chain files
 * A chromosome can have many chains.
 *
 * Created on April 28, 2011, 4:08 PM
 */

#include <boost/algorithm/string.hpp>
#include <vector>
#include <utility>
#include <stdlib.h>

#include "Chain.h"

#include "BuildConversion.h"

using boost::algorithm::split;
using boost::algorithm::is_any_of;
using std::vector;

namespace Knowledge {
namespace Liftover {

Chain::Chain()
	: _rChrom(""), _forward(true), _id(-1), _score(0),
		  _lLength(0), _lOffset(0), _rLength(0), _rOffset(0) {}

Chain::Chain(int id, long score, int lLength, int lOffset,
	const string& chrom, int rLength, int rOffset)
	: _rChrom(chrom), _forward(true), _id(-1), _score(0),
		  _lLength(0), _lOffset(0), _rLength(0), _rOffset(0) {}

int Chain::Count() {
	return _regions.size();
}

void Chain::EstimateConversion(int pos, set<BuildConversion>& conversionOptions) {
	BuildConversion c(_lChrom, pos, pos);
	if (EstimateConversion(pos, c))
		conversionOptions.insert(c);
}

bool Chain::EstimateConversion(int pos, BuildConversion& c_out) {

	return EstimateConversion(pos,pos,c_out);

	// TODO: Why is this necessary?
	int pos_expand = 5;
	BuildConversion test(c_out);

	EstimateConversion(pos - pos_expand, pos + pos_expand, test);
	if (test.getRemoteStop() - test.getRemoteStart() > -1){
		return EstimateConversion(pos,pos,c_out);
	}
	return false;


	/*
	RegionMap::const_iterator itr =	_regions.lower_bound(
			interval<uint>::closed(pos-pos_expand,pos+pos_expand));
	RegionMap::const_iterator end = _regions.upper_bound(
			interval<uint>::closed(pos-pos_expand,pos+pos_expand));

	//std::cerr<<"EstimateConversion("<<pos<<","<<c.lChrom<<":"<<pos-1<<"-"<<pos+1<<") on : "<<lChrom<<":"<<lOffset<<"-"<<lOffset+lLength<<"\n";
	set<RegionConverter> region_set;
	if (itr != region_set.end()) {
		//std::cerr<<"Accumulating region information:";
		while (itr != end) {
			region_set+=itr->second;
			++itr;
		}
		c.lChrom		= lChrom;
		c.rChrom		= rChrom;
		c.lStart		= pos;
		c.lStop		= pos;
		
		//Run a test-make sure the region we get back isn't too much smaller than the original
		BuildConversion test;
		test.lStart	= pos-5;
		test.lStop  = pos+5;
		region.EstimateConversion(direction==1, pos-5, pos+5, test);
		// TODO We need to be able to restrict these regions, but if this is too
		// restrictive, we'll miss edges... I had changed this value to 9 to
		// make large scale tests work more like what we see in liftover-however,
		// they caused errors to my tests here...
		if (test.rStop - test.rStart > 4)
			return region.EstimateConversion(direction == 1, pos, pos, c);
		return false;
	}
	//std::cerr<<"No region information available\n";
	return false;
	*/
}

void Chain::EstimateConversion(int start, int stop, std::set<BuildConversion>& conversionOptions) {
	BuildConversion c(_lChrom, start, stop);
	if (EstimateConversion(start, stop, c))
		conversionOptions.insert(c);
}

bool Chain::EstimateConversion(int start, int stop, BuildConversion& c) {
	RegionMap::const_iterator itr = _regions.lower_bound(interval<uint>::closed(start, stop));
	RegionMap::const_iterator end = _regions.upper_bound(interval<uint>::closed(start,stop));

	set<RegionConverter> region_set;
	if (itr != _regions.end()) {
		while (itr != end){
			region_set.insert(itr->second.begin(), itr->second.end());
			++itr;
		}

		c.setChroms(_lChrom, _rChrom);
		// I don't need this (I think); it should be set correctly already
		//c.setLocalStart(start);
		//c.setLocalStop(stop);
		return RegionConverter::estimateConversion(
				region_set.begin(), region_set.end(),_forward, start, stop, c);
	}
	return false;
}

std::string Chain::Parse(const string& data) {
	vector<string> lines;
	split(lines, data, is_any_of("\n"), boost::token_compress_on);

	vector<string>::iterator line_itr = lines.begin();
	vector<string>::const_iterator line_end = lines.end();

	vector<string> words;
	split(words, *line_itr, is_any_of(" \t\n"), boost::token_compress_on);

	_score = atoi(words[1].c_str());
	_lChrom = words[2];
	_forward = words[9]=="+";
	_id = atoi(words[12].c_str());

	_lOffset = atoi(words[5].c_str());
	_lLength = atoi(words[6].c_str()) - _lOffset;
	_rOffset = atoi(words[10].c_str());
	_rLength = atoi(words[11].c_str()) - _rOffset;

	if (!_forward){
		_rOffset = atoi(words[8].c_str()) - _rOffset;
	}

	_rChrom = words[7];

	int ldiff = _lOffset;
	int rdiff = _rOffset;
	int size = 0;
	int i=0;

	while (++line_itr != line_end) {
		words.clear();
		split(words,*line_itr,is_any_of(" \t\n"), boost::token_compress_on);
		size = atoi(words[0].c_str());

		//At this time, we are always moving left to right on the local side (that
		//is true of liftover as of the time of initial development)
		//So, the only thing that we really have to worry about being confusing
		//is the remote side.
		set<RegionConverter> new_set;
		new_set.insert(RegionConverter(ldiff,rdiff,size,++i));
		_regions += std::make_pair(
				interval<uint>::right_open(ldiff, ldiff + size),new_set);
		//std::cerr<<" "<<i<<" ( "<<lines[i]<<" ) "<<ldiff<<" -- "<<ldiff+size<<"\t"<<rdiff<<" - "<<rdiff+(size*direction)<<"\t"<<regions.size()<<"\n";
		if (words.size() == 3) {
			ldiff += size + atoi(words[1].c_str());
			rdiff += (size + atoi(words[2].c_str())) * Direction();
		}
	}

	return _lChrom;
}

}
}



#ifdef TEST_APP

#include <gtest/gtest.h>
using namespace LiftOver;

TEST(LoChainTest, ChainTestPositive) {
	std::string chunk;
	chunk = std::string("chain 788625 chr10 135374737 + 81241464 81249852 chr10 135534747 + 81251575 81259959 6147\n")
		   + "783     1       0\n"
			+ "883     23      33\n"
		   + "3247    1       0\n"
			+ "143     31      31\n"
		   + "957     15      3\n"
			+ "2304\n";

	/**
	 * This results in the following local / remote pairs:
	 *
	 * 81241464	81242247	-	81251575	81252358
	 * 81242248	81243131	-	81252358	81253241
	 * 81243154	81246401	-	81253274	81256521
	 * 81246402	81246545	-	81256521	81256664
	 * 81246576	81247533	-	81256695	81257652
	 * 81247548	81249852	-	81257655	81259959
	 */

	Chain c;
	c.Parse(chunk.c_str());
	//EXPECT_EQ(6, c.Count());

	EXPECT_EQ(6147, c.ID());
	EXPECT_EQ(8388, c.LLength());
	EXPECT_EQ(788625, c.Score());
	EXPECT_EQ(81241464, c.LStart());
	EXPECT_EQ(81251575, c.RStart());
	EXPECT_EQ(1, c.Direction());

	BuildConversion bc;

	//First just some boring inside by 30 on each end tests on each subchain
	EXPECT_TRUE(c.EstimateConversion(81241494, 81242217, bc));
	EXPECT_EQ(81241494, bc.lStart);
	EXPECT_EQ(81242217, bc.lStop);
	EXPECT_EQ(81251605, bc.rStart);
	EXPECT_EQ(81252328, bc.rStop);
	
	EXPECT_TRUE(c.EstimateConversion(81241494, 81241495, bc));
	EXPECT_EQ(81251605, bc.rStart);
	EXPECT_EQ(81251606, bc.rStop);
	
	EXPECT_TRUE(c.EstimateConversion(81241464, 81241465, bc));
	EXPECT_EQ(81251575, bc.rStart);
	EXPECT_EQ(81251576, bc.rStop);
	
	EXPECT_TRUE(c.EstimateConversion(81241464, bc));
	EXPECT_EQ(81251575, bc.rStart);
	EXPECT_EQ(81251575, bc.rStop);
	
	
	EXPECT_TRUE(c.EstimateConversion(81242278, 81243101, bc));
	EXPECT_EQ(81252388, bc.rStart);
	EXPECT_EQ(81253211, bc.rStop);

	EXPECT_TRUE(c.EstimateConversion(81241463, 81241465, bc));
	EXPECT_EQ(81251575, bc.rStart);
	EXPECT_EQ(81251576, bc.rStop);

	c.EstimateConversion(81243184, 81246371, bc);
	EXPECT_EQ(81253304, bc.rStart);
	EXPECT_EQ(81256491, bc.rStop);

	c.EstimateConversion(81246432, 81246515, bc);
	EXPECT_EQ(81256551, bc.rStart);
	EXPECT_EQ(81256634, bc.rStop);

	c.EstimateConversion(81246606, 81247503, bc);
	EXPECT_EQ(81256725, bc.rStart);
	EXPECT_EQ(81257622, bc.rStop);

	c.EstimateConversion(81247578, 81249822, bc);
	EXPECT_EQ(81257685, bc.rStart);
	EXPECT_EQ(81259929, bc.rStop);

	//Now let's mess around with boundaries
	//Starts before the first subchain. This should be truncated to the start
	//of the subchain
	EXPECT_TRUE(c.EstimateConversion(81000000, 81242217, bc));
	EXPECT_EQ(81241464, bc.lStart);
	EXPECT_EQ(81242217, bc.lStop);
	EXPECT_EQ(81251575, bc.rStart);
	EXPECT_EQ(81252328, bc.rStop);

	//And the other side...extends beyond the right side of the second subchain
	EXPECT_TRUE(c.EstimateConversion(81242278, 81243144, bc));
	EXPECT_EQ(81252388, bc.rStart);
	EXPECT_EQ(81253241, bc.rStop);

	//And a region the extends across first 3 subchains
	EXPECT_TRUE(c.EstimateConversion(81241494, 81246371, bc));
	EXPECT_EQ(81241494, bc.lStart);
	EXPECT_EQ(81256491, bc.rStop);

	//And the truncations where that extend beyond multiple subchains
	EXPECT_TRUE(c.EstimateConversion(81000000, 81246556, bc));
	c.EstimateConversion(81000000, 81246556, bc);
	EXPECT_EQ(81241464, bc.lStart);
	EXPECT_EQ(81246545, bc.lStop);
	EXPECT_EQ(81251575, bc.rStart);
	EXPECT_EQ(81256664, bc.rStop);
}

TEST(LoChainTest, ChainTestNegative) {
	std::string chunk;
	chunk = std::string("chain 870863 chr10_random 113275 + 9144 18556 chr10 135534747 - 46779061 46789268 4835\n")
		   + "301     2       0\n"
			+ "564     3       0\n"
		   + "802     1       0\n"
			+ "2624    0       112\n"
		   + "142     0       2\n"
		   + "128     100     787\n"
			+ "4745\n";

	/**
	 * This results in the following local / remote pairs:
				9144	9445		88755686	88755385
				9447	10011		88755385	88754821
				10014	10816		88754821	88754019
				10817	13441		88754019	88751395
				13441	13583		88751283	88751141
				13583	13711		88751139	88751011
				13811	18556		88750224	88745479
			*/


	Chain c;
	c.Parse(chunk.c_str());

	EXPECT_EQ(4835, c.ID());
	EXPECT_EQ(9412, c.LLength());
	EXPECT_EQ(870863, c.Score());
	EXPECT_EQ(9144, c.LStart());
	EXPECT_EQ(88755686, c.RStart());
	EXPECT_EQ(-1, c.Direction());

	BuildConversion bc;

	//First just some boring inside by 30 on each end tests on each subchain
	/*
				9174	9415		88755656	88755415
				9477	9981		88755355	88754851
				10044	10786		88754791	88754049
				10847	13411		88753989	88751425
				13471	13553		88751253	88751171
				13613	13681		88751109	88751041
				13841	18526		88750194	88745509
	 */
	EXPECT_TRUE(c.EstimateConversion(9174, 9415, bc));
	EXPECT_EQ(9174, bc.lStart);
	EXPECT_EQ(9415, bc.lStop);

	//The pairs are flipped when we are on opposite strand
	EXPECT_EQ(88755656, bc.rStop);
	EXPECT_EQ(88755415, bc.rStart);

	EXPECT_TRUE(c.EstimateConversion(9477, 9981, bc));
	EXPECT_EQ(88755355, bc.rStop);
	EXPECT_EQ(88754851, bc.rStart);

	EXPECT_TRUE(c.EstimateConversion(10044, 10786, bc));
	EXPECT_EQ(88754791, bc.rStop);
	EXPECT_EQ(88754049, bc.rStart);

	EXPECT_TRUE(c.EstimateConversion(10847, 13411, bc));
	EXPECT_EQ(88753989, bc.rStop);
	EXPECT_EQ(88751425, bc.rStart);

	EXPECT_TRUE(c.EstimateConversion(13471, 13553, bc));
	EXPECT_EQ(88751253, bc.rStop);
	EXPECT_EQ(88751171, bc.rStart);

	EXPECT_TRUE(c.EstimateConversion(13613, 13681, bc));
	EXPECT_EQ(88751109, bc.rStop);
	EXPECT_EQ(88751041, bc.rStart);

}


#endif
