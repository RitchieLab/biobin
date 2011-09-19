/* 
 * File:   converter.cpp
 * Author: torstees
 * 
 * Created on May 6, 2011, 1:25 PM
 */

#include "converter.h"
#include "utility/locus.h"
#include <utility>						// std::make_pair
#include "utility/exception.h"

namespace LiftOver {

Converter::Converter(const char *orig, const char *newb) : origBuild(orig), newBuild(newb) {}
Converter::Converter(const Converter& orig) : chains(orig.chains), origBuild(orig.origBuild), newBuild(orig.newBuild) {}
Converter::~Converter() {}
void Converter::AddChain(const char *chainDetails) {
	Chain c;
	std::string chrom = c.Parse(chainDetails);

	chrom = chrom.erase(0,3);
	uint chr				= Utility::ChromToInt(chrom.c_str());
	//std::cerr<<" Chromosome: "<<chrom<<"\t"<<chr<<"\n";
	if (chr > 0)
		chains.insert(std::make_pair<uint, Chain>(chr, c));
}

void Converter::ConvertDataset(const char *mapFilename, std::multimap<SNP, SNP>& converted) {
	std::ifstream file(mapFilename);
	if (!file.good())
		throw Utility::Exception::FileNotFound(mapFilename);

	char line[4096];
	SnpArray snps;
	
	while (file.good()) {
		file.getline(line, 4096);
		Utility::StringArray words = Utility::Split(line);
		
		if (words.size() > 3)
			snps.push_back(SNP(Utility::ChromToInt(words[0].c_str()), atoi(words[3].c_str()), words[1]));
	}

	ConvertDataset(snps, converted);
}
void Converter::ConvertDataset(const char *mapFilename, SnpArray& newBuild, std::ostream& droppedLog) {
	std::ifstream file(mapFilename);
	if (!file.good())
		throw Utility::Exception::FileNotFound(mapFilename);

	char line[4096];
	SnpArray snps;
	
	while (file.good()) {
		file.getline(line, 4096);
		Utility::StringArray words = Utility::Split(line);
		
		if (words.size() > 3)
			snps.push_back(SNP(Utility::ChromToInt(words[0].c_str()), atoi(words[3].c_str()), words[1]));
	}
	
	ConvertDataset(snps, newBuild, droppedLog);
}

void Converter::ConvertDataset(const SnpArray& originalBuild, SnpArray& newBuild, std::ostream& droppedLog) {
	std::multimap<SNP, SNP> converted;
	ConvertDataset(originalBuild, converted);
	std::multimap<LiftOver::SNP, LiftOver::SNP>::iterator itr = converted.begin();
	std::multimap<LiftOver::SNP, LiftOver::SNP>::iterator end = converted.end();

	std::stringstream missingSNPs;
	while (itr != end) {
		//Report any bad SNPS
		if (itr->second.pos == 0) 
			missingSNPs<<itr->first.RSID()<<"\t"<<Utility::ChromFromInt(itr->first.chrom)<<"\t"<<itr->first.pos<<"\n";
		else 
			newBuild.push_back(itr->second);
		itr++;
	}

	if (missingSNPs.str().length() > 0) {
		droppedLog<<"SNPs that weren't translated properly to the new build:\n";
		droppedLog<<"Chr/tPos/tRS ID\n"<<missingSNPs.str()<<"\n";
	}
}

void Converter::ConvertDataset(const SnpArray& originalBuild, std::multimap<Utility::Locus, Utility::Locus>& converted) {
	//If there are no chains, then there is nothing to convert and we'll assume we
	//are at the right build already

	
	SnpArray::const_iterator itr = originalBuild.begin();
	SnpArray::const_iterator end = originalBuild.end();

	Utility::Locus unconvertable;
	while (itr != end) {
		const Utility::Locus &s = *itr;
		if (chains.size() == 0)
			converted.insert(std::make_pair<SNP, SNP>(s, s));
		else {

			ChainMap::iterator citr = chains.lower_bound(itr->chrom);
			ChainMap::iterator cend = chains.upper_bound(itr->chrom);
			Chain::ConversionSet con;


			//Accumulate the various conversions
			while (citr != cend)
				citr++->second.EstimateConversion(s.pos, con);

			//The best should "float" to the top
			if (con.size() > 0) {
				Chain::ConversionSet::iterator bc = con.begin();
				SNP newSNP = *itr;
				newSNP.chrom = Utility::ChromToInt(bc->rChrom.c_str());
				newSNP.pos   = bc->rStart;
				newSNP.RSID(s.RSID().c_str());
				
				//std::cerr<<"ASDFASDF:"<<bc->rChrom.c_str()<<"\t"<<(int)newSNP.chr<<" "<<newSNP.pos<<" "<<newSNP.RSID()<<"\n";
				converted.insert(std::make_pair<SNP, SNP>(s, newSNP));
			} else {
				converted.insert(std::make_pair<SNP, SNP>(s, unconvertable));
			}
		}
		itr++;
	}
}

}




#ifdef TEST_APP

#include <gtest/gtest.h>
using namespace LiftOver;

TEST(LoConversionTest, ConversionBasic) {
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

	Converter cnv("old", "new");
	cnv.AddChain(chunk.c_str());
	chunk = std::string("chain 870863 chr10 113275 + 9144 18556 chr10 135534747 - 46779061 46789268 4835\n")
		   + "301     2       0\n"
			+ "564     3       0\n"
		   + "802     1       0\n"
			+ "2624    0       112\n"
		   + "142     0       2\n"
		   + "128     100     787\n"
			+ "4745\n";

	cnv.AddChain(chunk.c_str());
	SnpArray snps;

	int chr = Utility::ChromToInt("10");

	SNP s1(chr, 81241464, "rs1", 1);			// 81251575
	SNP s2(chr, 81241564, "rs2", 1);			// 81251675
	SNP s3(chr, 81242948, "rs3", 1);			// 81253058
	SNP s4(chr, 81246876, "rs4", 1);			// 81256995
	SNP s5(chr, 9244,     "rs5", 1);			// 88755586
	SNP s6(chr, 13941,	 "rs6", 1);			// 88751783

	snps.push_back(s1);
	snps.push_back(s2);
	snps.push_back(s3);
	snps.push_back(s4);
	snps.push_back(s5);
	snps.push_back(s6);
	
	std::multimap<SNP, SNP> conversions;
	cnv.ConvertDataset(snps, conversions);

	SnpArray::iterator itr = snps.begin();
	SnpArray::iterator end = snps.end();

	while (itr != end) {
		//std::cout<<"chr"<<(int)itr->chr<<"\t"<<itr->pos<<"\t"<<itr->RSID()<<"\n";
		SNP s = *itr;

		std::multimap<SNP, SNP>::iterator citr = conversions.lower_bound(s);
		std::multimap<SNP, SNP>::iterator cend = conversions.upper_bound(s);

		while (citr != cend) {
			//std::cout<<"\tchr"<<(int)citr->second.chr<<"\t"<<citr->second.pos<<"\t"<<citr->second.RSID()<<"\n";
			citr++;
		}
		itr++;
	}

}

#endif
