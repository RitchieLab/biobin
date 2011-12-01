/* 
 * File:   converter.cpp
 * Author: torstees
 * 
 * Created on May 6, 2011, 1:25 PM
 */

#include "Converter.h"

#include <utility>						// std::make_pair
#include <boost/algorithm/string/split.hpp>

// TODO: remove dependency on utility/exception here
#include "utility/exception.h"

#include "Chain.h"
#include "knowledge/Locus.h"
#include "BuildConversion.h"

using boost::algorithm::split;
using boost::algorithm::is_any_of;

namespace Knowledge {
namespace Liftover {

Converter::Converter(const string& orig, const string& newb) :
		_origBuild(orig), _newBuild(newb) {}

Converter::~Converter() {
	multimap<short, Chain*>::iterator itr = _chains.begin();
	multimap<short, Chain*>::const_iterator end = _chains.end();
	while(itr != end){
		delete itr->second;
		++itr;
	}
}
void Converter::addChain(const string& chainDetails) {
	Chain *c = new Chain();
	std::string chrom = c->Parse(chainDetails);

	short chr_idx = Knowledge::Locus::getChrom(chrom);
	//std::cerr<<" Chromosome: "<<chrom<<"\t"<<chr<<"\n";
	if (chr_idx != -1){
		_chains.insert(std::make_pair<uint, Chain*>(chr_idx, c));
	}else{
		std::cerr<<"Unknown Chromosome in chain file!";
	}
}

void Converter::ConvertDataset(const string& mapFilename, multimap<Locus, Locus>& converted) {

	vector<Locus> loci;
	readMapFile(mapFilename, loci);

	ConvertDataset(loci.begin(), loci.end(), converted);
}

template <class T_cont>
void Converter::ConvertDataset(const string& mapFilename, T_cont& newBuild, std::ostream& droppedLog) {
	
	vector<Locus> loci;
	readMapFile(mapFilename, loci);
	
	ConvertDataset(loci.begin(), loci.end(), newBuild, droppedLog);
}

template <class T_iter, class T_cont>
void Converter::ConvertDataset(T_iter itr, const T_iter& end, T_cont& newBuild, std::ostream& droppedLog) {
	multimap<Locus, Locus> converted;
	ConvertDataset(itr, end, converted);
	multimap<Locus, Locus>::const_iterator map_itr = converted.begin();
	multimap<Locus, Locus>::const_iterator map_end = converted.end();

	std::stringstream missingSNPs;
	while (map_itr != map_end) {
		//Report any bad SNPS
		if (map_itr->second.getPos() == 0){
			missingSNPs<<itr->first;
		}else{
			newBuild.insert(itr->second);
		}
		++map_itr;
	}

	if (missingSNPs.str().length() > 0) {
		droppedLog<<"SNPs that weren't translated properly to the new build:\n";
		droppedLog<<"Chr/tPos/tRS ID\n"<<missingSNPs.str()<<"\n";
	}
}

template <class T_iter>
void Converter::ConvertDataset(T_iter itr, const T_iter& end,
		multimap<Locus, Locus>& converted) {
	//If there are no chains, then there is nothing to convert and we'll assume we
	//are at the right build already

	
	//SnpArray::const_iterator itr = originalBuild.begin();
	//SnpArray::const_iterator end = originalBuild.end();

	Locus unconvertable(-1,0,"Unconvertable");
	while (itr != end) {
		const Locus &s = *itr;
		if (_chains.size() == 0)
			converted.insert(std::make_pair(s, s));
		else {

			multimap<short,Chain*>::const_iterator citr =
					_chains.lower_bound(itr->getChrom());
			multimap<short,Chain*>::const_iterator cend =
					_chains.upper_bound(itr->getChrom());

			set<BuildConversion> con;

			//Accumulate the various conversions
			while (citr != cend){
				citr->second->EstimateConversion(s.getPos(), con);
				++citr;
			}

			//The best should "float" to the top
			if (con.size() > 0) {
				set<BuildConversion>::iterator bc = con.begin();
				Locus newSNP(bc->getRemoteChrom(), bc->getRemoteStart(), s.getID());
				
				//std::cerr<<"ASDFASDF:"<<bc->rChrom.c_str()<<"\t"<<(int)newSNP.chr<<" "<<newSNP.pos<<" "<<newSNP.RSID()<<"\n";
				converted.insert(std::make_pair(s, newSNP));
			} else {
				converted.insert(std::make_pair(s, unconvertable));
			}
		}
		++itr;
	}
}

void Converter::readMapFile(const string& mapFilename, vector<Locus>& array_out) const{

	std::ifstream file(mapFilename.c_str());
	if (!file.good())
		throw Utility::Exception::FileNotFound(mapFilename.c_str());

	char line[4096];
	vector<string> words;

	while (file.good()) {
		words.clear();
		file.getline(line, 4096);
		split(words, line, is_any_of(" \n\t"));

		if (words.size() > 3){
			array_out.push_back(Locus(words[0], atoi(words[3].c_str()), words[1]));
		}
	}
}

} // namespace Lifover
} // namespace Knowledge




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
