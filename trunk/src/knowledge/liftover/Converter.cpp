/* 
 * File:   converter.cpp
 * Author: torstees
 * 
 * Created on May 6, 2011, 1:25 PM
 */

#include "Converter.h"

#include <iostream>
#include <fstream>
#include <utility>						// std::make_pair
#include <boost/algorithm/string.hpp>

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

void Converter::ConvertDataset(const string& mapFilename, multimap<Locus*, Locus*>& converted) {

	vector<Locus*> loci;
	readMapFile(mapFilename, loci);
	ConvertDataset(loci.begin(), loci.end(), converted);
}

template <class T_cont>
void Converter::ConvertDataset(const string& mapFilename, T_cont& newBuild, std::ostream& droppedLog) {
	
	vector<Locus*> loci;
	readMapFile(mapFilename, loci);
	
	ConvertDataset(loci.begin(), loci.end(), newBuild, droppedLog);
	vector<Locus*>::iterator itr = loci.begin();
	vector<Locus*>::iterator end = loci.end();
	while(itr != end){
		delete *itr;
		++itr;
	}
}

template <class T_iter, class T_cont>
void Converter::ConvertDataset(T_iter itr, const T_iter& end, T_cont& newBuild, std::ostream& droppedLog) {
	multimap<Locus*, Locus*> converted;
	ConvertDataset(itr, end, converted);
	multimap<Locus*, Locus*>::const_iterator map_itr = converted.begin();
	multimap<Locus*, Locus*>::const_iterator map_end = converted.end();

	std::stringstream missingSNPs;
	while (map_itr != map_end) {
		//Report any bad SNPS
		if (map_itr->second && map_itr->second->getPos() == 0){
			missingSNPs<<*(itr->first);
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
		multimap<Locus*, Locus*>& converted) {
	//If there are no chains, then there is nothing to convert and we'll assume we
	//are at the right build already

	
	//SnpArray::const_iterator itr = originalBuild.begin();
	//SnpArray::const_iterator end = originalBuild.end();

	Locus* not_found = NULL;
	while (itr != end) {
		Locus& s = **itr;
		if (_chains.size() == 0){
			converted.insert(std::make_pair(&s,
							new Locus(s.getChrom(), s.getPos(), s.isRare(), s.getID())));
		} else {

			multimap<short,Chain*>::const_iterator citr =
					_chains.lower_bound(s.getChrom());
			multimap<short,Chain*>::const_iterator cend =
					_chains.upper_bound(s.getChrom());

			set<BuildConversion> con;

			//Accumulate the various conversions
			while (citr != cend){
				citr->second->EstimateConversion(s.getPos(), con);
				++citr;
			}

			//The best should "float" to the top
			if (con.size() > 0) {
				set<BuildConversion>::iterator bc = con.begin();
				Locus* newSNP = new Locus(bc->getRemoteChrom(), bc->getRemoteStart(), s.isRare(), s.getID());
				
				//std::cerr<<"ASDFASDF:"<<bc->rChrom.c_str()<<"\t"<<(int)newSNP.chr<<" "<<newSNP.pos<<" "<<newSNP.RSID()<<"\n";
				converted.insert(std::make_pair(&s, newSNP));
			} else {
				converted.insert(std::make_pair(&s, not_found));
			}
		}
		++itr;
	}
}

void Converter::readMapFile(const string& mapFilename, vector<Locus*>& array_out) const{

	std::ifstream file(mapFilename.c_str());
	if (!file.good()){
		throw std::ios_base::failure("File " + mapFilename + " not found");
	}

	string line;
	vector<string> words;

	while (file.good()) {
		words.clear();
		getline(file, line);
		split(words, line, is_any_of(" \n\t"), boost::token_compress_on);

		if (words.size() > 3){
			array_out.push_back(new Locus(words[0], atoi(words[3].c_str()), false, words[1]));
		}
	}
}

} // namespace Lifover
} // namespace Knowledge


