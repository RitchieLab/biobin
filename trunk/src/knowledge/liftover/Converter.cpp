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


