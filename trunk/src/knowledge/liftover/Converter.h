/* 
 * File:   converter.h
 * Author: torstees
 *
 * Created on May 6, 2011, 1:25 PM
 */

#ifndef KNOWLEDGE_LIFTOVER_CONVERTER_H
#define	KNOWLEDGE_LIFTOVER_CONVERTER_H

#include <vector>
#include <string>
#include <ostream>
#include <map>

#include "Chain.h"
#include "knowledge/Locus.h"
#include "BuildConversion.h"

#include <iostream>
#include <fstream>
#include <utility>						// std::make_pair
#include <boost/algorithm/string.hpp>

using std::multimap;
using std::string;
using std::vector;
using boost::algorithm::split;
using boost::algorithm::is_any_of;

namespace Knowledge {


namespace Liftover{

using Knowledge::Locus;

class Converter {
public:
	Converter(const string& origBuild, const string& newBuild);
	virtual ~Converter();

	void addChain(const string& chainDetails);
	
	void ConvertDataset(const string& mapFilename, multimap<Locus*, Locus*>& converted);
	/**
	 * Takes a list (iterable) of Locus object ptrs, applies the chaining and
	 * returns the converted Locus objects in the supplied container.  Also
	 * logs any dropped Locus objects to the supplied ostream object.
	 */
	template <class T_iter, class T_cont>
	void ConvertDataset(T_iter itr, const T_iter& end, T_cont& newBuild, std::ostream& droppedLog);

	/**
	 * Takes a map filename and returns the new Locus objects in the supplied
	 * container.
	 */
	template <class T_cont>
	void ConvertDataset(const string& mapFilename, T_cont& newBuild, std::ostream& droppedLog);

	/**
	 * Takes a list (iterable) of Locus object ptrs and chains them together and
	 * returns the result as a multimap
	 */
	template <class T_iter>
	void ConvertDataset(T_iter itr, const T_iter& end,
			multimap<Locus*, Locus*>& converted);

	/**
	 * Loads all of the chain files from a given source (pure virtual)
	 */
	virtual int Load() = 0;


protected:
	multimap<short, Chain*> _chains;					///< Chrom -> chains
	string _origBuild;
	string _newBuild;

private:
	// No copying or assignment allowed
	Converter(const Converter& orig);
	Converter& operator=(const Converter& other);

	// Function to read a map filename and return the results in a vector of
	// Locus Objects
	void readMapFile(const string& mapFilename, vector<Locus*>& array_out) const;
};

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


}
}

#endif	/* CONVERTER_H */

