/* 
 * File:   converter.h
 * Author: torstees
 *
 * Created on May 6, 2011, 1:25 PM
 */

#ifndef KNOWLEDGE_LIFTOVER_CONVERTER_H
#define	KNOWLEDGE_LIFTOVER_CONVERTER_H

#include <utility>
#include <set>
#include <map>

#include "knowledge/Locus.h"

using Knowledge::Locus;
using std::pair;
using std::set;
using std::map;

namespace Knowledge {

namespace Liftover{

class Chain;

class Converter {
public:

	static const pair<short, pair<int, int> > FAILED_REGION;
	static const float MIN_MAPPING_FRACTION;

public:

	Converter(const string& origBuild);
	virtual ~Converter();

	/**
	 * Loads all of the chain files from a given source (pure virtual)
	 */
	virtual int Load() = 0;

	/**
	 * \brief Converts a region using the loaded chain files.
	 * This function will lift a given region on a chromosome to another region
	 * in a new build, using the loaded chains and segments.  If a conversion
	 * fails, this will return (-1, (-1,-1) ), or the FAILED_REGION.
	 *
	 * NOTE: It is imperative that the start and end positions are at least 1
	 * apart, as in the standard liftOver algorithm.  if you wish to lift a
	 * single position, call convertRegion(chr, pos, pos+1) and simply take
	 * the first element of the second element in the return value.
	 *
	 * \param chrom The old chromosome
	 * \param start The starting position of the region
	 * \param end The ending position of the region
	 *
	 * \return A pair containing the new chromosome and a pair of the
	 * new region's boundaries.  If unable to map, returns the FAILED_REGION
	 */
	pair<short, pair<int, int> > convertRegion(short chrom, int start, int end) const;

	/**
	 * \brief Converts an (iterable) set of Loci.
	 * This function takes two iterators that define a collection of Locus*
	 * objects and converts them using the chain files.  This function will
	 * create new Locus objects and insert them into the locus_map_out container
	 * which is assumed to be a map of Locus* to Locus* objects (or an object
	 * that conforms to the map specification).  The unmapped_out is a container
	 * that will hold the Locus* objects that were not mapped to any new Locus.
	 *
	 * \param itr The iterator to the beginning of the Locus* sequence
	 * \param end An iterator to the end of the Locus* sequence
	 * \param locus_map_out A mapping of old Locus* to newly converted Locus*
	 * \param unmapped_out A container of unmapped Locus* objects
	 */
	template <class T_iter, class T_map, class T_cont>
	void convertLoci(T_iter itr, const T_iter& end, T_map& locus_map_out, T_cont& unmapped_out) const;

protected:
	// A mapping of chromosome -> chains, ordered by score
	map<short, set<Chain*> > _chains;

	// The string describing the original build to lift from
	string _origBuild;

private:
	// No copying or assignment allowed
	Converter(const Converter& orig);
	Converter& operator=(const Converter& other);

};

template <class T_iter, class T_map, class T_cont>
void Converter::convertLoci(T_iter itr, const T_iter& end, T_map& locus_map_out, T_cont& unmapped_out) const{

	while (itr != end){
		pair<short, pair<int, int> > new_region =
				convertRegion((*itr)->getChrom(), (*itr)->getPos(), (*itr)->getPos()+1);

		if (new_region == FAILED_REGION){
			unmapped_out.insert(unmapped_out.end(), *itr);
		}else{
			Locus* converted =
					new Locus(new_region.first, new_region.second.first,
							(*itr)->isRare(), (*itr)->getID());
			converted->addAlleles((*itr)->beginAlleles(), (*itr)->endAlleles());
			locus_map_out[*itr] = converted;
		}

		++itr;
	}
}



}
}

#endif	/* CONVERTER_H */

