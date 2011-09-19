/* 
 * File:   chain.h
 * Author: torstees
 *
 * Represents a complete chain as described in the liftover chain files
 * A chromosome can have many chains.
 *
 * Created on April 28, 2011, 4:08 PM
 */

#ifndef CHAIN_H
#define	CHAIN_H

#include <boost/icl/interval_set.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <string>

#include "region.h"
#include "utility/strings.h"

namespace LiftOver {

class Chain {
public:

typedef boost::icl::split_interval_map<uint , RegionSet> RegionMap;
typedef boost::icl::discrete_interval<uint> Interval;
typedef std::set<BuildConversion> ConversionSet;

	Chain();
	Chain(int id, long score, int lLength, int lOffset, 
			const char *chrom, int rLength, int rOffset);
	Chain(const Chain& orig);
	virtual ~Chain();


	/**
	 * This is the entire chunk of data from a chain file associated
	 * to a single chain. Each component will be separated by new lines (\n)
	 */
	std::string Parse(const char *data);

	void EstimateConversion(int start, int stop, ConversionSet& conversionOptions);

	bool EstimateConversion(int start, int stop, BuildConversion& c);

	void EstimateConversion(int pos, ConversionSet& conversionOptions);
	bool EstimateConversion(int pos, BuildConversion& c);
	
	int LastPosition() { 
		return lOffset + lLength * direction;
	}


	int LLength();					///< Local length
	int RLength();					///< Remote Length

	int LStart();					///< Local start
	int RStart();					///< Remote start (this will be flipped if it's on - strain)

	int ID();						///< Chain ID
	int Direction();				///< direction of the remote segment
	int Score();					///< Chain's score
	
	std::string LChrom();		///< Local Chromosome
	std::string RChrom();		///< Remote chromosome

	int Count();					///< Return the number of subchains
private:
	RegionMap regions;			///< This represents the chaining data

	std::string lChrom;			///< Which chromosome we are starting from
	std::string rChrom;			///< Which chromosome this converts to
	int direction;					///< When -1, the SNPs in the region are transposed in the remote build

	int id;							///< ID associated with chain from the chain file
	/**
	 * Score is a percentage of the chain's score based on
	 * the amount of coverage a region represents (i.e. if we return a fragment
	 * that is 10% of the whole chain, we return score * 0.1
	 */
	long score;						///< Used to rank the conversions when multiple chains are hit
	int lLength;					///< Length of the local chain
	int lOffset;					///< Offset from beginning of the local chromosome
	
	int rLength;					///< Length of the remote chain
	int rOffset;					///< Offset from the beginning of the target chromosome
};


}

#endif	/* CHAIN_H */

