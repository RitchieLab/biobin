/* 
 * File:   chain.h
 * Author: torstees
 *
 * Represents a complete chain as described in the liftover chain files
 * A chromosome can have many chains.
 *
 * Created on April 28, 2011, 4:08 PM
 */

#ifndef KNOWLEDGE_CHAIN_H
#define	KNOWLEDGE_CHAIN_H

#include <boost/icl/interval_set.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <string>

#include <set>

#include "RegionConverter.h"

using std::set;
using boost::icl::interval;

namespace Knowledge {
namespace Liftover {

class BuildConversion;

class Chain {
public:

typedef boost::icl::split_interval_map<uint , set<RegionConverter> > RegionMap;

	Chain();
	Chain(int id, long score, int lLength, int lOffset, 
			const string& chrom, int rLength, int rOffset);
	virtual ~Chain() {}


	/**
	 * This is the entire chunk of data from a chain file associated
	 * to a single chain. Each component will be separated by new lines (\n)
	 */
	std::string Parse(const string& data);

	// Some of these should be private!!!
	void EstimateConversion(int start, int stop, set<BuildConversion>& conversionOptions);
	bool EstimateConversion(int start, int stop, BuildConversion& c);
	void EstimateConversion(int pos, set<BuildConversion>& conversionOptions);
	bool EstimateConversion(int pos, BuildConversion& c);
	
	int LastPosition() const{
		return _lOffset + _lLength * Direction();
	}


	int LLength() const {return _lLength;}					///< Local length
	int RLength() const {return _rLength;}					///< Remote Length
	int LStart() const {return _lOffset;}					///< Local start
	int RStart() const {return _rOffset;}					///< Remote start (this will be flipped if it's on - strain)
	int ID() const {return _id;}						///< Chain ID
	int Direction() const {return _forward ? 1 : -1;}				///< direction of the remote segment
	int Score() const {return _score;}					///< Chain's score
	const string& LChrom() const {return _lChrom;};		///< Local Chromosome
	const string& RChrom() const {return _rChrom;};		///< Remote chromosome

	int Count();					///< Return the number of subchains
private:
	// No copying or assingment!
	Chain(const Chain& orig);
	Chain& operator=(const Chain& orig);

	RegionMap _regions;			///< This represents the chaining data

	string _lChrom;			///< Which chromosome we are starting from
	string _rChrom;			///< Which chromosome this converts to
	bool _forward;					///< When -1, the SNPs in the region are transposed in the remote build

	int _id;							///< ID associated with chain from the chain file
	/**
	 * Score is a percentage of the chain's score based on
	 * the amount of coverage a region represents (i.e. if we return a fragment
	 * that is 10% of the whole chain, we return score * 0.1
	 */
	long _score;						///< Used to rank the conversions when multiple chains are hit
	int _lLength;					///< Length of the local chain
	int _lOffset;					///< Offset from beginning of the local chromosome
	
	int _rLength;					///< Length of the remote chain
	int _rOffset;					///< Offset from the beginning of the target chromosome
};

}// namespace Liftover
}// namespace Knowledge

#endif	/* CHAIN_H */

