/* 
 * File:   converter.h
 * Author: torstees
 *
 * Created on May 6, 2011, 1:25 PM
 */

#ifndef CONVERTER_H
#define	CONVERTER_H

#include <vector>
#include <string>
#include <ostream>

#include "snp.h"
#include "chain.h"

namespace LiftOver {

class Converter {
public:
	typedef std::multimap<uint, Chain> ChainMap;

	Converter(const char *origBuild, const char *newBuild);
	Converter(const Converter& orig);
	virtual ~Converter();

	void AddChain(const char *chainDetails);
	void ConvertDataset(const char *mapFilename, std::multimap<Utility::Locus, Utility::Locus>& converted);
	void ConvertDataset(const SnpArray& originalBuild, std::multimap<Utility::Locus, Utility::Locus>& converted);
	void ConvertDataset(const SnpArray& originalBuild, SnpArray& newBuild, std::ostream& droppedLog);
	
	void ConvertDataset(const char *mapFilename, SnpArray& newBuild, std::ostream& droppedLog);

protected:
	ChainMap chains;					///< Chrom -> chains
	std::string origBuild;
	std::string newBuild;
};



}

#endif	/* CONVERTER_H */

