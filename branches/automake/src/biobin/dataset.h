/* 
 * File:   dataset.h
 * Author: torstees
 *
 * Created on June 7, 2011, 10:41 AM
 * 
 * Convenience class 
 */

#ifndef DATASET_H
#define	DATASET_H

#include <boost/dynamic_bitset.hpp>
#include "utility/types.h"
#include "liftover/snp.h"
#include "liftover/chain.h"

namespace VCF {

class Dataset {
public:
	Dataset(Utility::StringArray& filenames);
	Dataset(const Dataset& orig);
	virtual ~Dataset();
	
	void RealizeDataset(uint totalIndividualCount, uint totalBinCount, uint totalGenotypeCount);
	
protected:
	
};



inline
Dataset::Dataset(Utility::StringArray& filenames) : filenames(filenames) {}

inline
Dataset::Dataset(const Dataset& orig) : filenames(orig.filenames) {}

inline
Dataset::~Dataset() {}




}
#endif	/* DATASET_H */

