/*
 * SKATUtils.h
 *
 *  Created on: Feb 23, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_SKATUTILS_H
#define BIOBIN_TEST_SKATUTILS_H

#include <vector>
#include <utility>
#include <string>

#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "biobin/PopulationManager.h"
#include "biobin/util/Phenotype.h"


namespace BioBin {

namespace Test {

class SKATUtils {
public:
	~SKATUtils() {}

private:
	SKATUtils() {}

public:
	// gets the genotype matrix and weight vector, returning the number of
	// SNPs in the gentoype matrix.  (Note: you really should check to see
	// if this is == 0, b/c it will fail!)
	static unsigned int getGenoWeights(const PopulationManager& pop_mgr,
			const Utility::Phenotype& pheno,
			const boost::dynamic_bitset<>& incl,
			const Bin& bin,
			const std::vector<std::pair<std::string, unsigned int> >& name_pos,
			gsl_matrix* &geno_wt);

	static double getPvalue(double Q, const gsl_matrix* W, double *accuracy);

        // configurable p-value calculation settings
        static double skat_matrix_threshold;
        static double skat_eigen_threshold;
        static double skat_pvalue_accuracy;
	static bool skat_raw_pvalues;

private:
	static unsigned char popcount(unsigned char v){
		v = v - ((v >> 1) & 0x55);                // put count of each 2 bits into those 2 bits
		v = (v & 0x33) + ((v >> 2) & 0x33); // put count of each 4 bits into those 4 bits
		return (v + ((v >> 4) & 0x0F)); // add the counts of each 4 bit
	}

	static boost::mutex _qfc_lock;
};

}

}

#endif /* SKATUTILS_H_ */
