/*
 * Regression.h
 *
 *  Created on: Feb 20, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_REGRESSION_H
#define BIOBIN_TEST_REGRESSION_H

#include <vector>
#include <string>
#include <utility>

#include <boost/dynamic_bitset.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "biobin/PopulationManager.h"
#include "biobin/util/Phenotype.h"

namespace BioBin {

namespace Test {

class Regression {
public:
	Regression() : _data(0), _phenos(0), _null_result(0), _willfail(false){
	}
	virtual ~Regression();

public:
	class Result{
	public:
		Result(gsl_vector* b, gsl_matrix* c) : beta(b), cov(c), resid(0), chisq(0), _conv(true) {}
		~Result(){
			if(beta){
				gsl_vector_free(beta);
			}
			if(cov){
				gsl_matrix_free(cov);
			}
			if(resid){
				gsl_vector_free(resid);
			}
		}

		std::vector<unsigned int> dropped_cols;
		gsl_vector* beta;
		gsl_matrix* cov;
		gsl_vector* resid;
		double chisq;
		bool _conv;
	};

//protected:

	void regressionSetup(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno);
	virtual Result* calculate(const gsl_vector& Y, const gsl_matrix& X) const = 0;
	virtual float getPhenotype(const PopulationManager& pop_mgr,
			const Utility::Phenotype& pheno, const std::string& samp) const = 0;

	//! The matrix of covariates + bin
	gsl_matrix* _data;

	//! The vector of phenotypes
	gsl_vector* _phenos;

	//! The null model result
	Result* _null_result;

	boost::dynamic_bitset<> _included;

	//! A vector of sample IDs containing completely non-missing covariates + phenotypes
	//! along with the index of each!
	std::vector<std::pair<std::string, unsigned int> > _samp_name;

	// set this in the init if we know that we will fail for some reason
	bool _willfail;

};

}

}

#endif /* REGRESSION_H_ */
