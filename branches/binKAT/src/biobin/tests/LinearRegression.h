/*
 * LinearRegression.h
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_LINEARREGRESSION_H
#define BIOBIN_TEST_LINEARREGRESSION_H

#include <vector>
#include <string>
#include <utility>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "Test.h"


namespace BioBin {

namespace Test {

class LinearRegression : public TestImpl<LinearRegression>{
public:
	class Result{
	public:
		Result(gsl_vector* b, gsl_matrix* c) : beta(b), cov(c), chisq(0) {}

		gsl_vector* beta;
		gsl_matrix* cov;
		double chisq;
	};

public:
	LinearRegression();
	virtual ~LinearRegression();

private:
	LinearRegression(const LinearRegression&);
	LinearRegression& operator=(const LinearRegression&);

public:
	virtual void init(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno);
	virtual double runTest(const Bin& bin) const;

private:

	static Result* calculate(const gsl_vector& Y, const gsl_matrix& X);

	static std::string testname;

	const PopulationManager* pop_mgr_ptr;
	const Utility::Phenotype* _pheno_ptr;
	const Knowledge::Information* _info;

	//! The matrix of covariates + bin
	gsl_matrix* _data;

	//! The vector of phenotypes
	gsl_vector* _phenos;

	//! The null model result
	Result* _null_result;

	//! A vector of sample IDs containing completely non-missing covariates + phenotypes
	//! along with the index of each!
	std::vector<std::pair<std::string, unsigned int> > _samp_name;


};

}

}

#endif /* LINEARREGRESSION_H_ */
