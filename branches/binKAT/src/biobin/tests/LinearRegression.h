/*
 * LinearRegression.h
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_LINEARREGRESSION_H
#define BIOBIN_TEST_LINEARREGRESSION_H

#include "Test.h"
#include "detail/Regression.h"

namespace BioBin {

namespace Test {

class LinearRegression : public TestImpl<LinearRegression>, public Regression {
	friend class SKATLinear;

public:
	LinearRegression() : TestImpl<LinearRegression>(testname), Regression() {}
	virtual ~LinearRegression() {}

	virtual Test* clone() const {return new LinearRegression(*this);}

protected:
	virtual void init(){regressionSetup(*_pop_mgr_ptr, *_pheno_ptr);}
	virtual double runTest(const Bin& bin) const;

	virtual Regression::Result* calculate(const gsl_vector& Y, const gsl_matrix& X) const;
	virtual float getPhenotype(const PopulationManager& pop_mgr,
				const Utility::Phenotype& pheno, const std::string& samp) const;

private:
	static std::string testname;

};

}

}

#endif /* LINEARREGRESSION_H_ */
