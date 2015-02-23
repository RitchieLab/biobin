/*
 * SKATLinear.h
 *
 *  Created on: Feb 18, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_SKATLINEAR_H
#define BIOBIN_TEST_SKATLINEAR_H

#include "Test.h"
#include "LinearRegression.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace BioBin {

namespace Test {

class SKATLinear : public TestImpl<SKATLinear> {
public:
	SKATLinear() : TestImpl<SKATLinear>(testname),
		_resid(0), resid_inv_var(1), _XT_X_inv(0){}
	virtual ~SKATLinear();

private:
	SKATLinear(const SKATLinear&);
	SKATLinear& operator=(const SKATLinear&);

protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

private:

	static std::string testname;

	LinearRegression _base_reg;

	gsl_vector* _resid;
	double resid_inv_var;

	// The matrix (X^T*X)^(-1), where X is the design matrix of covariates
	gsl_matrix* _XT_X_inv;

};

}

}

#endif /* SKATLINEAR_H_ */
