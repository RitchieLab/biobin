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
		resid_inv_var(1), X_svd_U(0), X_svd_S(0), X_svd_V(0), _willfail(false){}

	virtual ~SKATLinear();

//	virtual Test* clone() const {return new SKATLinear();}


protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

private:

	static std::string testname;

	LinearRegression _base_reg;

	double resid_inv_var;

	gsl_matrix* X_svd_U;
	gsl_vector* X_svd_S;
	gsl_matrix* X_svd_V;

	bool _willfail;

};

}

}

#endif /* SKATLINEAR_H_ */
