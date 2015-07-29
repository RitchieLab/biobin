/*
 * SKATLogistic.h
 *
 *  Created on: Feb 23, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_SKATLOGISTIC_H
#define BIOBIN_TEST_SKATLOGISTIC_H

#include "Test.h"
#include "LogisticRegression.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace BioBin {
namespace Test {

class SKATLogistic : public TestImpl<SKATLogistic> {
public:
	SKATLogistic() : TestImpl<SKATLogistic>(testname),
		_resid_wt(0), X_svd_U(0), X_svd_S(0), X_svd_V(0) {}

	virtual ~SKATLogistic();

	virtual Test* clone() const {return new SKATLogistic(*this);}

protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

private:

	// a vector of p*(1-p), where p is the predicted value of the residual
	gsl_vector* _resid_wt;

	gsl_matrix* X_svd_U;
	gsl_vector* X_svd_S;
	gsl_matrix* X_svd_V;

	LogisticRegression _base_reg;

	static std::string testname;

};

}

}

#endif /* SKATLOGISTIC_H_ */
