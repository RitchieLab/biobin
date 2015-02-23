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
		_XT_X_inv(0), _resid(0), _resid_wt(0) {}

	virtual ~SKATLogistic();

protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

private:
	gsl_matrix* _XT_X_inv;
	gsl_vector* _resid;
	// a vector of p*(1-p), where p is the predicted value of the residual
	gsl_vector* _resid_wt;

	LogisticRegression _base_reg;

	static std::string testname;

};

}

}

#endif /* SKATLOGISTIC_H_ */
