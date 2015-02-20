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
	LinearRegression();
	virtual ~LinearRegression();

private:
	LinearRegression(const LinearRegression&);
	LinearRegression& operator=(const LinearRegression&);

protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

	virtual Regression::Result* calculate(const gsl_vector& Y, const gsl_matrix& X) const;

private:
	static std::string testname;

};

}

}

#endif /* LINEARREGRESSION_H_ */
