/*
 * LogisticRegression.h
 *
 *  Created on: Feb 20, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_LOGISTICREGRESSION_H
#define BIOBIN_TEST_LOGISTICREGRESSION_H

#include "Test.h"
#include "detail/Regression.h"

#include <boost/array.hpp>

namespace BioBin {

namespace Test {

class LogisticRegression : public TestImpl<LogisticRegression> , public Regression{
public:
	friend class SKATLogistic;

	LogisticRegression() : TestImpl<LogisticRegression>(testname), Regression(){
	}

	virtual ~LogisticRegression() {}

	virtual Test* clone() const { return new LogisticRegression(*this);}

protected:

	// Inherited from Test
	virtual void init();
	virtual double runTest(const Bin& bin) const;

	// Inherited from Regression
	virtual Regression::Result* calculate(const gsl_vector& Y, const gsl_matrix& X) const;
	virtual float getPhenotype(const PopulationManager& pop_mgr,
			const Utility::Phenotype& pheno, const std::string& samp) const;

private:
	static std::string testname;

	static boost::array<double, 4> linkFunction(double v);

};

}

}

#endif /* LOGISTICREGRESSION_H_ */
