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
	LogisticRegression() : TestImpl<LogisticRegression>(testname), Regression(),
		_willfail(false) {
	}

	virtual ~LogisticRegression();

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

	// set this in the init if we know that we will fail for some reason
	bool _willfail;
};

}

}

#endif /* LOGISTICREGRESSION_H_ */
