/*
 * LogisticRegression.h
 *
 *  Created on: Feb 20, 2015
 *      Author: jrw32
 */

#ifndef LOGISTICREGRESSION_H_
#define LOGISTICREGRESSION_H_

#include "Test.h"
#include "detail/Regression.h"

namespace BioBin {

namespace Test {

class LogisticRegression : public TestImpl<LogisticRegression> , public Regression{
public:
	LogisticRegression() : TestImpl<LogisticRegression>(testname), Regression(),
		_willfail(false) {
	}

	virtual ~LogisticRegression();

protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

private:
	static std::string testname;

	// set this in the init if we know that we will fail for some reason
	bool _willfail;
};

}

}

#endif /* LOGISTICREGRESSION_H_ */
