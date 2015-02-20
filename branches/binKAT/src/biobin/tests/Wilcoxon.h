/*
 * Wilcoxon.h
 *
 *  Created on: Feb 19, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_WILCOXON_H
#define BIOBIN_TEST_WILCOXON_H

#include "Test.h"

namespace BioBin {

namespace Test {

class Wilcoxon : public TestImpl<Wilcoxon> {
public:
	Wilcoxon() : TestImpl<Wilcoxon>(testname){}
	virtual ~Wilcoxon() {}

protected:
	virtual void init();
	virtual double runTest(const Bin& bin) const;

private:
	static std::string testname;
};

}

}

#endif /* WILCOXON_H_ */
