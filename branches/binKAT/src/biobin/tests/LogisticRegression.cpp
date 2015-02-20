/*
 * LogisticRegression.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: jrw32
 */

#include "LogisticRegression.h"

namespace BioBin {

namespace Test {

LogisticRegression::~LogisticRegression() {
	// TODO Auto-generated destructor stub
}

void LogisticRegression::init(){
	if(_pheno_ptr->getStatus().first.count() == 0 ||
			_pheno_ptr->getStatus().second.count() == 0){
		std:cerr << "WARNING: Case or Control completely missing; "
				<< "Logistic regression will fail" << std::endl;
		_willfail = true;
	} else {
		regressionSetup(*_pop_mgr_ptr, *_pheno_ptr);
		if(_null_result->_conv){
			std:cerr << "WARNING: Logistic Regression null model diverged; "
							<< "Logistic regression will likely fail" << std::endl;
			_willfail = true;
		}
	}
}

double LogisticRegression::runTest(const Bin& bin) const {
	// If we are guaranteed to fail, bail out early
	if (_willfail) {
		return 1;
	}



	return 1;
}

}
