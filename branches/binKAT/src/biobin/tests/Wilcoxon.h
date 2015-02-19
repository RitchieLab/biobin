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
	Wilcoxon() : TestImpl<Wilcoxon>(testname), _pop_mgr_ptr(0), _pheno_ptr(0) {}
	virtual ~Wilcoxon() {}

public:
	virtual void init(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno);
	virtual double runTest(const Bin& bin) const;

private:
	static std::string testname;

	const PopulationManager* _pop_mgr_ptr;
	const Utility::Phenotype* _pheno_ptr;
};

}

}

#endif /* WILCOXON_H_ */
