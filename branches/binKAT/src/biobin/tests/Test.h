/*
 * Test.h
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_TEST_H
#define BIOBIN_TEST_TEST_H

#include <string>

#include "biobin/Bin.h"
#include "biobin/PopulationManager.h"

#include "biobin/util/Phenotype.h"

#include "TestFactory.h"

namespace BioBin {

namespace Test{

/*!
 * A generic base class for running any stistical test, given a bin (and the
 * PopulationManager to get the phenotype/covariate/genetic data, of course)
 */
class Test {
public:
	Test() {}
	virtual ~Test() {}

	virtual const std::string& getName() const  = 0;

	template<class Bin_ptr_cont, class Pval_cont>
	void runAllTests(const PopulationManager& pop_mgr,
			const Utility::Phenotype& pheno,
			const Bin_ptr_cont& bins, Pval_cont& pvals_out);

protected:
	virtual void init(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno) = 0;
	virtual double runTest(const Bin& bin) const = 0;
};

template <class T>
class TestImpl : public virtual Test {
public:
	TestImpl(const std::string& name)
		: Test(), _name(name) {};

	static Test* create(){return new T();}

	virtual const std::string& getName() const {return _name;}

protected:
	static const std::string& doRegister(const std::string& key_in);

protected:
	std::string _name;


};

template<class Bin_ptr_cont, class Pval_cont>
void Test::runAllTests(const PopulationManager& pop_mgr,
		const Utility::Phenotype& pheno, const Bin_ptr_cont& bins,
		Pval_cont& pvals_out) {
	init(pop_mgr, pheno);
	pvals_out.clear();
	typename Bin_ptr_cont::const_iterator bin_itr = bins.begin();
	while(bin_itr != bins.end()){
		pvals_out.push_back(runTest(**bin_itr));
		++bin_itr;
	}
}

template<typename T>
const std::string& TestImpl<T>::doRegister(const std::string& key_in){
	return TestFactory::getFactory().RegisterTest(key_in, &T::create);
}

}
}

#endif /* TEST_H_ */
