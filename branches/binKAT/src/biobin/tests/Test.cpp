/*
 * Test.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#include "Test.h"
#include "biobin/PopulationManager.h"
#include "biobin/Bin.h"

namespace BioBin {

namespace Test{

void Test::setup(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno){
	_pop_mgr_ptr = &pop_mgr;
	_pheno_ptr = &pheno;
	init();
}

}

}
