/*
 * Regression.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: jrw32
 */

#include "Regression.h"
#include "MatrixUtils.h"

#include <vector>
#include <iostream>

#include <gsl/gsl_blas.h>

using std::vector;

using BioBin::PopulationManager;
using BioBin::Utility::Phenotype;


namespace BioBin {
namespace Test {

Regression::~Regression() {
	if(_data){
		gsl_matrix_free(_data);
	}

	if(_phenos){
		gsl_vector_free(_phenos);
	}

	if(_null_result){
		delete _null_result;
	}
}

void Regression::regressionSetup(const PopulationManager& pop_mgr, const Phenotype& pheno){

	// set up matrix of non-missing covariates (note the +2 is for the intercept + bin)

	gsl_matrix* data_tmp = gsl_matrix_alloc(pop_mgr.getNumSamples(), pop_mgr.getNumCovars() + 2);
	gsl_vector* pheno_tmp = gsl_vector_alloc(pop_mgr.getNumSamples());
	_included.resize(pop_mgr.getNumSamples(),false);

	unsigned int i=0;
	unsigned int s_idx=0;

	PopulationManager::const_sample_iterator s_end = pop_mgr.endSample();

	for(PopulationManager::const_sample_iterator si = pop_mgr.beginSample();
			si != s_end; si++){
		bool missing=false;
		float status=getPhenotype(pop_mgr, pheno, *si);
		const vector<float>& covars(pop_mgr.getCovariates(*si));

		// Note the empty loop here; the exit statement will exit iff we get to
		// the end of the covariates OR we see a missing value (ie nan).  If we
		// see a missing value, missing will be TRUE
		for(vector<float>::const_iterator ci = covars.begin();
				ci!=covars.end() && !(missing |= isnan(*ci)); ci++);

		missing |= isnan(status);

		if(!missing){
			gsl_vector_set(pheno_tmp, i, status);
			gsl_matrix_set(data_tmp, i, 0, 1);

			for(unsigned int j=0; j<covars.size(); j++){
				gsl_matrix_set(data_tmp, i, j+1, covars[j]);
			}

			_samp_name.push_back(std::make_pair(*si, s_idx));
			_included.set(s_idx, true);
			++i;
		}
		++s_idx;

	}

	if(i == 0){
		std::cerr << "WARNING: All samples dropped due to missingness, regression will fail!" << std::endl;
		gsl_matrix_free(data_tmp);
		gsl_vector_free(pheno_tmp);
		_willfail = true;

		// bail out now!
		return;
	}

	// OK, now create the _data and _pheno vars and copy the first i rows
	_data = gsl_matrix_alloc(i, data_tmp->size2);
	_phenos = gsl_vector_alloc(i);

	gsl_matrix_const_view data_tmp_view = gsl_matrix_const_submatrix(data_tmp, 0, 0, i, data_tmp->size2);
	gsl_vector_const_view pheno_tmp_view = gsl_vector_const_subvector(pheno_tmp, 0, i);

	gsl_matrix_memcpy(_data, &data_tmp_view.matrix);
	gsl_vector_memcpy(_phenos, &pheno_tmp_view.vector);

	gsl_matrix_free(data_tmp);
	gsl_vector_free(pheno_tmp);

	// first, get a view of the data excluding the missing column
	gsl_matrix_const_view X_view = gsl_matrix_const_submatrix(_data, 0, 0, _data->size1, _data->size2-1);

	// test for more columns than rows - guaranteed to fail!
	if (_data->size1 < _data->size2){
		std::cerr << "WARNING: Not enough samples; regression will fail!" << std::endl;
		_willfail = true;
		return;
	}

	_null_result = calculate(*_phenos, X_view.matrix);

	// if we have colinearity in the null result, let's just drop those
	// columns entirely
	if(_null_result && _null_result->dropped_cols.size() > 0){

		// allocate a new matrix to be used for _data
		gsl_matrix* data_new = gsl_matrix_alloc(_data->size1, _data->size2 - _null_result->dropped_cols.size());

		// permute the columns in _data
		gsl_permutation* permu = MatrixUtils::getPermutation(_null_result->dropped_cols, _data->size2);
		MatrixUtils::applyPermutation(_data, permu);

		// get a matrix view of the 1st collumns of _data (corresponding to dimensions of data_new)
		gsl_matrix_const_view indep_view = gsl_matrix_const_submatrix(_data, 0, 0, data_new->size1, data_new->size2);
		gsl_matrix_memcpy(data_new, &indep_view.matrix);

		// switch out the data_new and _data pointers
		std::swap(data_new, _data);

		// free the data_new structure (which was the "old" _data)
		gsl_matrix_free(data_new);
		gsl_permutation_free(permu);
		_null_result->dropped_cols.clear();
	}

	// If the null result didn't converge, something went horribly wrong!
	if(!_null_result || !_null_result->_conv){
		_willfail = true;
	}

}

}

}
