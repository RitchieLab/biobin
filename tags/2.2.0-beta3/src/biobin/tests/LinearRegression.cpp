/*
 * LinearRegression.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#include "LinearRegression.h"

#include "detail/MatrixUtils.h"

#include <limits>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_errno.h>


using std::string;
using std::vector;

using BioBin::Utility::Phenotype;

namespace BioBin {

namespace Test {

string LinearRegression::testname = LinearRegression::doRegister("linear");

double LinearRegression::runTest(const Bin& bin) const{

	if (_willfail) {
		return 1;
	}

	for(unsigned int i=0; i<_samp_name.size(); i++){
		gsl_matrix_set(_data, i, _data->size2 - 1,
				_pop_mgr_ptr->getTotalIndivContrib(bin, _samp_name[i].second, *_pheno_ptr));
	}

	// run the model now
	Result* r = calculate(*_phenos, *_data);

	double pval = 1;
	if(r && r->_conv){
	// Get the p-value of the last term
		double se = sqrt(gsl_matrix_get(r->cov, r->cov->size1 - 1, r->cov->size2 - 1));
		double c = gsl_vector_get(r->beta, r->beta->size - 1);
		pval = 2*gsl_cdf_tdist_Q(fabs(c / se),_data->size1 - _data->size2 + 1);

		delete r;
	}

	return pval;
}

float LinearRegression::getPhenotype(const PopulationManager& pop_mgr,
		const Utility::Phenotype& pheno, const std::string& samp) const{
	return pop_mgr.getPhenotypeVal(samp, pheno);
}

Regression::Result* LinearRegression::calculate(const gsl_vector& Y, const gsl_matrix& X) const{

	// these are the variables I need
	gsl_vector* beta = gsl_vector_alloc(X.size2);
	gsl_vector_set_all(beta, std::numeric_limits<double>::quiet_NaN());
	gsl_matrix* cov_mat = gsl_matrix_calloc(X.size2, X.size2);
	double chisq;
	Result* r = new Result(beta, cov_mat);
	int errcode = GSL_SUCCESS;

	gsl_permutation* permu = gsl_permutation_alloc(X.size2);

	// permutation matrix checking for colinearity
	unsigned int n_colinear = MatrixUtils::checkColinear(&X, permu);

	if(n_colinear == static_cast<unsigned int>(-1)){
		// if # colinear columns is -1, then something BAD happened while checking
		// for colinearity, set the errcode, please!
		errcode |= GSL_EFAILED;
	} else if(n_colinear > 0){
		// If we have colinearity, find them
		gsl_vector* idx_vec = gsl_vector_alloc(X.size2);
		for(unsigned int i=0; i<idx_vec->size; i++){
			gsl_vector_set(idx_vec, i, i);
		}
		errcode |= gsl_permute_vector(permu, idx_vec);

		r->dropped_cols.reserve(n_colinear);
		for(unsigned int i=0; i<n_colinear; i++){
			r->dropped_cols.push_back(gsl_vector_get(idx_vec, idx_vec->size - i - 1));
		}

		// sort this from smallest to largest so I'll get a consistent
		// permutation vector when needed
		std::sort(r->dropped_cols.begin(), r->dropped_cols.end());
		gsl_vector_free(idx_vec);
	}

	unsigned int n_indep = X.size2 - n_colinear;

	// set A = X*P
	gsl_matrix* A = gsl_matrix_alloc(X.size1, X.size2);

	if (errcode == GSL_SUCCESS){
		errcode |= gsl_matrix_memcpy(A, &X);
		errcode |= MatrixUtils::applyPermutation(A, permu);
	}

	if(errcode == GSL_SUCCESS){

		// get a matrix view of the first n-#colinear columns
		gsl_matrix_const_view A_v = gsl_matrix_const_submatrix(A, 0,0, X.size1, n_indep);
		gsl_vector_view bv = gsl_vector_subvector(r->beta, 0, n_indep);

		gsl_matrix_view cov_mat_view = gsl_matrix_submatrix(r->cov, 0, 0, n_indep, n_indep);

		// run the regression now
		gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(X.size1, n_indep);
		errcode |= gsl_multifit_linear(&A_v.matrix, &Y, &bv.vector, &cov_mat_view.matrix, &chisq, ws);
		gsl_multifit_linear_free(ws);

		// set the residuals here
		if(errcode == GSL_SUCCESS){
			r->resid = gsl_vector_alloc(A_v.matrix.size1);
			errcode |= gsl_multifit_linear_residuals(&A_v.matrix, &Y, &bv.vector, r->resid);
		}
	}


	if(errcode == GSL_SUCCESS){
		// Note: to unpermute, multiply by P transpose!
		// Also, we need to unpermute both the rows AND columns of cov_mat
		// permute columns
		errcode |= MatrixUtils::applyInversePermutation(r->cov, permu, true);
		// permute rows
		errcode |=MatrixUtils::applyInversePermutation(r->cov, permu, false);
		errcode |= gsl_permute_vector_inverse(permu, r->beta);
	}

	r->chisq = chisq;

	// And free everything... no memory leaks please!!
	gsl_matrix_free(A);
	gsl_permutation_free(permu);

	if(errcode != GSL_SUCCESS){
		r->_conv = false;
	}

	return r;
}

}

}
