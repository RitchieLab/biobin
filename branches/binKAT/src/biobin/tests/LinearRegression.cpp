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

using std::string;
using std::vector;

using BioBin::Utility::Phenotype;

namespace BioBin {

namespace Test {

string LinearRegression::testname = LinearRegression::doRegister("linear");

LinearRegression::LinearRegression() : TestImpl<LinearRegression>(testname),
		Regression(){
}

LinearRegression::~LinearRegression() {
	if(_data){
		gsl_matrix_free(_data);
	}
	if(_phenos){
		gsl_vector_free(_phenos);
	}
	if(_null_result){
		gsl_matrix_free(_null_result->cov);
		gsl_vector_free(_null_result->beta);
		delete(_null_result);
	}
}

void LinearRegression::init(){

	regressionSetup(*_pop_mgr_ptr, *_pheno_ptr);
}

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

	// Get the p-value of the last term
	double se = sqrt(gsl_matrix_get(r->cov, r->cov->size1 - 1, r->cov->size2 - 1));
	double c = gsl_vector_get(r->beta, r->beta->size - 1);
	double pval = 2*gsl_cdf_tdist_Q(fabs(c / se),_data->size1 - _data->size2 + 1);

	delete r;

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

	// permutation matrix checking for colinearity
	gsl_matrix* P = gsl_matrix_alloc(X.size2, X.size2);

	unsigned int n_colinear = MatrixUtils::checkColinear(&X, P);
	// If we have colinearity, find them
	if(n_colinear > 0){
		gsl_vector* idx_vec = gsl_vector_alloc(X.size2);
		for(unsigned int i=0; i<idx_vec->size; i++){
			gsl_vector_set(idx_vec, i, i);
		}
		gsl_vector* idx_res = gsl_vector_calloc(idx_vec->size);
		gsl_blas_dgemv(CblasNoTrans, 1.0, P, idx_vec, 0, idx_res);
		r->dropped_cols.reserve(n_colinear);
		for(unsigned int i=0; i<n_colinear; i++){
			r->dropped_cols.push_back(gsl_vector_get(idx_res, idx_res->size - i - 1));
		}

		// sort this from smallest to largest so I'll get a consistent
		// permutation vector when needed
		std::sort(r->dropped_cols.begin(), r->dropped_cols.end());
		gsl_vector_free(idx_vec);
		gsl_vector_free(idx_res);
	}

	unsigned int n_indep = X.size2 - n_colinear;

	// set A = X*P
	gsl_matrix* A = gsl_matrix_alloc(X.size1, X.size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &X, P, 0.0, A);

	// get a matrix view of the first n-#colinear columns
	gsl_matrix_const_view A_v = gsl_matrix_const_submatrix(A, 0,0, X.size1, n_indep);
	gsl_vector_view bv = gsl_vector_subvector(beta, 0, n_indep);

	gsl_matrix_view cov_mat_view = gsl_matrix_submatrix(cov_mat, 0, 0, n_indep, n_indep);

	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(X.size1, n_indep);

	// run the regression now
	gsl_multifit_linear(&A_v.matrix, &Y, &bv.vector, &cov_mat_view.matrix, &chisq, ws);
	// set the residuals here
	r->resid = gsl_vector_alloc(A_v.matrix.size1);
	gsl_multifit_linear_residuals(&A_v.matrix, &Y, &bv.vector, r->resid);

	// Note: to unpermute, multiply by P transpose!
	// Also, we need to unpermute both the rows AND columns of cov_mat
	//bv = gsl_vector_view_array(r->beta_vec, n_cols);
	gsl_matrix* _cov_work = gsl_matrix_calloc(X.size2, X.size2);
	// permute columns
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, cov_mat, P, 0.0, _cov_work);
	// permute rows
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, _cov_work, 0.0, cov_mat);
	gsl_matrix_free(_cov_work);

	gsl_vector* _bv_work = gsl_vector_alloc(X.size2);
	gsl_vector_memcpy(_bv_work, beta);
	gsl_blas_dgemv(CblasTrans, 1.0, P, _bv_work, 0.0, beta);
	gsl_vector_free(_bv_work);

	r->chisq = chisq;

	// And free everything... no memory leaks please!!
	gsl_matrix_free(P);
	gsl_matrix_free(A);
	gsl_multifit_linear_free(ws);

	return r;
}

}

}
