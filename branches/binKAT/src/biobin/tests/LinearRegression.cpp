/*
 * LinearRegression.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#include "LinearRegression.h"

#include "MatrixUtils.h"

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
		pop_mgr_ptr(0), _pheno_ptr(0), _data(0), _phenos(0), _null_result(0) {
	// TODO Auto-generated constructor stub

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

void LinearRegression::init(const PopulationManager& pop_mgr, const Phenotype& pheno){
	pop_mgr_ptr = &pop_mgr;
	_pheno_ptr = &pheno;

	// set up matrix of non-missing covariates (note the +2 is for the intercept + bin)

	gsl_matrix* data_tmp = gsl_matrix_alloc(pop_mgr.getNumSamples(), pop_mgr.getNumCovars() + 2);
	gsl_vector* pheno_tmp = gsl_vector_alloc(pop_mgr.getNumSamples());

	unsigned int i=0;
	unsigned int s_idx=0;

	for(PopulationManager::const_sample_iterator si = pop_mgr.beginSample(); si != pop_mgr.endSample(); si++){
		bool missing=false;
		float status=pop_mgr.getPhenotypeVal(*si, pheno);
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
			++i;
		}
		++s_idx;

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

	_null_result = calculate(*_phenos, X_view.matrix);
}

double LinearRegression::runTest(const Bin& bin) const{

	for(unsigned int i=0; i<_samp_name.size(); i++){
		gsl_matrix_set(_data, i, _data->size2 - 1,
				pop_mgr_ptr->getTotalIndivContrib(bin, _samp_name[i].second, *_pheno_ptr));
	}

	// run the model now
	Result* r = calculate(*_phenos, *_data);

	// Get the p-value of the last term
	double se = sqrt(gsl_matrix_get(r->cov, r->cov->size1 - 1, r->cov->size2 - 1));
	double c = gsl_vector_get(r->beta, r->beta->size - 1);
	double pval = 2*gsl_cdf_tdist_Q(fabs(c / se),_data->size1 - _data->size2 + 1);

	gsl_matrix_free(r->cov);
	gsl_vector_free(r->beta);
	delete r;

	return pval;
}

LinearRegression::Result* LinearRegression::calculate(const gsl_vector& Y, const gsl_matrix& X){

	// these are the variables I need
	gsl_vector* beta = gsl_vector_alloc(X.size2);
	gsl_vector_set_all(beta, std::numeric_limits<double>::quiet_NaN());
	gsl_matrix* cov_mat = gsl_matrix_calloc(X.size2, X.size2);
	double chisq;
	Result* r = new Result(beta, cov_mat);

	// permutation matrix checking for colinearity
	gsl_matrix* P = gsl_matrix_alloc(X.size2, X.size2);

	unsigned int n_colinear = MatrixUtils::checkColinear(&X, P);
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
