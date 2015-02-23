/*
 * SKATLogistic.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: jrw32
 */

#include "SKATLogistic.h"

#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

#include "detail/SKATUtils.h"

using std::string;

namespace BioBin{
namespace Test{

string SKATLogistic::testname = SKATLogistic::doRegister("SKAT-logistic");

SKATLogistic::~SKATLogistic(){
	if(_XT_X_inv){
		gsl_matrix_free(_XT_X_inv);
	}
	if(_resid){
		gsl_vector_free(_resid);
	}
	if(_resid_wt){
		gsl_vector_free(_resid_wt);
	}
}

void SKATLogistic::init(){
	_base_reg.setup(*_pop_mgr_ptr, *_pheno_ptr);

	if(_base_reg._willfail){
		std::cerr << "WARNING: Base Logistic regression will fail, so will SKAT logistic." << std::endl;
	}

	// get the residual vector from the null model
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data, 0,0,
			_base_reg._data->size1, _base_reg._data->size2-1);
	if(_resid){
		gsl_vector_free(_resid);
	}
	_resid = gsl_vector_alloc(_base_reg._data->size1);
	gsl_multifit_linear_residuals(&X_v.matrix, _base_reg._phenos,
			_base_reg._null_result->beta, _resid);

	// now, get the _XT_X_inv matrix
	// it's a little more complex than in the SKAT_linear - now we need to do a
	// weighted regression and get the covariance matrix from THAT
	if(_resid_wt){
		gsl_vector_free(_resid_wt);
	}

	_resid_wt = gsl_vector_alloc(_resid->size);
	for(unsigned int i=0; i<_resid->size; i++){
		double val = gsl_vector_get(_resid, i);
		gsl_vector_set(_resid_wt, i, val*(1-val));
	}

	if(_XT_X_inv){
		gsl_matrix_free(_XT_X_inv);
	}
	_XT_X_inv = gsl_matrix_calloc(_base_reg._null_result->cov->size1, _base_reg._null_result->cov->size2);

	double chisq;
	gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(_resid->size, _XT_X_inv->size2);
	gsl_vector* tmp_beta = gsl_vector_alloc(_XT_X_inv->size2);

	// note that this gives me exactly what I want in _XT_X_inv
	gsl_multifit_wlinear(&X_v.matrix, _resid_wt, _base_reg._phenos, tmp_beta, _XT_X_inv, &chisq, ws);

	gsl_vector_free(tmp_beta);
	gsl_multifit_linear_free(ws);

}

double SKATLogistic::runTest(const Bin& bin) const{

	// check for guaranteed failure to begin with...
	if(_base_reg._willfail){
		return 1;
	}

	// first things first, let's set up the genotype matrix

	gsl_matrix* GW;

	unsigned int n_snp = SKATUtils::getGenoWeights(*_pop_mgr_ptr, *_pheno_ptr,
			bin, _base_reg._samp_name, GW);
	if(n_snp == 0){
		// clean up and return 1
		// note: GW will come out unallocated!
		return 1;
	}

	// now, get the Q statistic, defined to be
	// r*GKG*r, with K == Weights
	// temporary vectors - we want to keep everyting BLAS lv. 2 at this point
	gsl_vector* tmp_nsnp = gsl_vector_calloc(GW->size2);
	gsl_vector* tmp_nind = gsl_vector_calloc(GW->size1);

	// At this point, I'm done with weights and the permuted genotype matrix!
	// as well as the permutation matrix

	// Now, tmp_nsnp = (GW)^T * resid
	gsl_blas_dgemv(CblasTrans, 1.0, GW, _resid, 0, tmp_nsnp);
	// tmp_nind = (GW) * tmp_nsnp2
	gsl_blas_dgemv(CblasNoTrans, 1.0, GW, tmp_nsnp, 0, tmp_nind);

	// and finally, Q = _resid^T * tmp_nsnp
	double Q;
	gsl_blas_ddot(_resid, tmp_nind, &Q);

	// now divide by var(residuals) and divide by 2
	// acutally, multiply by the inverse of the above for a hint of extra speed
	Q *= 0.5;

	// get rid of the temporary vectors
	gsl_vector_free(tmp_nsnp);
	gsl_vector_free(tmp_nind);

	// And now we get the matrix for calculating the p-value
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data,
			0,0,_base_reg._data->size1, _base_reg._data->size2 - 1);

	// First, set Z = Z * pi_1 (or in our parlance GW_w = GW * _resid_wt)
	// We do this by scaling each row appropriately
	gsl_matrix* GW_w = gsl_matrix_alloc(GW->size1, GW->size2);
	gsl_matrix_memcpy(GW_w, GW);
	for(unsigned int i=0; i<GW->size2; i++){
		gsl_vector_view GW_row = gsl_matrix_row(GW_w, i);
		gsl_vector_scale(&GW_row.vector, gsl_vector_get(_resid_wt, i));
	}

	// Now, calculate (GW) * (GW)^T - (GW)^T * X * (X^T * X)^(-1) * X^T * (GW)
	gsl_matrix* tmp_ss = gsl_matrix_calloc(GW->size2, GW->size2);
	gsl_matrix* tmp_nv = gsl_matrix_calloc(GW->size1, X_v.matrix.size2);
	gsl_matrix* tmp_sv = gsl_matrix_calloc(GW->size2,X_v.matrix.size2);
	gsl_matrix* tmp_sn = gsl_matrix_calloc(GW->size2, GW->size1);

	// 1st lets get Z^T*Z
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GW, GW_w, 0, tmp_ss);
	// also, be sure to scale the rows by 1/resid_wt[i]

	// tmp_nv = X * (X^T * X)^(-1)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &X_v.matrix, _XT_X_inv, 0, tmp_nv);
	// tmp_sv = GW^T * tmp_nv
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, GW_w, tmp_nv, 0, tmp_sv);
	// tmp_sn = tmp_sv * X^T
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tmp_sv, &X_v.matrix, 0, tmp_sn);
	// tmp_ss = tmp_ss - tmp_sn * (GW)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, tmp_sn, GW_w, 1, tmp_ss);

	gsl_matrix_free(tmp_sv);
	gsl_matrix_free(tmp_sn);
	gsl_matrix_free(tmp_nv);
	gsl_matrix_free(GW);
	gsl_matrix_free(GW_w);

	// get the p-value from tmp_ss and the Q statistic
	double pval = SKATUtils::getPvalue(Q, tmp_ss);

	gsl_matrix_free(tmp_ss);

	return pval;
}

}
}

