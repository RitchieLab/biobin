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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#include "detail/SKATUtils.h"

using std::string;

namespace BioBin{
namespace Test{

string SKATLogistic::testname = SKATLogistic::doRegister("SKAT-logistic");

SKATLogistic::~SKATLogistic(){
	if(_resid_wt){
		gsl_vector_free(_resid_wt);
	}
	if(X_svd_U){
		gsl_matrix_free(X_svd_U);
	}
	if(X_svd_S){
		gsl_vector_free(X_svd_S);
	}
	if(X_svd_V){
		gsl_matrix_free(X_svd_V);
	}
}

void SKATLogistic::init(){
	_base_reg.setup(*_pop_mgr_ptr, *_pheno_ptr);

	int errcode = GSL_SUCCESS;

	if(_base_reg._willfail){
		std::cerr << "WARNING: Base Logistic regression will fail, so will SKAT logistic." << std::endl;
		return;
	}

	// get the residual vector from the null model
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data, 0,0,
			_base_reg._data->size1, _base_reg._data->size2-1);

	// now, get the _XT_X_inv matrix
	// it's a little more complex than in the SKAT_linear - now we need to do a
	// weighted regression and get the covariance matrix from THAT
	if(_resid_wt){
		gsl_vector_free(_resid_wt);
	}

	// First, get the predicted value
	_resid_wt = gsl_vector_calloc(_base_reg._phenos->size);
	errcode |= gsl_vector_memcpy(_resid_wt, _base_reg._phenos);
	errcode |= gsl_blas_daxpy(-1, _base_reg._null_result->resid, _resid_wt);

	// now, val * (1-val) for each entry in _resid_wt
	for(unsigned int i=0; i<_resid_wt->size; i++){
		// val is the predicted value
		double val = gsl_vector_get(_resid_wt, i);
		gsl_vector_set(_resid_wt, i, val*(1-val));
	}

	if(X_svd_U){
		gsl_matrix_free(X_svd_U);
	}
	if(X_svd_S){
		gsl_vector_free(X_svd_S);
	}
	if(X_svd_V){
		gsl_matrix_free(X_svd_V);
	}
	// let's get the SVD of X here
	X_svd_U = gsl_matrix_alloc(X_v.matrix.size1, X_v.matrix.size2);
	X_svd_S = gsl_vector_alloc(X_v.matrix.size2);
	X_svd_V = gsl_matrix_alloc(X_v.matrix.size2, X_v.matrix.size2);
	errcode |= gsl_matrix_memcpy(X_svd_U, &X_v.matrix);

	// now, weight the X_svd_U matrix by sqrt(_resid_wt)
	for(unsigned int i=0; i<_resid_wt->size; i++){
		gsl_vector_view X_row = gsl_matrix_row(X_svd_U, i);
		errcode |= gsl_vector_scale(&X_row.vector, sqrt(gsl_vector_get(_resid_wt, i)));
	}

	// Now, perform the SVD
	gsl_vector* svd_vec_ws = gsl_vector_alloc(X_v.matrix.size2);
	gsl_matrix* svd_mat_ws = gsl_matrix_alloc(X_v.matrix.size2, X_v.matrix.size2);
	errcode |= gsl_linalg_SV_decomp_mod(X_svd_U, svd_mat_ws, X_svd_V, X_svd_S, svd_vec_ws);
	gsl_vector_free(svd_vec_ws);
	gsl_matrix_free(svd_mat_ws);

	if(errcode != GSL_SUCCESS){
		_willfail = true;
	}


}

double SKATLogistic::runTest(const Bin& bin) const{

	// check for guaranteed failure to begin with...
	if(_willfail || _base_reg._willfail){
		return 1;
	}

	int errcode = GSL_SUCCESS;

	// first things first, let's set up the genotype matrix

	gsl_matrix* GW;

	unsigned int n_snp = SKATUtils::getGenoWeights(*_pop_mgr_ptr, *_pheno_ptr, _base_reg._included,
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

	// Now, tmp_nsnp = (GW)^T * resid
	errcode |= gsl_blas_dgemv(CblasTrans, 1.0, GW, _base_reg._null_result->resid, 0, tmp_nsnp);

	double Q;
	// taking t(tmp_nsnp) %*% tmp_nsnp gives:
	// resid^T * (GW) * (GW)^T * resid
	errcode |= gsl_blas_ddot(tmp_nsnp, tmp_nsnp, &Q);

	// now divide by 2
	// acutally, multiply by the inverse of the above for a hint of extra speed
	Q *= 0.5;

	// get rid of the temporary vectors
	gsl_vector_free(tmp_nsnp);

	// And now we get the matrix for calculating the p-value
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data,
			0,0,_base_reg._data->size1, _base_reg._data->size2 - 1);

	// Now, calculate (GW) * (GW)^T - (GW)^T * X * (X^T * X)^(-1) * X^T * (GW)
	gsl_matrix* tmp_ss = gsl_matrix_calloc(GW->size2, GW->size2);
	gsl_matrix* tmp_vs = gsl_matrix_calloc(X_v.matrix.size2, GW->size2);
	gsl_matrix* tmp_sv = gsl_matrix_calloc(GW->size2, X_v.matrix.size2);

	// First, set Z = Z * pi_1 (or in our parlance GW_w = GW * _resid_wt)
	// We do this by scaling each row appropriately
	gsl_matrix* GW_w = gsl_matrix_alloc(GW->size1, GW->size2);
	errcode |= gsl_matrix_memcpy(GW_w, GW);
	for(unsigned int i=0; i<GW->size1; i++){
		gsl_vector_view GW_row = gsl_matrix_row(GW_w, i);
		errcode |= gsl_vector_scale(&GW_row.vector, gsl_vector_get(_resid_wt, i));
	}

	/*
	gsl_vector* gw_norm = gsl_vector_alloc(GW->size2);
	for (unsigned int i=0; i<GW->size2; i++){
		gsl_vector_const_view GW_col = gsl_matrix_const_column(GW, i);
		gsl_vector_set(gw_norm, i, gsl_blas_dnrm2(&GW_col.vector));
	}
	gsl_vector* gww_norm = gsl_vector_alloc(GW->size1);
	for (unsigned int i=0; i<GW->size2; i++){
		gsl_vector_const_view GW_col = gsl_matrix_const_column(GW_w, i);
		gsl_vector_set(gww_norm, i, gsl_blas_dnrm2(&GW_col.vector));
	}
	*/

	// 1st lets get Z^T*Z
	errcode |= gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, GW, GW_w, 0, tmp_ss);
	// and get (GW)^T * X
	errcode |= gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GW_w, &X_v.matrix, 0, tmp_sv);

	// Now, I'm going to need GW scaled by sqrt(_resid_wt) to use the SVD here
	errcode |= gsl_matrix_memcpy(GW_w, GW);
	for(unsigned int i=0; i<GW->size1; i++){
		gsl_vector_view GW_row = gsl_matrix_row(GW_w, i);
		errcode |= gsl_vector_scale(&GW_row.vector, sqrt(gsl_vector_get(_resid_wt, i)));
	}

	// Now, let's get (X^T X)^-1 X^T GW using the SVD
	for(unsigned int i=0; i<GW_w->size2; i++){
		gsl_vector_const_view GW_col = gsl_matrix_const_column(GW_w, i);
		gsl_vector_view tmp_vs_col = gsl_matrix_column(tmp_vs, i);
		errcode |= gsl_linalg_SV_solve(X_svd_U, X_svd_V, X_svd_S, &GW_col.vector, &tmp_vs_col.vector);
	}

	// now, tmp_ss = tmp_ss - tmp_sv * tmp_vs
	errcode |= gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, tmp_sv, tmp_vs, 1, tmp_ss);

	double pval;
	if(errcode == GSL_SUCCESS){
		// get the p-value from tmp_ss and the Q statistic
		pval = SKATUtils::getPvalue(Q, tmp_ss);
	} else {
		pval = 10;
	}

	gsl_matrix_free(tmp_sv);
	gsl_matrix_free(tmp_vs);
	gsl_matrix_free(GW);
	gsl_matrix_free(GW_w);

	gsl_matrix_free(tmp_ss);

	return pval;
}

}
}

