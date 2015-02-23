/*
 * SKATLinear.cpp
 *
 *  Created on: Feb 18, 2015
 *      Author: jrw32
 */

#include "SKATLinear.h"

#include "detail/MatrixUtils.h"
#include "detail/SKATUtils.h"

#include <vector>
#include <set>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

using std::string;
using std::vector;
using std::set;

namespace BioBin {
namespace Test {

string SKATLinear::testname = SKATLinear::doRegister("SKAT-linear");

SKATLinear::~SKATLinear() {
	if(_resid){
		gsl_vector_free(_resid);
	}

	if(_XT_X_inv){
		gsl_matrix_free(_XT_X_inv);
	}

}

void SKATLinear::init(){
	_base_reg.setup(*_pop_mgr_ptr, *_pheno_ptr);


	// get the residual vector from the null model
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data, 0,0,
			_base_reg._data->size1, _base_reg._data->size2-1);
	if(_resid){
		gsl_vector_free(_resid);
	}
	_resid = gsl_vector_alloc(_base_reg._data->size1);
	gsl_multifit_linear_residuals(&X_v.matrix, _base_reg._phenos,
			_base_reg._null_result->beta, _resid);

	// now, get the _XT_X_inv matrix from the "cov" matrix of the null regression
	if(_XT_X_inv){
		gsl_matrix_free(_XT_X_inv);
	}
	_XT_X_inv = gsl_matrix_calloc(_base_reg._null_result->cov->size1, _base_reg._null_result->cov->size2);

	// calculate the (inverse variance) of the residuals
	resid_inv_var = 1 / gsl_stats_variance(_resid->data, _resid->stride, _resid->size);
	// set _XT_X_inv = inv_var * cov
	gsl_matrix* I = gsl_matrix_alloc(_XT_X_inv->size1, _XT_X_inv->size2);
	gsl_matrix_set_identity(I);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, resid_inv_var, _base_reg._null_result->cov, I, 0.0, _XT_X_inv);
	gsl_matrix_free(I);
}

double SKATLinear::runTest(const Bin& bin) const{
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
	Q *= 0.5*resid_inv_var;

	// get rid of the temporary vectors
	gsl_vector_free(tmp_nsnp);
	gsl_vector_free(tmp_nind);

	// And now we get the matrix for calculating the p-value
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data,
			0,0,_base_reg._data->size1, _base_reg._data->size2 - 1);

	// Now, calculate (GW) * (GW)^T - (GW)^T * X * (X^T * X)^(-1) * X^T * (GW)
	// or, written another way,
	gsl_matrix* tmp_ss = gsl_matrix_calloc(GW->size2, GW->size2);
	gsl_matrix* tmp_nv = gsl_matrix_calloc(GW->size1, X_v.matrix.size2);
	gsl_matrix* tmp_sv = gsl_matrix_calloc(GW->size2,X_v.matrix.size2);
	gsl_matrix* tmp_sn = gsl_matrix_calloc(GW->size2, GW->size1);

	// 1st lets get Z^T*Z
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GW, GW, 0, tmp_ss);

	// tmp_nv = X * (X^T * X)^(-1)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &X_v.matrix, _XT_X_inv, 0, tmp_nv);
	// tmp_sv = GW^T * tmp_nv
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, GW, tmp_nv, 0, tmp_sv);
	// tmp_sn = tmp_sv * X^T
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tmp_sv, &X_v.matrix, 0, tmp_sn);
	// tmp_ss = tmp_ss - tmp_sn * (GW)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, tmp_sn, GW, 1, tmp_ss);

	gsl_matrix_free(tmp_sv);
	gsl_matrix_free(tmp_sn);
	gsl_matrix_free(tmp_nv);
	gsl_matrix_free(GW);

	// get the p-value from tmp_ss and the Q statistic
	double pval = SKATUtils::getPvalue(Q, tmp_ss);

	gsl_matrix_free(tmp_ss);

	return pval;
}


}

}
