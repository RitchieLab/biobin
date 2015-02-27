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

#include <gsl/gsl_linalg.h>
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

	// calculate the (inverse variance) of the residuals
	resid_inv_var = 1 / gsl_stats_variance(_resid->data, _resid->stride, _resid->size);

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
	gsl_matrix_memcpy(X_svd_U, &X_v.matrix);
	gsl_vector* svd_vec_ws = gsl_vector_alloc(X_v.matrix.size2);
	gsl_matrix* svd_mat_ws = gsl_matrix_alloc(X_v.matrix.size2, X_v.matrix.size2);
	gsl_linalg_SV_decomp_mod(X_svd_U, svd_mat_ws, X_svd_V, X_svd_S, svd_vec_ws);
	gsl_vector_free(svd_vec_ws);
	gsl_matrix_free(svd_mat_ws);
}

double SKATLinear::runTest(const Bin& bin) const{
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
	// r^T*GKG^T*r, with K == Weights
	// temporary vectors - we want to keep everyting BLAS lv. 2 at this point
	gsl_vector* tmp_nsnp = gsl_vector_calloc(GW->size2);

	// Now, tmp_nsnp = (GW)^T * resid
	gsl_blas_dgemv(CblasTrans, 1.0, GW, _resid, 0, tmp_nsnp);

	double Q;
	// taking t(tmp_nsnp) %*% tmp_nsnp gives:
	// resid^T * (GW) * (GW)^T * resid
	gsl_blas_ddot(tmp_nsnp, tmp_nsnp, &Q);

	// now divide by var(residuals) and divide by 2
	// acutally, multiply by the inverse of the above for a hint of extra speed
	Q *= 0.5*resid_inv_var;

	// get rid of the temporary vectors
	gsl_vector_free(tmp_nsnp);

	// And now we get the matrix for calculating the p-value
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data,
			0,0,_base_reg._data->size1, _base_reg._data->size2 - 1);

	// Now, calculate (GW) * (GW)^T - (GW)^T * X * (X^T * X)^(-1) * X^T * (GW)
	// or, written another way,
	gsl_matrix* tmp_ss = gsl_matrix_calloc(GW->size2, GW->size2);
	gsl_matrix* tmp_vs = gsl_matrix_calloc(X_v.matrix.size2,GW->size2);
	gsl_matrix* tmp_sv = gsl_matrix_calloc(GW->size2,X_v.matrix.size2);

	// 1st lets get Z^T*Z
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GW, GW, 0, tmp_ss);

	// and get (GW)^T * X
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GW, &X_v.matrix, 0, tmp_sv);

	// Now, solve X * M = (GW) using the SVD of X
	for(unsigned int i=0; i<GW->size2; i++){
		gsl_vector_const_view GW_col = gsl_matrix_const_column(GW, i);
		gsl_vector_view tmp_vs_col = gsl_matrix_column(tmp_vs, i);
		gsl_linalg_SV_solve(X_svd_U, X_svd_V, X_svd_S, &GW_col.vector, &tmp_vs_col.vector);
	}

	// now, tmp_ss = tmp_ss - tmp_sv * tmp_vs
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1, tmp_sv, tmp_vs, 1, tmp_ss);

	gsl_matrix_free(tmp_sv);
	gsl_matrix_free(tmp_vs);
	gsl_matrix_free(GW);

	// get the p-value from tmp_ss and the Q statistic
	double pval = SKATUtils::getPvalue(Q, tmp_ss);

	gsl_matrix_free(tmp_ss);

	return pval;
}


}

}
