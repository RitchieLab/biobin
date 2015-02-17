/*
 * MatrixUtils.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#include "MatrixUtils.h"

#include <vector>
#include <set>
#include <limits>
#include <cmath>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using std::vector;
using std::fabs;
using std::set;

namespace BioBin{
namespace Test{

unsigned int MatrixUtils::checkColinear(const gsl_matrix* X, gsl_matrix* P){
	unsigned int n_rows = X->size1;
	unsigned int n_cols = X->size2;

	gsl_matrix* A = gsl_matrix_alloc(n_rows, n_cols);
	gsl_matrix* V = gsl_matrix_alloc(n_cols, n_cols);
	gsl_vector* S = gsl_vector_calloc(n_cols);
	gsl_vector* __ws_v = gsl_vector_alloc(n_cols);
	gsl_matrix* __ws_m = gsl_matrix_alloc(n_cols, n_cols);

	gsl_matrix_memcpy(A, X);

	gsl_linalg_SV_decomp_mod (A, __ws_m, V, S, __ws_v);

	// OK, now we look for all singular values less than machine epsilon
	// (we'll use float epsilon for double precision - lots of willge room there)
	unsigned int n_indep = n_cols;
	set<unsigned int> idx_set;
	while(gsl_vector_get(S,n_indep - 1) < std::numeric_limits<float>::epsilon()){
		// If we're here, the "n_indep - 1" column of V is a linear combination
		// of columns that adds to 0, so we want the last non-zero coefficient

		unsigned int p_idx = n_cols;

		// At the end of this (empty) loop, p_idx will be the index of the column
		// to drop
		while(--p_idx > 0){
			double val = fabs(gsl_matrix_get(V,p_idx,n_indep-1));
			if(val > std::numeric_limits<float>::epsilon() && idx_set.count(p_idx) == 0){
				break;
			}
		}

		idx_set.insert(p_idx);
		--n_indep;
	}

	vector<unsigned int> idx_permu(idx_set.begin(), idx_set.end());

	gsl_matrix* _P_ret = gsl_matrix_alloc(n_cols,n_cols);

	// P is our permutation matrix
	gsl_matrix_set_identity(_P_ret);

	if(n_indep != n_cols){

		// OK, let's construct a permutation matrix by walking from the front of
		// the idx_permu vector and creating a temporary permutation matrix
		// and left-multiplying by P (which is now the identity)
		gsl_matrix* P_tmp = gsl_matrix_alloc(n_cols, n_cols);
		gsl_matrix* _P_work = gsl_matrix_calloc(n_cols, n_cols);

		for(unsigned int i=0; i<idx_permu.size(); i++){
			unsigned int p_idx = idx_permu[i];
			gsl_matrix_set_identity(P_tmp);

			// Now, we want to transpose the (p_idx) and (n_cols - i - 1) columns
			// so we set (p_idx, p_idx) = (n_cols-i-1,n_cols-i-1) = 0
			// and (p_idx,n_cols-i-1) = (n_cols-i-1,p_idx) = 1
			gsl_matrix_set(P_tmp, p_idx,p_idx, 0);
			gsl_matrix_set(P_tmp, n_cols-i-1,n_cols-i-1, 0);
			gsl_matrix_set(P_tmp, n_cols-i-1,p_idx, 1);
			gsl_matrix_set(P_tmp, p_idx,n_cols-i-1, 1);

			// Now, set P = P * P_tmp
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, _P_ret, P_tmp, 0.0, _P_work);
			std::swap(_P_ret, _P_work);
		}
		gsl_matrix_free(P_tmp);
		gsl_matrix_free(_P_work);
	}

	// copy _P_ret into our return variable
	gsl_matrix_memcpy(P, _P_ret);
	gsl_matrix_free(_P_ret);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(__ws_v);
	gsl_matrix_free(__ws_m);

	return n_cols - n_indep;
}


}
}
