/*
 * SKATLinear.cpp
 *
 *  Created on: Feb 18, 2015
 *      Author: jrw32
 */

#include "SKATLinear.h"

#include "MatrixUtils.h"

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

	if(_XT_X_inv){
		gsl_matrix_free(_XT_X_inv);
	}

}

void SKATLinear::init(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno){
	_base_reg.init(pop_mgr, pheno);

	_pop_mgr_ptr = &pop_mgr;
	_pheno_ptr = &pheno;

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
	double inv_var = 1 / gsl_stats_variance(_resid->data, _resid->stride, _resid->size);
	// set _XT_X_inv = inv_var * cov
	gsl_matrix* I = gsl_matrix_alloc(_XT_X_inv->size1, _XT_X_inv->size2);
	gsl_matrix_set_identity(I);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, inv_var, _base_reg._null_result->cov, I, 0.0, _XT_X_inv);
	gsl_matrix_free(I);
}

double SKATLinear::runTest(const Bin& bin) const{
	// first things first, let's set up the genotype matrix
	gsl_matrix* geno = gsl_matrix_alloc(_resid->size, bin.getVariantSize());
	gsl_matrix* weights = gsl_matrix_calloc(geno->size2, geno->size2);

	// get the average genotype (respecting the encoding)
	vector<float> avg_geno;
	avg_geno.reserve(bin.getVariantSize());
	Bin::const_locus_iterator ci=bin.variantBegin();
	unsigned int idx=0;
	while(ci != bin.variantEnd()){
		avg_geno.push_back(_pop_mgr_ptr->getAvgGenotype(**ci));
		gsl_matrix_set(weights, idx, idx, _pop_mgr_ptr->getLocusWeight(**ci, *_pheno_ptr));
		++ci;
		++idx;
	}

	// set up the amount of missingness for each SNP
	vector<unsigned int> missing(avg_geno.size());

	// vector ued to track the amount of variation in each SNP
	vector<unsigned char> n_genos(avg_geno.size());

	ci = bin.variantBegin();
	for(unsigned int j=0; j<geno->size2; j++){
		for(unsigned int i=0; i<geno->size1; i++){
			unsigned char g = _pop_mgr_ptr->getIndivGeno(**ci,_base_reg._samp_name[i].second);
			if(g == _pop_mgr_ptr->missing_geno){
				missing[j]++;
				gsl_matrix_set(geno,i,j,avg_geno[j]);
			} else {
				n_genos[j] |= (1 << g);
				gsl_matrix_set(geno,i,j,g);
			}
		}
		++ci;
	}

	vector<unsigned int> bad_idx;
	// max of 5% missing - perhaps customizeable some day?
	float missing_thresh = 0.05 *geno->size1;
	// OK, now let's check the missingness and variation requirements
	for(unsigned int i=0; i<missing.size(); i++){
		if(missing[i] > missing_thresh){
			bad_idx.push_back(i);
		} else if( popcount(n_genos[i]) <= 0){
			bad_idx.push_back(i);
		}
	}

	gsl_matrix* P = gsl_matrix_alloc(geno->size2, geno->size2);
	MatrixUtils::getPermuMatrix(bad_idx, P);

	// OK, now move the "bad" rows to the end, please!
	gsl_matrix* G_P = gsl_matrix_alloc(geno->size1, geno->size2);
	// permute the columns of geno into G_P
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, geno, P, 0.0, G_P);
	gsl_matrix_const_view G_v = gsl_matrix_const_submatrix(G_P, 0, 0,
			G_P->size1, G_P->size2 - bad_idx.size());

	// I'm officially done with "geno" now
	gsl_matrix_free(geno);

	// permute both rows + columns of weights
	gsl_matrix* wt_tmp = gsl_matrix_alloc(weights->size1, weights->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, weights, P, 0.0, wt_tmp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P, wt_tmp, 0.0, weights);
	gsl_matrix_free(wt_tmp);

	gsl_matrix_const_view W_v = gsl_matrix_const_submatrix(weights, 0, 0,
			weights->size1 - bad_idx.size(), weights->size2 - bad_idx.size());

	// now, get the Q statistic, defined to be
	// r*GKG*r, with K == Weights
	// temporary vectors - we want to keep everyting BLAS lv. 2 at this point
	gsl_vector* tmp_nsnp = gsl_vector_calloc(G_v.matrix.size2);
	gsl_vector* tmp_nsnp2 = gsl_vector_calloc(G_v.matrix.size2);
	gsl_vector* tmp_nind = gsl_vector_calloc(G_v.matrix.size1);

	// multiplying from right to left
	// tmp_nsnp = G^T * r
	gsl_blas_dgemv(CblasTrans, 1.0, &G_v.matrix, _resid, 0, tmp_nsnp);
	// tmp_nsnp2 = W * tmp_nsnp
	gsl_blas_dgemv(CblasNoTrans, 1.0, &W_v.matrix, tmp_nsnp, 0, tmp_nsnp2);
	// tmp_nind = G * tmp_nsnp2
	gsl_blas_dgemv(CblasNoTrans, 1.0, &G_v.matrix, tmp_nsnp2, 0, tmp_nind);
	// and finally, Q = _resid^T * tmp_nsnp
	double Q;
	gsl_blas_ddot(_resid, tmp_nind, &Q);

	// get rid of the temporary vectors
	gsl_vector_free(tmp_nsnp);
	gsl_vector_free(tmp_nsnp2);
	gsl_vector_free(tmp_nind);

	// And now we get the matrix for calculating the p-value
	// first, we're going to need some space to do an SVD
	gsl_matrix* U = gsl_matrix_alloc(G_v.matrix.size1, G_v.matrix.size2);
	gsl_matrix* V = gsl_matrix_alloc(G_v.matrix.size2, G_v.matrix.size2);
	gsl_vector* S = gsl_vector_alloc(G_v.matrix.size2);

	gsl_vector* svd_work_v = gsl_vector_alloc(G_v.matrix.size2);
	gsl_matrix* svd_work_m = gsl_matrix_alloc(G_v.matrix.size2, G_v.matrix.size2);

	// Copy G_v into U in preparation for SVD
	gsl_matrix_memcpy(U, &G_v.matrix);

	// perform SVD
	gsl_linalg_SV_decomp_mod(U, svd_work_m, V, S, svd_work_v);

	// Now, calculate Z^T * Z -

	// given that, we now compute the quadratic form

	return 1;
}


}

}
