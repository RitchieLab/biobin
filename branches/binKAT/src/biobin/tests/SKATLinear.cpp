/*
 * SKATLinear.cpp
 *
 *  Created on: Feb 18, 2015
 *      Author: jrw32
 */

#include "SKATLinear.h"

#include "detail/qfc.h"
#include "detail/MatrixUtils.h"

#include <vector>
#include <set>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

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
	gsl_matrix* geno = gsl_matrix_alloc(_resid->size, bin.getVariantSize());
	gsl_matrix* weights = gsl_matrix_calloc(geno->size2, geno->size2);

	// get the average genotype (respecting the encoding)
	vector<float> avg_geno;
	avg_geno.reserve(bin.getVariantSize());
	Bin::const_locus_iterator ci=bin.variantBegin();
	unsigned int idx=0;
	while(ci != bin.variantEnd()){
		avg_geno.push_back(_pop_mgr_ptr->getAvgGenotype(**ci));
		gsl_matrix_set(weights, idx, idx, _pop_mgr_ptr->getLocusWeight(**ci, *_pheno_ptr, bin.getRegion()));
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
			if(g > 2){
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
		// too much missingness
		if(missing[i] > missing_thresh){
			bad_idx.push_back(i);
		// not polymorphic
		} else if( popcount(n_genos[i]) <= 0){
			bad_idx.push_back(i);
		}
	}

	// We have no SNPS!
	if(bad_idx.size() == geno->size2){
		// clean up and return 1
		gsl_matrix_free(geno);
		gsl_matrix_free(weights);
		return 1;
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
	gsl_vector* tmp_nind = gsl_vector_calloc(G_v.matrix.size1);

	// Let's calculate the matrix G*W
	gsl_matrix* GW = gsl_matrix_calloc(G_v.matrix.size1, G_v.matrix.size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &G_v.matrix, &W_v.matrix, 0.0, GW);

	// At this point, I'm done with weights and the permuted genotype matrix!
	// as well as the permutation matrix
	gsl_matrix_free(weights);
	gsl_matrix_free(G_P);
	gsl_matrix_free(P);

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

	// Now, calculate (GW) * (GW)^T - (GW)^T * X * (X^T * X)^(-1) * X^T * (GW)^T
	// or, written another way,
	// (GW)^T * (I - X (X^T X)^-1 X^T) * (GW)^T
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

	// OK, now we have to take the eigenvalues of tmp_ss
	gsl_eigen_symm_workspace* eigen_w = gsl_eigen_symm_alloc(tmp_ss->size1);
	gsl_vector* eval = gsl_vector_alloc(tmp_ss->size1);
	gsl_eigen_symm(tmp_ss, eval, eigen_w);
	gsl_eigen_symm_free(eigen_w);

	// at this point, I'm doe with EVERYTHING except eval
	gsl_matrix_free(tmp_ss);
	gsl_matrix_free(GW);

	// Now, sort the eigenvalues in descending order
	std::sort(eval->data, eval->data + eval->size, std::greater<double>());
	// TODO: find where these eigenvalues go to zero and only take those
	int n_eval = eval->size;

	// I don't feel like doing memory management, so use a vector instead of an array
	std::vector<double> nct(n_eval, 0);
	std::vector<int> df(n_eval, 1);
	std::vector<double> qfc_detail(7);
	int qfc_err;
	double pval;
	int lim=10000;
	double acc=0.0001;
	double sigma=0;

	qfc(eval->data, &nct[0], &df[0], &n_eval, &sigma, &Q, &lim, &acc, &qfc_detail[0], &qfc_err, &pval);

	// clean up, please!
	gsl_vector_free(eval);

	return pval;
}


}

}
