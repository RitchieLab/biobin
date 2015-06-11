/*
 * SKATUtils.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: jrw32
 */

#include "SKATUtils.h"

#include <algorithm>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_permute.h>

#include "qfc.h"
#include "MatrixUtils.h"

using std::vector;
using std::string;
using std::pair;


namespace BioBin{
namespace Test{

boost::mutex SKATUtils::_qfc_lock;

unsigned int SKATUtils::getGenoWeights(const PopulationManager& pop_mgr,
		const Utility::Phenotype& pheno,
		const bm::bvector<>& incl,
		const Bin& bin,
		const vector<pair<string, unsigned int> >& name_pos,
		gsl_matrix* &geno){

	unsigned int n_col = bin.getVariantSize();
	unsigned int n_row = name_pos.size();

	gsl_matrix* geno_tmp = gsl_matrix_alloc(n_row, n_col);
	gsl_vector* wt_tmp = gsl_vector_alloc(n_col);

	// get the average genotype (respecting the encoding)
	vector<double> avg_geno;
	avg_geno.reserve(bin.getVariantSize());
	Bin::const_locus_iterator ci=bin.variantBegin();
	unsigned int idx=0;
	while(ci != bin.variantEnd()){
		avg_geno.push_back(pop_mgr.getAvgGenotype(**ci, &incl));
		gsl_vector_set(wt_tmp, idx, pop_mgr.getLocusWeight(**ci, pheno, bin.getRegion()));
		++ci;
		++idx;
	}

	// set up the amount of missingness for each SNP
	vector<unsigned int> missing(avg_geno.size());

	// vector ued to track the amount of variation in each SNP
	vector<unsigned char> n_genos(avg_geno.size());

	ci = bin.variantBegin();
	for(unsigned int j=0; j<n_col; j++){
		for(unsigned int i=0; i<n_row; i++){
			unsigned char g = pop_mgr.getIndivGeno(**ci,name_pos[i].second);
			if(g > 2){
				missing[j]++;
				gsl_matrix_set(geno_tmp,i,j,avg_geno[j]);
			} else {
				n_genos[j] |= (1 << g);
				gsl_matrix_set(geno_tmp,i,j,g);
			}
		}
		++ci;
	}

	unsigned int n_miss = 0;
	unsigned int n_mono = 0;

	vector<unsigned int> bad_idx;
	// max of 15% missing - perhaps customizeable some day?
	float missing_thresh = 0.15 * n_row;
	// OK, now let's check the missingness and variation requirements
	for(unsigned int i=0; i<missing.size(); i++){
		// too much missingness
		if(missing[i] > missing_thresh){
			++n_miss;
			bad_idx.push_back(i);
			for(unsigned int j=0; j<geno_tmp->size1; j++){
				gsl_matrix_set(geno_tmp, j, i, std::numeric_limits<double>::quiet_NaN());
			}
		// not polymorphic
		} else if( popcount(n_genos[i]) <= 1){
			++n_mono;
			bad_idx.push_back(i);
			for(unsigned int j=0; j<geno_tmp->size1; j++){
				gsl_matrix_set(geno_tmp, j, i, std::numeric_limits<double>::quiet_NaN());
			}
		}
	}

	// We have no SNPS!  bail out!
	if(bad_idx.size() == n_col){
		// clean up and return 1
		gsl_matrix_free(geno_tmp);
		gsl_vector_free(wt_tmp);
		geno = 0;
		return 0;
	} else if(bad_idx.size() > 0){
		// get rid of the "bad" indexes

		gsl_permutation* permu = MatrixUtils::getPermutation(bad_idx, n_col);
		MatrixUtils::applyPermutation(geno_tmp, permu, true);
		gsl_permute_vector(permu, wt_tmp);
		gsl_permute(permu->data, &avg_geno[0], 1, avg_geno.size());
		// done with permutation, please clean up!
		gsl_permutation_free(permu);
		gsl_matrix_const_view G_v = gsl_matrix_const_submatrix(geno_tmp, 0, 0,
			geno_tmp->size1, geno_tmp->size2 - bad_idx.size());

		geno = gsl_matrix_alloc(n_row, n_col - bad_idx.size());
		gsl_matrix_memcpy(geno, &G_v.matrix);
		gsl_matrix_free(geno_tmp);

	} else {
		// nothing more to do here, just return as-is!
		geno = geno_tmp;
	}

	// OK, now go through the columns of geno and scale them by the
	// corresponding weight.  What will return will be the matrix (GW)
	for(unsigned int i=0; i<wt_tmp->size - bad_idx.size(); i++){
		gsl_vector_view gc = gsl_matrix_column(geno, i);
		gsl_vector_scale(&gc.vector, gsl_vector_get(wt_tmp, i));
	}

	gsl_vector_free(wt_tmp);

	return n_col - bad_idx.size();
}

double SKATUtils::getPvalue(double Q, const gsl_matrix* W){
	gsl_matrix* W_tmp = gsl_matrix_calloc(W->size1, W->size2);

	// find columns (and rows) that are essentially 0 to remove them
	// they can cause problems in the eigenvalue calculations
	vector<unsigned int> bad_idx;
	for(unsigned int i=0; i<W->size2; i++){
		gsl_vector_const_view W_col = gsl_matrix_const_column(W, i);
		if(gsl_blas_dasum(&W_col.vector) < std::numeric_limits<float>::epsilon()){
			bad_idx.push_back(i);
		}
	}
	if(bad_idx.size() > 0){

		gsl_matrix_memcpy(W_tmp, W);

		gsl_permutation* permu = MatrixUtils::getPermutation(bad_idx, W->size2);
		MatrixUtils::applyPermutation(W_tmp, permu, true);
		MatrixUtils::applyPermutation(W_tmp, permu, false);

		gsl_matrix* subW = gsl_matrix_alloc(W->size1 - bad_idx.size(), W->size2 - bad_idx.size());
		gsl_matrix_const_view W_v = gsl_matrix_const_submatrix(W_tmp, 0,0, W->size1 - bad_idx.size(), W->size2 - bad_idx.size());
		gsl_matrix_memcpy(subW, &W_v.matrix);
		gsl_matrix_free(W_tmp);
		W_tmp = subW;
	} else {
		gsl_matrix_memcpy(W_tmp, W);
	}

	// first, divide W_tmp by 2
	gsl_matrix_scale(W_tmp, 0.5);

	// OK, now we have to take the eigenvalues of tmp_ss
	gsl_vector* eval = gsl_vector_alloc(W_tmp->size1);
	gsl_eigen_symm_workspace* eigen_w = gsl_eigen_symm_alloc(W_tmp->size1);
	gsl_eigen_symm(W_tmp, eval, eigen_w);
	gsl_eigen_symm_free(eigen_w);
	// Now, sort the eigenvalues in descending order
	std::sort(eval->data, eval->data + eval->size, std::greater<double>());

	int n_eval = 0;
	while(gsl_vector_get(eval, n_eval) > std::numeric_limits<float>::epsilon()
			&& static_cast<unsigned int>(++n_eval) < eval->size);

	// try the SVD instead
	//gsl_matrix* V_tmp = gsl_matrix_alloc(W_tmp->size2, W_tmp->size2);
	//gsl_vector* eval_sq = gsl_vector_alloc(W_tmp->size2);
	//gsl_vector* svd_ws = gsl_vector_alloc(W_tmp->size2);

	//gsl_linalg_SV_decomp (W_tmp, V_tmp, eval_sq, svd_ws);

	//gsl_matrix_free(V_tmp);
	//gsl_vector_free(svd_ws);

	// TODO: find where these eigenvalues go to zero and only take those
	//while(static_cast<unsigned int>(n_eval) < eval_sq->size && gsl_vector_get(eval_sq, n_eval) > std::numeric_limits<float>::epsilon()){
	//	gsl_vector_set(eval_sq, n_eval, sqrt(gsl_vector_get(eval_sq, n_eval)));
	//	++n_eval;
	//}

	gsl_matrix_free(W_tmp);


	// I don't feel like doing memory management, so use a vector instead of an array
	std::vector<double> nct(n_eval, 0);
	std::vector<int> df(n_eval, 1);
	std::vector<double> qfc_detail(7);
	int qfc_err;
	double pval;
	int lim=10000;
	double acc=std::numeric_limits<float>::epsilon();
	double sigma=0;


	// qfc is NOT thread-safe (or even reentrant!)
	_qfc_lock.lock();
	qfc(eval->data, &nct[0], &df[0], &n_eval, &sigma, &Q, &lim, &acc, &qfc_detail[0], &qfc_err, &pval);
	_qfc_lock.unlock();
	//qfc(eval_sq->data, &nct[0], &df[0], &n_eval, &sigma, &Q, &lim, &acc, &qfc_detail[0], &qfc_err, &pval);


	// clean up, please!
	gsl_vector_free(eval);
	//gsl_vector_free(eval_sq);

	return 1-pval;
}

}
}
