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

#include "qfc.h"
#include "MatrixUtils.h"

using std::vector;
using std::string;
using std::pair;

namespace BioBin{
namespace Test{

unsigned int SKATUtils::getGenoWeights(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno,
			const Bin& bin, const vector<pair<string, unsigned int> >& name_pos,
			gsl_matrix* &geno){

	unsigned int n_col = bin.getVariantSize();
	unsigned int n_row = name_pos.size();

	gsl_matrix* geno_tmp = gsl_matrix_alloc(n_row, n_col);
	gsl_vector* wt_tmp = gsl_vector_alloc(n_col);

	// get the average genotype (respecting the encoding)
	vector<float> avg_geno;
	avg_geno.reserve(bin.getVariantSize());
	Bin::const_locus_iterator ci=bin.variantBegin();
	unsigned int idx=0;
	while(ci != bin.variantEnd()){
		avg_geno.push_back(pop_mgr.getAvgGenotype(**ci));
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
		for(unsigned int i=0; i<n_col; i++){
			unsigned char g = pop_mgr.getIndivGeno(**ci,name_pos[i].second);
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
	float missing_thresh = 0.05 * n_row;
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

	// We have no SNPS!  bail out!
	if(bad_idx.size() == n_col){
		// clean up and return 1
		gsl_matrix_free(geno_tmp);
		gsl_vector_free(wt_tmp);
		geno = 0;
		return 0;
	} else if(bad_idx.size() > 0){
		// get rid of the "bad" indexes

		gsl_matrix* P = gsl_matrix_alloc(n_col, n_col);
		MatrixUtils::getPermuMatrix(bad_idx, P);

		// OK, now move the "bad" rows to the end, please!
		gsl_matrix* G_P = gsl_matrix_alloc(n_row, n_col);
		// permute the columns of geno into G_P
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, geno, P, 0.0, G_P);
		gsl_matrix_const_view G_v = gsl_matrix_const_submatrix(G_P, 0, 0,
			G_P->size1, G_P->size2 - bad_idx.size());

		// I'm officially done with "geno" now
		gsl_matrix_free(geno_tmp);
		geno = gsl_matrix_alloc(n_row, n_col = bad_idx.size());
		gsl_matrix_memcpy(geno, &G_v.matrix);
		gsl_matrix_free(G_P);

		// and do the same with the weights
		gsl_vector* wt_P = gsl_vector_alloc(n_col);
		gsl_blas_dgemv(CblasNoTrans, 1.0, P, wt_tmp, 0.0, wt_P);
		std::swap(wt_P, wt_tmp);
		gsl_vector_free(wt_P);

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
	// OK, now we have to take the eigenvalues of tmp_ss
	gsl_vector* eval = gsl_vector_alloc(W->size1);
	gsl_eigen_symm_workspace* eigen_w = gsl_eigen_symm_alloc(W->size1);
	gsl_matrix* W_tmp = gsl_matrix_alloc(W->size1, W->size2);
	gsl_matrix_memcpy(W_tmp, W);
	gsl_eigen_symm(W_tmp, eval, eigen_w);
	gsl_eigen_symm_free(eigen_w);
	gsl_matrix_free(W_tmp);

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
