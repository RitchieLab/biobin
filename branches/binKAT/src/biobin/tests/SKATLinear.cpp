/*
 * SKATLinear.cpp
 *
 *  Created on: Feb 18, 2015
 *      Author: jrw32
 */

#include "SKATLinear.h"

#include <vector>
#include <set>

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
}

void SKATLinear::init(const PopulationManager& pop_mgr, const Utility::Phenotype& pheno){
	_base_reg.init(pop_mgr, pheno);

	pop_mgr_ptr = &pop_mgr;
	_pheno_ptr = &pheno;

	// get the residual vector from the null model
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(_base_reg._data, 0,0,
			_base_reg._data->size1, _base_reg._data->size2-1);
	_resid = gsl_vector_alloc(_base_reg._data->size1);
	gsl_multifit_linear_residuals(&X_v.matrix, _base_reg._phenos,
			_base_reg._null_result->beta, _resid);
}

double SKATLinear::runTest(const Bin& bin) const{
	// first things first, let's set up the genotype matrix
	gsl_matrix* geno = gsl_matrix_alloc(_resid.size, bin.getVariantSize());

	// get the average genotype (respecting the encoding)
	vector<float> avg_geno;
	avg_geno.reserve(bin.getVariantSize());
	Bin::const_locus_iterator ci=bin.variantBegin();
	while(ci != bin.variantEnd()){
		avg_geno.push_back(_pop_mgr_ptr->getAvgGenotype(**ci));
		++ci;
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

	gsl_matrix* weights = gsl_matrix_calloc(G_P->size2, G_P->size2);



	// now, get the Q statistic, defined to be
	// r*GKG*r

	return 1;
}


}

}
