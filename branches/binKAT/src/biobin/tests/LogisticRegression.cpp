/*
 * LogisticRegression.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: jrw32
 */

#include "LogisticRegression.h"

#include <iostream>
#include <limits>
#include <cmath>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

#include "detail/MatrixUtils.h"

using std::numeric_limits;
using std::string;

using boost::array;

using BioBin::Test::Regression;

namespace BioBin {

namespace Test {

string LogisticRegression::testname = LogisticRegression::doRegister("logistic");

LogisticRegression::~LogisticRegression() {
	// TODO Auto-generated destructor stub
}

void LogisticRegression::init(){
	if(_pheno_ptr->getStatus().first.count() == 0 ||
			_pheno_ptr->getStatus().second.count() == 0){
		std::cerr << "WARNING: Case or Control completely missing; "
				<< "Logistic regression will fail" << std::endl;
		_willfail = true;
	} else {
		regressionSetup(*_pop_mgr_ptr, *_pheno_ptr);
		if(!_null_result->_conv){
			std::cerr << "WARNING: Logistic Regression null model diverged; "
							<< "Logistic regression will likely fail" << std::endl;
			_willfail = true;
		}
	}
}

double LogisticRegression::runTest(const Bin& bin) const {
	// If we are guaranteed to fail, bail out early
	if (_willfail) {
		return 1;
	}

	for(unsigned int i=0; i<_samp_name.size(); i++){
		gsl_matrix_set(_data, i, _data->size2 - 1,
				_pop_mgr_ptr->getTotalIndivContrib(bin, _samp_name[i].second, *_pheno_ptr));
	}

	// run the model now
	Regression::Result* r = calculate(*_phenos, *_data);

	// Get the p-value of the last term
	double se = sqrt(gsl_matrix_get(r->cov, r->cov->size1 - 1, r->cov->size2 - 1));
	double c = gsl_vector_get(r->beta, r->beta->size - 1);
	double pval = 2*gsl_cdf_tdist_Q(fabs(c / se),_data->size1 - _data->size2 + 1);

	delete r;

	return pval;

}

float LogisticRegression::getPhenotype(const PopulationManager& pop_mgr,
		const Utility::Phenotype& pheno, const std::string& samp) const{

	float val = std::numeric_limits<float>::quiet_NaN();

	// get 0 or 1 based on the phenotypic status
	unsigned int pos = pop_mgr.getSamplePosition(samp);
	if(pos <= pheno.getStatus().first.size()){
		if(pheno.getStatus().first[pos]){
			val = 0;
		} else if(pheno.getStatus().second[pos]){
			val = 1;
		}
	}

	return val;

}

Regression::Result* LogisticRegression::calculate(const gsl_vector& Y, const gsl_matrix& X) const {

	static const unsigned int maxIterations = 30;

	// val is the value of the logit function
	// deriv is the derivative of the logit
	// log_val and log_val_c are log(val) and log(1-val), respectively.
	// NOTE: the reason we do them here is to prevent underflow; for large values
	// of the exponent, we need to use approximations for log_val and log_val_c
	//double val, deriv, log_val, log_val_c;

	// This is the current estimate of the parameters
	// Note: position 0 is reserved for the intercept

	gsl_vector* beta = gsl_vector_calloc(X.size2);
	gsl_matrix* r_cov = gsl_matrix_calloc(X.size2, X.size2);

	Regression::Result* r = new Regression::Result(beta, r_cov);

	// Add up all the values in Y
	// note that since Y holds 0 or 1, we can just sum the absolute values here
	double sum_Y = gsl_blas_dasum(&Y);

	// get the exponent (and derivative, log and 1-log values) for the
	// null model (i.e., best fit of the intercept parameter)
	double beta_0 = log(sum_Y / (Y.size - sum_Y));
	array<double, 4> null_v = linkFunction(beta_0);

	// check for an intercept-only model here
	if(X.size2 == 1){
		gsl_vector_const_view X_col1 = gsl_matrix_const_column(&X, 0);
		if(gsl_blas_dasum(&X_col1.vector) == X.size1){
			// if we're here, it's an intercept-only model - trivial to get
			gsl_vector_set(beta, 0, beta_0);
			gsl_matrix_set(r_cov, 0, 0, null_v[1]);
			r->chisq = null_v[1]; // I don't know what to put here - sounds good
			r->_conv = true;
			r->resid = gsl_vector_alloc(Y.size);
			double pred_val = - null_v[0];
			// set the residual at first to -\mu, where \mu == 1/(1+exp(-\beta_0))
			gsl_vector_set_all(r->resid, pred_val);
			// and set resid = Y + resid (i.e., resid = Y + (-\mu)
			gsl_blas_daxpy(1, &Y, r->resid);

			return r;
		}
	}

	gsl_vector* weight = gsl_vector_calloc(Y.size);
	gsl_matrix* P = gsl_matrix_alloc(X.size2,X.size2);
	// First, let's check for colinearity!
	unsigned int n_drop = MatrixUtils::checkColinear(&X, P);
	unsigned int n_indep = X.size2 - n_drop;

	// Right-hand side of the IRLS equation.  Defined to be X*w_t + S_t^-1*(y-mu_t)
	// Or, in our parlance: rhs_i = (X*beta_t)_i + 1/deriv * (y_i - val)
	gsl_vector* rhs = gsl_vector_alloc(Y.size);

	// I need these to work with the default values
	gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(Y.size, n_indep);
	//gsl_matrix_set_all(cov_mat, -1);
	gsl_matrix_view cov_view = gsl_matrix_submatrix(r->cov, 0, 0, n_indep, n_indep);
	gsl_matrix* cov = &cov_view.matrix;
	gsl_matrix* A = gsl_matrix_calloc(Y.size, X.size2);
	double tmp_chisq;

	// Let's perform our permutation ans set A = data * P
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &X, P, 0.0, A);

	// this is the previous beta vector, for checking convergence
	gsl_vector* b_prev = gsl_vector_alloc(n_indep);
	// make this b_prev nonzero to begin
	gsl_vector_set_all(b_prev, 1.0);
	gsl_vector_view b = gsl_vector_subvector(r->beta, 0, n_indep);
	gsl_matrix_const_view X_v = gsl_matrix_const_submatrix(A, 0, 0, Y.size, n_indep);

	double LLp = std::numeric_limits<double>::infinity(); // stores previous value of LL to check for convergence
	double LLn, LL = 0;

	LLn = 0;

	for(unsigned int i=0; i<Y.size; i++){
		LLn -= 2 *(gsl_vector_get(&Y, i) * null_v[2] + (1-gsl_vector_get(&Y, i)) * null_v[3]);
	}

	unsigned int numIterations = 0;

	// set the tolerances in single precision, but do work in double precision!
	double TOL= numeric_limits<float>::epsilon();
	double MAX_vec = numeric_limits<float>::max();
	// 2*numeric_limits<double>::epsilon() ??

	// check for the following:
	// 1) convergence of the likelihood
	// 1a) convergence of the beta vector (distance stored in b_prev)
	// 2) divergence of one (or more) coefficients
	// 3) maximum number of iterations
	while ((fabs(LLp - LL) > TOL*LLn || gsl_blas_dasum(b_prev) > TOL*n_indep*gsl_blas_dasum(&b.vector)) &&
		   gsl_blas_dasum(&b.vector) < MAX_vec &&
		   ++numIterations < maxIterations ) {

		// save the old beta vector in b_prev
		gsl_vector_memcpy(b_prev, &b.vector);

		// First, let's initialize the RHS to X*beta_t (rhs = 1 * X * b + 0* rhs)
		gsl_blas_dgemv(CblasNoTrans, 1, &X_v.matrix, &b.vector, 0, rhs);

		LLp = LL;
		LL = 0;

		// add to LL for each row
		for (unsigned int i = 0; i < Y.size; i++) {

			// calculate the value of the exponent for the individual
			double v = gsl_vector_get(rhs, i);

			// At this point, v is the value of the exponent

			array<double, 4> v_arr = linkFunction(v);

			// calculate LL for this ind and add to running total
			LL -= 2 *(gsl_vector_get(&Y, i) * v_arr[2] + (1-gsl_vector_get(&Y, i)) * v_arr[3]);

			// get the weight and update the rhs for IRLS
			//weight[i] = v_arr[1];
			gsl_vector_set(weight, i, v_arr[1] );
			gsl_vector_set(rhs, i, v + 1/v_arr[1] * (gsl_vector_get(&Y, i) - v_arr[0]));

		}

		// Look, magic!
		gsl_multifit_wlinear(&X_v.matrix, weight, rhs, &b.vector, cov, &tmp_chisq, ws);

		// check for NaNs here
		if(std::isfinite(gsl_blas_dasum(&b.vector))){
			// get the difference between the old beta and the new beta
			gsl_vector_sub(b_prev, &b.vector);
		} else {
			// terminate the iteration, giving us the previous beta
			gsl_vector_memcpy(&b.vector, b_prev);
			// and set the "difference" to 0
			gsl_vector_set_zero(b_prev);
			// set loglikelihood to NaN
			LL = std::numeric_limits<double>::quiet_NaN();
		}

		LL += 0;

	} // complete iteration

	// set the residuals here
	r->resid = gsl_vector_alloc(X_v.matrix.size1);
	gsl_multifit_linear_residuals(&X_v.matrix, rhs, &b.vector, r->resid);

	// nonconvergence happens if:
	// -Log likelihood is not finite (inf or NaN)
	// too many iteratons
	// The current log likelihood is less than the null model
	if(!std::isfinite(LL) ||
	   numIterations >= maxIterations ||
	   LL-LLn > 0 ){
		r->_conv = false;
	}

	r->chisq = tmp_chisq;

	gsl_vector_free(weight);
	gsl_vector_free(b_prev);
	gsl_vector_free(rhs);
	gsl_multifit_linear_free(ws);
	gsl_matrix_free(A);

	// OK, now time to unpermute everything!
	// Note: to unpermute, multiply by P transpose!
	// Also, we need to unpermute both the rows AND columns of cov_mat
	//b = gsl_vector_view_array(r->beta_vec, n_cols);
	gsl_matrix* _cov_work = gsl_matrix_calloc(X.size2, X.size2);
	// permute columns
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, r->cov, P, 0.0, _cov_work);
	// permute rows
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, _cov_work, 0.0, r->cov);
	gsl_matrix_free(_cov_work);

	gsl_vector* _bv_work = gsl_vector_alloc(X.size2);
	gsl_vector_memcpy(_bv_work, r->beta);
	gsl_blas_dgemv(CblasTrans, 1.0, P, _bv_work, 0.0, r->beta);

	// create a set of all of the removed indices
	// I'm going to re-use _bv_work from earlier to save a few bytes of memory
	gsl_vector_free(_bv_work);

	gsl_matrix_free(P);


	return r;
}

array<double, 4> LogisticRegression::linkFunction(double v){
	// returns val, deriv (val * 1-val), log(val), 1-log(val)

	array<double, 4> retval;
	// max_val is the maximum value of v above that will not result in
	// loss of precision. (with a factor of 2 in there for good luck)
	static const double max_val = -log(numeric_limits<double>::epsilon());

	if (v > max_val) {
		retval[0] = 1;
		retval[1] = exp(-v);
		// log(f(x)) obtained by Taylor series expansion of log(x) at x=1
		// note that f(x) - 1 = -exp(-x)/(1+exp(-x)) ~= -exp(-x)
		// Also, that shows the log(1-f(x)) ~= exp(-x)
		retval[2] = -exp(-v);
		retval[3] = -v;
	} else {
		// we won't underflow here
		retval[0] = 1 / (1 + exp(-v));
		retval[2] = log(retval[0]);
		// however, we might underflow when calculating derivatives and log
		if (-v > max_val) {
			// the traditional derivative WILL underflow
			retval[1] = exp(v);
			retval[3] = -exp(v);
		} else {
			retval[1] = retval[0] * (1 - retval[0]);
			retval[3] = log(1 - retval[0]);
		}
	}
	return retval;
}


}
}
