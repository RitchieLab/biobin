/*
 * MatrixUtils.h
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_TEST_MATRIXUTILS_H
#define BIOBIN_TEST_MATRIXUTILS_H

#include<gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include <vector>

namespace BioBin {

namespace Test {

class MatrixUtils {
public:
	~MatrixUtils() {}
private:
	MatrixUtils() {}

public:

	/*
	 * \brief A function to check the colinearity of the columns of a given matrix
	 * This function will check if the columns in the given matrix are colinear
	 * or not.  Also, as a secondary parameter, it will return a permutation
	 * matrix that will move all linearly dependent columns to the end of the
	 * matrix.
	 *
	 * NOTE: copied from PLATO!
	 *
	 * Note that X must have at least as many rows as columns, and P should
	 * be a square matrix with same number of columns as X
	 *
	 * In order to get the indices of the colinear columns, you can take the vector
	 * v = [0 1 2 3 ...] and calculate a = Pv.  Then the last k entries of a are
	 * the indices of the colinear columns (here k is the return value of this
	 * function)
	 *
	 * \param X the matrix to check
	 * \param P the permutation matrix to be returned
	 * \return the number of colinear columns
	 */
	static unsigned int checkColinear(const gsl_matrix*, gsl_matrix* &P);

	// gets a permutation matrix obtained by moving the columns identified
	// by idx to the end.  We guarantee that the final idx.size() columns
	// will be those identified by the indices in idx (in no particular order)
	static void getPermuMatrix(const std::vector<unsigned int>& idx, gsl_matrix* &P);

	static void applyPermutation(gsl_matrix* mat, const gsl_permutation* permu, bool byCol = true);
	static gsl_permutation* getPermutation(const std::vector<unsigned int>& idx_permu, unsigned int size);

};

}

}

#endif /* MATRIXUTILS_H_ */
