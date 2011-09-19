/*
 * dgeev.h
 * ($Revision: $)
  This file contains a wrapper around the LAPACK dgeev routine,
  used to calculate eigenvalues and associated eigenvectors
  for a square asymmetric matrix H. Since the matrix is asymmetrix,
  both real and imaginary components of the eigenvalues are returned.
  For real, symmetric matricies, the eigenvalues will be real, and
  hence the complex terms equal to zero.

  There are two function calls defined in this header, of the
  forms

    void dgeev(double **H, int n, double *Er, double *Ei)
    void dgeev(double **H, int n, double *Er, double *Ei, double **Evecs)

    H: the n by n matrix that we are solving.
    n: the order of the square matrix H.
    Er: an n-element array to hold the real parts of the eigenvalues.
    Ei: an n-element array to hold the imaginary parts of the eigenvalues.
    Evecs: an n by n matrix to hold the eigenvectors.
*/

#include <cmath>

void dgeev(double **H, int n, double *Er, double *Ei);
void dgeev(double **H, int n, double *Er, double *Ei, double **Evecs);

double *dgeev_ctof(double **in, int rows, int cols);
void dgeev_ftoc(double *in, double **out, int rows, int cols);
void dgeev_sort(double *Er, double *Ei, int N);
void dgeev_sort(double *Er, double *Ei, double **Evecs, int N);
 
extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a,
		       int *lda, double *wr, double *wi, double *vl,
		       int *ldvl, double *vr, int *ldvr,
		       double *work, int *lwork, int *info);


void dgeev(double **H, int n, double *Er, double *Ei)
{
  char jobvl, jobvr;
  int lda,  ldvl, ldvr, lwork, info;
  double *a, *vl, *vr, *work;
  
  jobvl = 'N'; // V/N to calculate/not calculate the left eigenvectors of the matrix H.
  jobvr = 'N'; // As above, but for the right eigenvectors.

  lda = n; // The leading dimension of the matrix a.
  a = dgeev_ctof(H, n, lda); // Convert the matrix H from double pointer C form to single pointer Fortran form.

  /* Whether we want them or not, we need to define the matrices
     for the eigenvectors, and give their leading dimensions.
     We also create a vector for work space. */

  ldvl = n;
  vl = new double[n*n];
  ldvr = n;
  vr = new double[n*n];
  work = new double[4*n];
  lwork = 4*n;

  dgeev_(&jobvl, &jobvr, &n, a, &lda, Er, Ei, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  dgeev_sort(Er, Ei, n); //Sort the results by eigenvalue in decreasing magnitude.

  delete [] a;
  delete [] vl;
  delete [] vr;
  delete [] work;
}


void dgeev(double **H, int n, double *Er, double *Ei, double **Evecs)
{
  char jobvl, jobvr;
  int lda,  ldvl, ldvr, lwork, info;
  double *a, *vl, *vr, *work;
  
  jobvl = 'N';
  jobvr = 'V';
  lda = n;
  a = dgeev_ctof(H, n, lda);

  ldvl = n;
  vl = new double[n*n];
  ldvr = n;
  vr = new double[n*n];
  work = new double[4*n];
  lwork = 4*n;

  dgeev_(&jobvl, &jobvr, &n, a, &lda, Er, Ei, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  dgeev_ftoc(vr, Evecs, n, ldvr);
  dgeev_sort(Er, Ei, Evecs, n);

  delete [] a;
  delete [] vl;
  delete [] vr;
  delete [] work;
}


double* dgeev_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++)
	  for (j=0; j<cols; j++)
		  out[i+j*cols] = in[i][j];
  return(out);
}


void dgeev_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			out[i][j] = in[i+j*cols];
}


void dgeev_sort(double *Er, double *Ei, int N)
{
	double temp, *E2;
	int i, j;

	E2 = new double[N];
	for (i=0; i<N; i++)
		E2[i] = Er[i]*Er[i]+Ei[i]*Ei[i];

	for (j=0; j<N; j++)
		for (i=0; i<N-1; i++)
			if (fabs(E2[i])<fabs(E2[i+1]))
			{
				temp = E2[i]; E2[i] = E2[i+1]; E2[i+1] = temp;
				temp = Er[i]; Er[i] = Er[i+1]; Er[i+1] = temp;
				temp = Ei[i]; Ei[i] = Ei[i+1]; Ei[i+1] = temp;
			}

	delete [] E2;
}


void dgeev_sort(double *Er, double *Ei, double **Evecs, int N)
{
	double temp, *E2;
	int i, j, k;

	E2 = new double[N];
	for (i=0; i<N; i++)
		E2[i] = Er[i]*Er[i]+Ei[i]*Ei[i];

	for (j=0; j<N; j++)
		for (i=0; i<N-1; i++)
			if (fabs(E2[i])<fabs(E2[i+1]))
			{
				temp = E2[i]; E2[i] = E2[i+1]; E2[i+1] = temp;
				temp = Er[i]; Er[i] = Er[i+1]; Er[i+1] = temp;
				temp = Ei[i]; Ei[i] = Ei[i+1]; Ei[i+1] = temp;

				for (k=0; k<N; k++)
				{
					temp = Evecs[k][i];
					Evecs[k][i] = Evecs[k][i+1];
					Evecs[k][i+1] = temp;
				}
			}

	delete [] E2;
}
