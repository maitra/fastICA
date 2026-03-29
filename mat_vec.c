/* mat_vec.c
 * 
 * Basic vector and matrix functions
 * 
 * Authors: Rouben Rostamian<rostamian@umbc.edu> and Ranjan Maitra<maitra@iastate.edu>
 * Fall 1996
 * Revised January 1999
 * Revised October 2000
 * Thoroughly revised March 2005
*/

#include "mat_vec.h"

/* Allocates memory and creates an mxn Hilbert matrix of double */
double **dhilbert(SIZE_T m, SIZE_T n)
{
	SIZE_T i, j;
	double **h;

	MAKE_MATRIX(h, m, n);
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			h[i][j] = 1.0 / ( 1.0 + i + j);

	return h;
}

/* Multiplies matrices a and b and puts the result in c which should not be
 * pre-allocated.   exit() will be called if a and b are incompatible
*/
double **multiply(double **a, SIZE_T arows, SIZE_T acols,
		  double **b, SIZE_T brows, SIZE_T bcols)
{
	SIZE_T i, j, k;
	double **c;

	assert(acols == brows);
  
	MAKE_MATRIX(c, arows, bcols);
	
	for (i = 0; i < arows; i++)
		for (j = 0; j < bcols; j++) {
			c[i][j] = 0;
			for (k = 0; k < acols; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	return  c;
}

/* same as in multiply above, except that long doubles are used */

long double **longmultiply(long double **a, SIZE_T arows, SIZE_T acols,
		       long double **b, SIZE_T brows, SIZE_T bcols)
{
	SIZE_T i, j, k;
	long double **c;

	assert(acols == brows);
  
	MAKE_MATRIX(c, arows, bcols);
	
	for (i = 0; i < arows; i++)
		for (j = 0; j < bcols; j++) {
			c[i][j] = 0;
			for (k = 0; k < acols; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	return  c;
}



/* Multiplies matrices a and b and puts the result in c which should be
 * pre-allocated.   exit() will be called if a and b are incompatible
*/
void matxmat(double **a, SIZE_T arows, SIZE_T acols,
	     double **b, SIZE_T brows, SIZE_T bcols, double **c)
{
	SIZE_T i, j, k;

	assert(acols == brows);
  
	for (i = 0; i < arows; i++)
		for (j = 0; j < bcols; j++) {
			c[i][j] = 0;
			for (k = 0; k < acols; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
}

/**
 * Multiply A^t * B.
 *
 * @param a first matrix (n x p)
 * @param b second matrix (n x q)
 * @param n number of rows of two matrices
 * @param p number of columns of a
 * @param q number of columns of b
 * @return resulting matrix
 */
double **aprimeb(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q)
{
	SIZE_T i, j, k;
	double **c;

	MAKE_MATRIX(c, p, q);

	if (!c)
		return NULL;
	
	for (i = 0; i < p; i++)
		for (j = 0; j < q; j++) {
			c[i][j] = 0;
			for (k = 0; k < n; k++)
				c[i][j] += a[k][i] * b[k][j];
		}
	
	return c;
} /* aprimeb */

/**
 * Multiply A^t * B.
 *
 * @param a first matrix (n x p)
 * @param b second matrix (n x q)
 * @param n number of rows of two matrices
 * @param p number of columns of a
 * @param q number of columns of b
 * @param c contains result in p x q matrix (pre-assigned)
 */
void AprimeB(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q, double **c)
{
	SIZE_T i, j, k;

	for (i = 0; i < p; i++)
		for (j = 0; j < q; j++) {
			c[i][j] = 0;
			for (k = 0; k < n; k++)
				c[i][j] += a[k][i] * b[k][j];
		}
} /* AprimeB */

/**
 * Multiple A^t * B^t.
 *
 * @param a first matrix (n x p)
 * @param b second matrix (q x n)
 * @param n number of rows of two matrices
 * @param p number of columns of a
 * @param q number of columns of b
 * @return resulting matrix
 */
double **aprimebprime(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q)
{
	SIZE_T i, j, k;
	double **c;

	MAKE_MATRIX(c, p, q);

	if (!c)
		return NULL;
	
	for (i = 0; i < p; i++)
		for (j = 0; j < q; j++) {
			c[i][j] = 0;
			for (k = 0; k < n; k++)
				c[i][j] += a[k][i] * b[j][k];
		}
	
	return c;
} /* aprimebprime */

/**
 * Multiple C = A^t * B^t and return C^t.
 *
 * @param a first matrix (n x p)
 * @param b second matrix (q x n)
 * @param n number of rows of two matrices
 * @param p number of columns of a
 * @param q number of columns of b
 * @return resulting matrix
 */
double **aprimebprimeprime(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q)
{
	SIZE_T i, j, k;
	double **c;

	MAKE_MATRIX(c, q, p);

	if (!c)
		return NULL;
	
	for (i = 0; i < p; i++)
		for (j = 0; j < q; j++) {
			c[j][i] = 0;
			for (k = 0; k < n; k++)
				c[j][i] += a[k][i] * b[j][k];
		}
	
	return c;
} /* aprimebprimeprime */

/**
 * Multiple C = A * B^t and return C^t.
 *
 * @param a first matrix (p x n)
 * @param b second matrix (q x n)
 * @param n number of rows of two matrices
 * @param p number of columns of a
 * @param q number of columns of b
 * @return resulting matrix (q x p)
 */
double **abprimeprime(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q)
{
	SIZE_T i, j, k;
	double **c;

	MAKE_MATRIX(c, q, p);

	if (!c)
		return NULL;
	
	for (i = 0; i < p; i++)
		for (j = 0; j < q; j++) {
			c[j][i] = 0;
			for (k = 0; k < n; k++)
				c[j][i] += a[i][k] * b[j][k];
		}
	
	return c;
} /* abprimeprime */

/* Multiplies matrix a and vector x and puts the result in y which should not be
 * pre-allocated.   exit() will be called if a and x are incompatible
*/
double *matxvec(double **a, SIZE_T arows, SIZE_T acols,
		double *x, SIZE_T xrows)
{
	SIZE_T i, k;
	double *y;

	assert(acols == xrows);
	
	MAKE_VECTOR(y, arows);
	
	for (i = 0; i < arows; i++) {
		y[i] = 0;
		for (k = 0; k < acols; k++)
			y[i] += a[i][k] * x[k];
	}
	return y;
}

/* same as matxvec above, but using long doubles */

long double *longmatxvec(long double **a, SIZE_T arows, SIZE_T acols,
		long double *x, SIZE_T xrows)
{
	SIZE_T i, k;
	long double *y;

	assert(acols == xrows);
	
	MAKE_VECTOR(y, arows);
	
	for (i = 0; i < arows; i++) {
		y[i] = 0;
		for (k = 0; k < acols; k++)
			y[i] += a[i][k] * x[k];
	}
	return y;
}



/* Multiplies matrix a and vector x and puts the result in y which should be
 * pre-allocated.   exit() will be called if a and x are incompatible
*/
void matXvec(double **a, SIZE_T arows, SIZE_T acols, double *x, SIZE_T xrows,
	double *y)
{
	SIZE_T i, k;

	assert(acols == xrows);
	
	for (i = 0; i < arows; i++) {
		y[i] = 0;
		for (k = 0; k < acols; k++)
			y[i] += a[i][k] * x[k];
	}
}



/* PrSIZE_Ts matrix with a spefified format */
void print_dmatrix(double **a, SIZE_T rows, SIZE_T cols, const char *format)
{
	SIZE_T i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++)
			printf(format, a[i][j]);
		putchar('\n');
	}
	printf("\n");
}

/* PrSIZE_Ts vector with a spefified format */
void print_dvector(double *a, SIZE_T rows, const char *format)
{
	SIZE_T i;
	for (i=0; i<rows; i++)
		printf(format, a[i]);
	printf("\n");
}

void matrpose(double **a, SIZE_T rows, SIZE_T cols,double **aT)
{
	SIZE_T i,j;
	for (i=0;i<rows;i++)
		for (j=0;j<cols;j++)
			aT[j][i]=a[i][j];
}

/* Calculates and returns the Euclidean norm of a vector.
 * Modeled after BLAS's norm2 routine.
 *
 * The `sum' variable accummulates (x_i/scale)^2, where `scale'
 * is the largest |x_i| seen so far.
*/
double dEnorm(double *x, SIZE_T n)
{
	double scale = 0.0;
	double sum = 1.0;
	double b, r;
	SIZE_T i;

	if (n < 1)
		return 0.0;
	else if (n == 1)
		return fabs(x[0]);

	for (i=0; i<n; i++) {

		if (x[i]==0.0)
			continue;

		b = fabs(x[i]);

		if (b < scale) {
			r = b/scale;
			sum += r*r;
		} else {
			r = scale/b;
			sum = sum*r*r + 1.0;
			scale = b;
		}
	}

	return scale * sqrt(sum);
}

/*copy matrix A to matrix B*/
void cpy(double **a, SIZE_T nrows, SIZE_T ncols,double **b)
{
	SIZE_T i,j;
	for(i=0;i<nrows;i++)
		for (j=0;j<ncols;j++)
			b[i][j]=a[i][j];
}

double quadratic(double **A, double *x, SIZE_T p) 
{                  /*calculates x'Ax where A is a p x p-dimensional matrix */
	SIZE_T i,j;
	double qform=0.0;
	for(i=0;i<p;i++)
		for(j=0;j<p;j++)
			qform+=x[i]*x[j]*A[i][j];
	return(qform);
}

double ltquadratic(double *ltA, double *x, SIZE_T p) 
{                  /*calculates x'Ax where A is a p x p-dimensional symmetric 
		     matrix in packed lower-triangular form*/
	SIZE_T i,j;
	double qform=0.0;
	for(i=0;i<p;i++) {
		qform+=x[i]*x[i]*ltA[((i+1)*(i+2))/2-1];
		for(j=0;j<i;j++)
			qform+=2.*x[i]*x[j]*ltA[i*(i+1)/2+j];
	}
	return(qform);
}


SIZE_T ar(double **A, SIZE_T n, double rho) /* sets up a AR-1 correlation matrix 
				       with rho */
{
	SIZE_T i,j;
	double *x;
	MAKE_VECTOR(x, n);
	x[0]=1.0;
	for(i=1;i<n;i++)
		x[i]=x[i-1]*rho;
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			A[i][j]=x[abs(i-j)];
	FREE_VECTOR(x);
	return 0;
}

SIZE_T arinv(double **A, SIZE_T n,double rho) /* sets up the inverse of a AR-1
					  correlation matrix  with rho */
{
  SIZE_T i,j;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      if (((i==0) && (j==0)) || ((i==(n-1)) && (j==(n-1)))) {
	A[i][i]=1./(1-rho*rho);
      }
      else {
	if (i==j) {
	  A[i][i]=(1+rho*rho)/(1-rho*rho);
	}
	else {
	  if (abs(i-j)==1) {
	    A[i][j]=-rho/(1-rho*rho);
	  }
	  else {
	    A[i][j]=0.;
	  }
	}
      }
    }
  }
  return 0;
}

/* Function to perform the Cholesky decomposition of the n x n symmetric 
   matrix A. The matrix A is in lower-triangular form as also the output.  
   The code is written following the details in:

   http://www.math-linux.com/spip.php?article43

   though I also verified it by deriving it. A fuller note will follow and 
   perhaps be included in the class notes.

   Note that Wikipedia calls this the Cholesky-Crout algorithm.

   The original Cholesky decomposition as in Golub and van Loan is a bit more 
   computer-SIZE_Tensive and more stable. 

   Copyright Ranjan Maitra <maitra@iastate.edu>, Ames, IA, USA 2012/12/23.
*/

/* The Cholesky-Crout algorithm returns 1 if the input matrix A is not
   positive-definite, The lower-triangle of the Cholesky decomposition 
   is returned in L **only** when the return value is zero. */   

SIZE_T cholesky(double *A, SIZE_T n, double *L)
{
	double tmp;
	SIZE_T i, j, k, l, ind = 0;

	/* now for the jth column where j runs from 1 through n */
	
	for (j = 0; (j < n && !ind); j++) {
		/* first obtain the jth diagonal */
		double sum = 0;
		k = LTINDEX(j, j);
		for (l = 0; l < j; l++) 
			sum += SQ(L[LTINDEX(j, l)]);
		tmp = A[k] - sum;

		if (tmp > 0) {
			L[k] = sqrt(tmp);
			/* now get the other elements in lower triangle of the
			   jth column */
			for (i = (j + 1); i < n; i++) {
				double sum = 0;
				k = LTINDEX(i, j);
				for (l = 0; l < j ; l++) 
					sum += L[LTINDEX(i, l)] * L[LTINDEX(j, l)];
				L[k] = (A[k] - sum) / L[LTINDEX(j, j)];
			}
		}
		else 
			ind = 1; /* not a positive-definite matrix */
	}
	return ind;
}

/* 
   A faster version of Cholesky decomposes a positive definite matrix A SIZE_To
   A = LDL' where L is lower lower triangle (lower triangle without the 
   diagonal elements which are all unity and hence do not need to be stored) 
   and the diagonal matrix D. It is faster because it does not involve taking
   a square root operation.
   
   Note that here, L is a vector of n * (n - 1)/2 elements, and D is a vector
   of n elements.

   Runtime experiments indicate that this approach is substantially faster in
   terms of system time, but marginally slower in terms of user time. Since
   the bulk of CPU time is the user time, it is slower overall than the other
   method. So we would perhaps prefer the other method (cholesky).

   Copyright Ranjan Maitra <maitra@iastate.edu>, Ames, IA, USA 2012/12/23.

 */
		


SIZE_T faster_cholesky(double *A, SIZE_T n, double *L, double *D)
{
	SIZE_T i, j, k, l, ind = 0;

	/* for the jth column where j runs from 1 through n */
	
	for (j = 0; (j < n && !ind); j++) {
		/* first obtain the jth diagonal */
		double sum = 0;
		k = LTINDEX(j, j);
		for (l = 0; l < j; l++) 
			sum += SQ(L[LLTINDEX(j, l)]) * D[l];
		D[j] = A[k] - sum;

		if (D[j] > 0) {
			/* now get the other elements in lower triangle of the
			   jth column */
			for (i = (j + 1); i < n; i++) {
				double sum = 0;
				for (l = 0; l < j ; l++) 
					sum += L[LLTINDEX(i, l)] * L[LLTINDEX(j, l)] * D[l];

				L[LLTINDEX(i, j)] = (A[LTINDEX(i, j)] - sum) / D[j];
			}
		}
		else 
			ind = 1; /* not a positive-definite matrix */
	}
	return ind;
}

		


/* This function calculates the determinant of a lower triangular matrix 
   without using any external library. It is assumed that the matrix is 
   positive definite whereby the function returns 0, otherwise the function
   returns a non-zero value and the logdet value is to be ignored.

   Copyright Ranjan Maitra <maitra@iastate.edu> Ames, IA 50011-1210, USA.
*/

SIZE_T logLTdeterminant(double *LTA, SIZE_T p, double *logdet)
{
	double *L;
	SIZE_T i, k, j = LTINDEX(p, p);
	
	*logdet = 0;

	MAKE_VECTOR(L, j);
	k = cholesky(LTA, p, L);
	if (!k) {
		for (i = 0; i < p; i++) 
			*logdet += log(L[LTINDEX(i, i)]);
	}
	FREE_VECTOR(L);
	return k;
}

double L2norm(SIZE_T n, double *y)
{
	SIZE_T i;
	double sum = 0;

	for (i = 0; i < n; i++) 
		sum += SQ(y[i]);
	
	return sqrt(sum);
}

/* same as L2norm but uses double precision arithmetic */

long double longL2norm(SIZE_T n, long double *y)
{
	SIZE_T i;
	long double sum = 0;

	for (i = 0; i < n; i++) 
		sum += SQ(y[i]);
	
	return sqrt(sum);
}






/* Multiplies matrix a and vector x and puts the result in y which should be
 * pre-allocated.   exit() will be called if a and x are incompatible
*/
void ltmatxvec(double *a, SIZE_T p, double *x, double *ax)
{
	SIZE_T i, k;
	
	for (i = 0; i < p; i++) {
		ax[i] = 0;
		for (k = 0; k < p; k++)
			ax[i] += a[LTINDEX(i, k)] * x[k];
	}
}
/* Multiplies matrix a and vector x and puts the result in y which should be
 * pre-allocated.   exit() will be called if a and x are incompatible
 * same as ltmatvec, but uses double precision.
*/
void longltmatxvec(long double *a, SIZE_T p, long double *x, long double *ax)
{
	SIZE_T i, k;
	
	for (i = 0; i < p; i++) {
		ax[i] = 0;
		for (k = 0; k < p; k++)
			ax[i] += a[LTINDEX(i, k)] * x[k];
	}
}





double vecxvec(double *x, SIZE_T p, double *y)
{
	double sum = 0;
	SIZE_T i;
	for (i = 0; i < p; i++)
		sum += x[i] * y[i];
	return sum;
}

/* same as vecxvec, but uses double precision */

long double longvecxvec(long double *x, SIZE_T p, long double *y)
{
	long double sum = 0;
	SIZE_T i;
	for (i = 0; i < p; i++)
		sum += x[i] * y[i];
	return sum;
}


/**
 * Multiple $X^tX$ for a matrix $X$.
 *
 * @param x the matrix
 * @param n the number of rows in the matrix
 * @param p the number of columns in the matirx
 * @return a lower triangular matrix with result
 */
double *XprimeX(double **x, SIZE_T n, SIZE_T p)
{
	SIZE_T i, j, l;		/* indices */
	double *ltxtx;		/* lower triangle matrix to be filled */
	SIZE_T dim = p*(p+1)/2;	/* expensive to recompute in for loop */

	MAKE_VECTOR(ltxtx, dim);

	if (!ltxtx)
		return NULL;

	for (i = 0; i < dim; i++)
		ltxtx[i] = 0.;

	for (i = 0; i < n; i++)
		for (j = 0; j < p; j++) {
			dim = j*(j+1)/2;
			for (l = 0; l <= j; l++) 
				ltxtx[dim + l] +=  x[i][j] * x[i][l];
		}

	return ltxtx;
} /* XprimeX */


/**
 * Compute just a portion of $X^tX$ for a matrix $X$.  This function is
 * identical to XprimeX() except it computes the product for a subset of rows
 * in $X$ as identified by argument #index.
 *
 * @param x the matrix
 * @param n the number of rows in the matrix
 * @param p the number of columns in the matirx
 * @param index identifies n rows of a larger matrix x to use in calculations
 * @return a lower triangular matrix with result
 */
double *indexed_XprimeX(double **x, SIZE_T n, SIZE_T p, SIZE_T *index)
{
	SIZE_T i, j, l;		/* indices */
	double *ltxtx;		/* lower triangle matrix to be filled */
	SIZE_T dim = p*(p+1)/2;	/* expensive to recompute in for loop */

	MAKE_VECTOR(ltxtx, dim);

	if (!ltxtx)
		return NULL;

	for (i = 0; i < dim; i++)
		ltxtx[i] = 0.;

	for (i = 0; i < n; i++)
		for (j = 0; j < p; j++) {
			dim = j*(j+1)/2;
			for (l = 0; l <= j; l++) 
				ltxtx[dim + l] +=  x[index[i]][j] * x[index[i]][l];
		}

	return ltxtx;
} /* indexed_XprimeX */

double *Linverse(SIZE_T p, double *L)
{
	/* function to calculate and return the inverse of a p x p lower 
	   triangular matrix stored in packed lower-triangular format. 
	   returns the inverse of L in lower packed triangular format.
	*/
	SIZE_T i, j, k;
	double *Linv;
	
	MAKE_VECTOR(Linv, p*(p+1)/2);
	for (i = 0; i < p*(p+1)/2; i++)
		Linv[i] = L[i];

	for (i = 0; i < p; i++) {
		Linv[LTINDEX(i, i)] = 1/L[LTINDEX(i, i)];
		for (j = i + 1; j < p; j++) {
			double sum=0.0;
			for (k = i; k < j; k++) sum -= Linv[LTINDEX(j, k)] * Linv[LTINDEX(k, i)];
			Linv[LTINDEX(j, i)] = sum / L[LTINDEX(j, j)];
		}
	}
	return Linv;
}

double *LL_multiply(SIZE_T p, double *L, double *J)
{
	/*function to calculate the product of two lower triangular matrices 
	  (input and return output in packed lower-triangular form) */
	SIZE_T i, j, k;
	double *M;

	MAKE_VECTOR(M, p*(p+1)/2);

	for (i = 0; i < p*(p+1)/2; i++) 
		M[i] = 0;
	for (i = 0; i < p; i++)
		for (j = 0; j <= i; j++)
			for (k = j; k <= i; k++)
				M[LTINDEX(i, j)] += L[LTINDEX(i,k)] * J[LTINDEX(k,j)];
	return M;
}


double *LLT(double *l, SIZE_T p)
{
	/* function to calculate LL' input and output both in lower-triangular 
	   form 
	*/
	SIZE_T i, j, k = p * (p + 1)/2, n;
	double *m;

	MAKE_VECTOR(m, k);
	for (i = 0; i < k; i++)
		m[i] = 0;
	for (i = 0; i < p; i++) {
		for (j = 0; j <= i; j++) {
			for (n = 0; n <= j; n++)
				m[LTINDEX(i, j)] += l[LTINDEX(i, n)] * l[LTINDEX(j, n)];
		}
	}
	return m;
}

/* Given a symmetric positive definite matrix in packed lower triangular form
   this function computes the inverse of the matrix using a Cholesky 
   decomposition. If the input matrix is not positive definite, the return is
   nonnull (1) and the output is meaningless. */

SIZE_T cholesky_inverse(double *A, SIZE_T n, double *Ainv)
{
	double *L;
	MAKE_VECTOR(L, n * (n + 1)/2);

	if (!cholesky(A, n, L)) {
		double *Linv, *ainv;
		SIZE_T i;

		Linv = Linverse(n, L);
		ainv = LLT(Linv, n);
		FREE_VECTOR(Linv);

		for (i = 0; i < n*(n + 1)/2; i++)
			Ainv[i] = ainv[i];
		FREE_VECTOR(ainv);
		FREE_VECTOR(L);
		return 0;
	}
	else {
		FREE_VECTOR(L);
		return 1;
	}
}


