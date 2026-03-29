#ifndef _MATVEC_H
#define _MATVEC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "array.h"
#include "constants.h"
#include "util.h"

SIZE_T cholesky(double *A, SIZE_T n, double *L);

SIZE_T faster_cholesky(double *A, SIZE_T n, double *L, double *D);

double **dhilbert(SIZE_T m, SIZE_T n);
double **identmatrix(SIZE_T n);
double dEnorm(double v[], SIZE_T n);
double **multiply(double **a, SIZE_T arows, SIZE_T acols,
		  double **b, SIZE_T brows, SIZE_T bcols);
long double **longmultiply(long double **a, SIZE_T arows, SIZE_T acols,
			   long double **b, SIZE_T brows, SIZE_T bcols);
void matxmat(double **a, SIZE_T arows, SIZE_T acols,
	     double **b, SIZE_T brows, SIZE_T bcols, double **c);
void print_dmatrix(double **a, SIZE_T rows, SIZE_T cols, const char *format);

void print_dvector(double *a, SIZE_T rows, const char *format);

double *matxvec(double **a, SIZE_T arows, SIZE_T acols, double *x, SIZE_T xrows);
long double *longmatxvec(long double **a, SIZE_T arows, SIZE_T acols, long double *x, SIZE_T xrows);
void matXvec(double **a, SIZE_T arows, SIZE_T acols, double *x, SIZE_T xrows, double *y);

void matrpose(double **a,SIZE_T rows,SIZE_T cols,double **aT);
SIZE_T arinv(double **A,SIZE_T n,double rho); /* sets up the inverse of a AR-1
					  correlation matrix  with rho */
SIZE_T ar(double **A,SIZE_T n,double rho); /* sets up a AR-1 correlation matrix 
				       with rho */
double quadratic(double **A,double *x,SIZE_T p); 
     /*calculates x'Ax where A is a p x p-dimensional matrix */
double ltquadratic(double *A,double *x,SIZE_T p); 
     /*calculates x'Ax where A is a p x p-dimensional symmetric matrix in 
      packed lower-triangular form*/
void cpy(double **a,SIZE_T nrows,SIZE_T ncols,double **b);

double L2norm(SIZE_T n, double *y);
long double longL2norm(SIZE_T n, long double *y);

void ltmatxvec(double *a, SIZE_T p, double *x, double *ax);
void longltmatxvec(long double *a, SIZE_T p, long double *x, long double *ax);

double vecxvec(double *x, SIZE_T p, double *y);
long double longvecxvec(long double *x, SIZE_T p, long double *y);
double *XprimeX(double **x, SIZE_T n, SIZE_T p);
double *indexed_XprimeX(double **x, SIZE_T n, SIZE_T p, SIZE_T *index);
double **aprimeb(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q);
void AprimeB(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q, double **c);
double **aprimebprime(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q);
double **aprimebprimeprime(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q);
double **abprimeprime(double **a, double **b, SIZE_T n, SIZE_T p, SIZE_T q);

SIZE_T logLTdeterminant(double *LTA, SIZE_T p, double *logdet);

double *Linverse(SIZE_T p, double *L);

double *LL_multiply(SIZE_T p, double *L, double *J);

double *LLT(double *l, SIZE_T p);

SIZE_T cholesky_inverse(double *A, SIZE_T n, double *Ainv);


#endif /* MATVEC_H */

