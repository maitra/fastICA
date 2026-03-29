#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "mat_vec.h"
#include "util.h"

SIZE_T first_eigen(SIZE_T p, long double *A, long double *Y, SIZE_T maxiter, 
		   long double eps, long double *lambda);

SIZE_T kth_eigen(SIZE_T p, SIZE_T k, long double *A, long double **evecs, 
		 SIZE_T maxiter, long double eps, long double *evals);

void first_k_eigens(SIZE_T p, SIZE_T k, long double *A, long double **evecs, 
		    SIZE_T maxiter, long double eps, long double *evals);

SIZE_T first_rank_eigens(SIZE_T p, SIZE_T rankmax, long double *A, 
			 long double **evecs, long double *evals);

SIZE_T top_ranked_eigens(SIZE_T p, SIZE_T rankmax, double *A, double **evecs,
			 double *evals);
