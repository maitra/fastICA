#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "mat_vec.h"
#include "util.h"

SIZE_T first_eigen(SIZE_T p, double *A, double *Y, SIZE_T maxiter, double eps,
		   double *lambda);

SIZE_T kth_eigen(SIZE_T p, SIZE_T k, double *A, double **evecs, SIZE_T maxiter, 
		 double eps, double *evals);

void first_k_eigens(SIZE_T p, SIZE_T k, double *A, double **evecs, 
		    SIZE_T maxiter, double eps, double *evals);

SIZE_T first_rank_eigens(SIZE_T p, SIZE_T rankmax, double *A, double **evecs,
			 double *evals);

SIZE_T symm_Moore_Penrose_inverse(SIZE_T p, SIZE_T rankmax, double *A, double *Ainv);
