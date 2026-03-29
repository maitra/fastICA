#include "power_eigens_long.h"

/* Calculates the dominant (meaning largest in magnitude) eigenvalue and 
   the corresponding eigenvector using the power method. 

Discussion:

For a given lower triangular symmetric n x n matrix A (with lower triangular 
stored in the form of a vector of length n*(n+1)/2, the power method 
iterates through a series of estimates to converge at the largest (in magnitude)
eigenvalue: the corresponding vector is the eigenvector. 


The iteration repeats the following steps
AY     = A * Y
LAMBDA = || AY ||
Y      = AY / LAMBDA

If the matrix A has a single real eigenvalue of maximum modulus, then this 
iteration will generally produce a good estimate for that eigenvalue and its
corresponding eigenvector.

If there are multiple distinct eigenvalues of the same modulus, perhaps two 
values of opposite sign, or complex eigenvalues, then the situation is more 
complicated and this method should not be pursued.

Separate issues:

* when estimating the value of LAMBDA, we use the Rayleigh quotient,
LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
bottom of the fraction is 1.  Using this estimate allows us to
easily capture the sign of LAMDBA.  Using the eucldean norm
instead, for instance, would always give a positive value.

* If the dominant eigenvalue is negative, then the iteration 
as given will produce eigenvector iterates that alternate in sign.  
   
* It is worth knowing whether the successive eigenvector estimates
are tending to some value.  Since an eigenvector is really a direction,
we need to normalize the vectors, and we need to somehow treat both
a vector and its negative as holding the same information.  This
means that the proper way of measuring the difference between two
eigenvector estimates is to normalize them both, and then compute
the cosine between them as y1'y2, followed by the sine, which is
sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
are "close" in the sense of direction.

Parameters:
    
    Input, SIZE_T p, the order of the matrix.
    
    Input, long double A, the matrix in lower-triangular form.

    Input/output, long double Y, the estimate for the eigenvector.

    Input, int maxiter, the maximum number of iterations allowed.

    Input, long double eps, an error tolerance.

    Output, long double *lambda, the estimate for the eigenvalue.

    Value, 0 indicates convergence, 1 otherwise

Copyright Ranjan Maitra, 4 January 2013.

*/
/******************************************************************************/

SIZE_T first_eigen_long(SIZE_T p, long double *A, long double *Y, 
			SIZE_T maxiter, long double eps, long double *lambda)
{
	long double *AY, siny1y2, norm, dif = 1e+100, *Ycurr, lambdacurr;
	SIZE_T i, iter;

	MAKE_VECTOR(Ycurr, p);
	MAKE_VECTOR(AY, p);
	
	for (iter = 0; ((iter < maxiter) && (dif > eps)); iter++) { //main loop
		lambdacurr = *lambda;
		for ( i = 0; i < p; i++) 
			Ycurr[i] = Y[i];
		longltmatxvec(A, p, Y, AY); 

		*lambda = longvecxvec(Y, p, AY);

//		printf("l = %f", *lambda);

		norm = longL2norm(p, AY);
		for (i = 0; i < p; i++)
			Y[i] = AY[i] / norm;
		if (*lambda < 0) {
			for ( i = 0; i < p; i++ )
				Y[i] = - Y[i];
		}
		dif = fabs(*lambda - lambdacurr);
		siny1y2 = sqrt(1 - SQ(longvecxvec (Y, p, Ycurr)));
//		printf("%e \n", siny1y2);
	}
	
	for (i = 0; i < p; i++ )
		Y[i] = AY[i] / *lambda;
  
	FREE_VECTOR(AY);
	FREE_VECTOR(Ycurr);
	if (iter >= maxiter) 
		return 1;
	else 
		return 0;
}

SIZE_T kth_eigen_long(SIZE_T p, SIZE_T k, long double *A, long double **evecs, 
		      SIZE_T maxiter, long double eps, long double *evals)
{
	/* Given (k - 1) eigenvalue-eigenvector combinations, this function
	   provides the kth dominant eigenvalue-eigenvector combination of a 
	   symmetric matrix A stored in lower-triangular form.

	   The kth eigenvalue is the dominant eigenvalue of the matrix
	   A - \sum_{l=1}^{k-1} lambda_l v_l v_l' where v_l is the eigenvector
	   corresponding to the lth eigenvalue.	  
	   
	   INPUT: 
	   p: dimension of the matrix
	   k: the eigenvalue-eigenvector combination needed
	   A: the lower triangle of the symmetric matrix A
	   evecs: the matrix of eigenvectors (ith column storing the
	          eigenvector corresponding to the ith eigenvalue). Storage
		  of p pointers to pointers with at least k elements is
		  provided. In matrix parlance, we have a matrix of p rows
		  and at least k columns, of which the first k-1 columns 
		  contain the eigenvectors corresponding to the (k-1) dominant
		  (largest in magnitude) eigenvalues.
	   evals: the vector of eigenvalues (pointer to at least k elements) 
	          with the first (k-1) elements containing the (k-1) dominant
		  eigenvalues.
           maxiter, eps, are as in first_eigen.
	   
	   OUTPUT:
	   evecs: (in matrix parlance, the kth column has the eigenvector 
	           corresponding to the kth eigenvalue).
	   evals: the kth cell returns the kth dominant eigenvalue.
	*/

	long double *y, nrm;
	SIZE_T i, j, l, ind;
	
/*	for (i = 0; i < p * (p + 1)/2; i++)
		printf("A = %f ", A[i]);
		printf("\n");*/

	MAKE_VECTOR(y, p);

	for (i = 0; i < p; i++) 
		y[i] = A[LTINDEX(i, k - 1)]; /* initial guess for eigenvector y, the
				     kth column of A*/
	nrm = longL2norm(p, y); 
	for (i = 0; i < p; i++) /*normalize the initial guess*/
		y[i] /= nrm;
	evals[k - 1] = nrm; /* initial guess for the eigenvalue*/
	
	if (k == 1)  
		ind = first_eigen_long(p, A, y, maxiter, eps, &evals[0]);
	else 	{
		long double *a;
		MAKE_VECTOR(a, p * (p + 1)/2);
		
		for (i = 0; i < p * (p + 1)/2; i++)
			a[i] = A[i];

		for (i = 0; i < p; i++) 
			for (j = 0; j <= i; j++) 
				for (l = 0; l < (k - 1); l++) 
					a[LTINDEX(i, j)] -= evals[l] * evecs[i][l] * evecs[j][l];

/*		for (i = 0; i < p * (p + 1)/2; i++)
			printf("a=%f ",a[i]);
			printf("\n");*/
		
		ind = first_eigen_long(p, a, y, maxiter, eps, &evals[k-1]);
		
		FREE_VECTOR(a);
	}
	
	for (i = 0; i < p; i++) 
		evecs[i][k-1] = y[i];

	FREE_VECTOR(y);
	return ind;
}

void first_k_eigens_long(SIZE_T p, SIZE_T k, long double *A, 
			 long double **evecs, SIZE_T maxiter, long double eps, 
			 long double *evals)
{
	int i;
	for (i = 1; i <= k; i++)
		kth_eigen_long(p, i, A, evecs, maxiter, eps, evals);
}

SIZE_T first_rank_eigens_long(SIZE_T p, SIZE_T rankmax, long double *A, 
			      long double **evecs, long double *evals)
{
	/* this function calculates the positive eigenvalues (defined in terms
	   of those that explain 99.5% of the sum of the diagonals) of 
	   a p x p nonnegative definite matrix in lower triangular form. It
	   does not make much sense to use for other matrices. */
	long double trc = 0, proptr = 0;
	SIZE_T rankmatrix;

	for (rankmatrix = 0; rankmatrix < p; rankmatrix++) 
		trc += A[LTINDEX(rankmatrix, rankmatrix)];

	for (rankmatrix = 0; ((proptr <= 0.995) && (rankmatrix < MIN(p, rankmax))) ; rankmatrix++) {  
		/* dmaxmin is limited to dproj for computational practicality */
		if (kth_eigen_long(p, rankmatrix + 1, A, evecs, 1000000, 1e-16, evals))
			/* don't do anything about it*/
			printf("eigenvalue approximation did not converge\n");
		proptr += fabs(evals[rankmatrix]) / trc;
	}	
	return rankmatrix;
}

SIZE_T top_ranked_eigens(SIZE_T p, SIZE_T rankmax, double *A, double **evecs,
			 double *evals)
{
	/* this function calculates the positive eigenvalues (defined in terms
	   of those that explain 99.5% of the sum of the diagonals) of 
	   a p x p nonnegative definite matrix in lower triangular form. It
	   does not make much sense to use for other matrices. */
	long double *a, **ev, *d;
	SIZE_T i, j, rank;

	MAKE_VECTOR(a, p * (p + 1)/2);
	for (i = 0; i < p * (p + 1)/2; i++)
		a[i] = A[i];
	
	MAKE_MATRIX(ev, p, MIN(p, rankmax));
	MAKE_VECTOR(d, MIN(p, rankmax));
	
	rank = first_rank_eigens_long(p, rankmax, a, ev, d);

	for (i = 0; i < rank; i++) 
		evals[i] = d[i];
	
	for (i = 0; i < p; i++)
		for (j = 0; j < rank; j++)
			evecs[i][j] = ev[i][j];
	FREE_MATRIX(ev);
	FREE_VECTOR(d);
	FREE_VECTOR(a);
	return rank;
}

