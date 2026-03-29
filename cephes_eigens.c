/* cephes_eigens.c
 *
 * This file was orignially called eigens.c and was found in the cephese
 * computational library.  I edited it to bring it up to the ANSI C
 * standard and renamed it to cephes_eigens.c.
 *
 * --Rouben Rostamian
 * July 2002
 *
 * Comment added: The eigenvalues are surprisingly output in any order.
 * 
 * --Ranjan Maitra
 * April 2004
 *
 *
 * Changed to allow for reverse order of eigenvalues new function is called
 * ordered_eigens
 * 
 * --Ranjan Maitra
 * August 2011
 *
 * Eigenvalues and eigenvectors of a real symmetric matrix
 *
 *
 *
 * SYNOPSIS:
 *
 * SIZE_T n;
 * double A[n*(n+1)/2], EV[n*n], E[n];
 * void eigens( A, EV, E, n );
 *
 *
 *
 * DESCRIPTION:
 *
 * The algorithm is due to J. vonNeumann.
 *                   -     -
 * A[] is a symmetric matrix stored in lower triangular form.
 * That is, A[ row, column ] = A[ (row*row+row)/2 + column ]
 * or equivalently with row and column interchanged.  The
 * indices row and column run from 0 through n-1.
 *
 * EV[] is the output matrix of eigenvectors stored columnwise.
 * That is, the elements of each eigenvector appear in sequential
 * memory order.  The jth element of the ith eigenvector is
 * EV[ n*i+j ] = EV[i][j].
 *
 * E[] is the output matrix of eigenvalues.  The ith element
 * of E corresponds to the ith eigenvector (the ith row of EV).
 *
 * On output, the matrix A will have been diagonalized and its
 * orginal contents are destroyed.
 *
 * ACCURACY:
 *
 * The error is controlled by an internal parameter called RANGE
 * which is set to 1e-10.  After diagonalization, the
 * off-diagonal elements of A will have been reduced by
 * this factor.
 *
 * EEVOR MESSAGES:
 *
 * None.
 *
 */

/*
Copyright 1973, 1991 by Stephen L. Moshier
Copyleft version.
*/

#include <math.h>
#include "cephes_eigens.h"
#include "order.h"
#include "array.h"

void eigens(double *A, double *EV, double *E, SIZE_T N)
{
	SIZE_T IND, L, LL, LM, M, MM, MQ, I, J, IA, LQ;
	SIZE_T IQ, IM, IL, NLI, NMI;
	double ANORM, ANORMX, AIA, THR, ALM, ALL, AMM, X, Y;
	double SINX, SINX2, COSX, COSX2, SINCS, AIL, AIM;
	double RLI, RMI;
	static double RANGE = 1.0e-10;	/* 3.0517578e-5; */

	/* Initialize identity matrix in EV[] */
	for (J = 0; J < N * N; J++)
		EV[J] = 0.0;

	MM = 0;
	for (J = 0; J < N; J++) {
		EV[MM + J] = 1.0;
		MM += N;
	}

	ANORM = 0.0;
	for (I = 0; I < N; I++) {
		for (J = 0; J < N; J++) {
			if (I != J) {
				IA = I + (J * J + J) / 2;
				AIA = A[IA];
				ANORM += AIA * AIA;
			}
		}
	}

	if (ANORM <= 0.0)
		goto done;

	ANORM = sqrt(ANORM + ANORM);
	ANORMX = ANORM * RANGE / N;
	THR = ANORM;

	while (THR > ANORMX) {
		THR = THR / N;

		do {		/* while IND != 0 */
			IND = 0;

			for (L = 0; L < N - 1; L++) {

				for (M = L + 1; M < N; M++) {
					MQ = (M * M + M) / 2;
					LM = L + MQ;
					ALM = A[LM];
					if (fabs(ALM) < THR)
						continue;

					IND = 1;
					LQ = (L * L + L) / 2;
					LL = L + LQ;
					MM = M + MQ;
					ALL = A[LL];
					AMM = A[MM];
					X = (ALL - AMM) / 2.0;
					Y = -ALM / sqrt(ALM * ALM + X * X);
					if (X < 0.0)
						Y = -Y;
					SINX = Y/sqrt(2.0*(1.0+sqrt(1.0-Y*Y)));
					SINX2 = SINX * SINX;
					COSX = sqrt(1.0 - SINX2);
					COSX2 = COSX * COSX;
					SINCS = SINX * COSX;

					/* ROTATE L AND M COLUMNS */
					for (I = 0; I < N; I++) {
						IQ = (I * I + I) / 2;
						if ((I != M) && (I != L)) {
							if (I > M)
								IM = M + IQ;
							else
								IM = I + MQ;

							if (I >= L)
								IL = L + IQ;
							else
								IL = I + LQ;

							AIL = A[IL];
							AIM = A[IM];
							X = AIL * COSX - AIM * SINX;
							A[IM] = AIL * SINX + AIM * COSX;
							A[IL] = X;
						}

						NLI = N * L + I;
						NMI = N * M + I;
						RLI = EV[NLI];
						RMI = EV[NMI];
						EV[NLI] = RLI * COSX - RMI * SINX;
						EV[NMI] = RLI * SINX + RMI * COSX;
					}

					X = 2.0 * ALM * SINCS;
					A[LL] = ALL * COSX2 + AMM * SINX2 - X;
					A[MM] = ALL * SINX2 + AMM * COSX2 + X;
					A[LM] = (ALL - AMM) * SINCS + ALM * (COSX2 - SINX2);
				}	/* for M=L+1 to N-1 */
			}	/* for L=0 to N-2 */

		} while (IND != 0);
	}			/* while THR > ANORMX */

done:

	/* Extract eigenvalues from the reduced matrix */
	L = 0;
	for (J = 1; J <= N; J++) {
		L = L + J;
		E[J - 1] = A[L - 1];
	}
	return;
}

void ordered_eigens(double *a, double **evec, double *eval, SIZE_T p)
{
	double *u, *d;
	SIZE_T i, j;
	SIZE_T *ind;

	MAKE_VECTOR(u, p * p);
	MAKE_VECTOR(d, p);

	eigens(a, u, d, p);
	
	ind = orderDouble(d, p);

	/*for (i = (p - 1); i >=0; i--) {*/
	for (i = 0; i < p; i++) {
		eval[p - i - 1] = d[ind[i]];
		for (j = 0; j < p; j++) 
			evec[j][p - i - 1] = u[j + p * ind[i]];
	}
	
	FREE_VECTOR(ind);

	FREE_VECTOR(d);	
	FREE_VECTOR(u);
	
}
	
SIZE_T top_few_eigens(SIZE_T p, SIZE_T rankmax, double *lta, double **evecs, 
		      double *evals)
{
	double trc = 0, proptr = 0;
	SIZE_T rankmatrix;

	ordered_eigens(lta, evecs, evals, p); 

	for (rankmatrix = 0; rankmatrix < p; rankmatrix++) 
		trc += lta[LTINDEX(rankmatrix, rankmatrix)];
	
	for (rankmatrix = 0; ((proptr <= 0.995) && (rankmatrix < MIN(p, rankmax))) ; rankmatrix++) 
		proptr += fabs(evals[rankmatrix]) / trc;
	return rankmatrix;
}
