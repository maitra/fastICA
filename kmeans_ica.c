#include "fastICA.h"

/** Implement Onoda et al (2012) method for K-means initialization (Ma09).
 *
 * @param x data
 * @param n number of observations
 * @param p number of dimensions
 * @param nclus number of desired clusters
 * @param return seeds (nclus x p)
 * @param alpha, set to 1
 * @param deflate, 1 (for deflation only option when nclus > m, 0 for 
 *        symmetry.
 */

void kmeans_init_ica(double **x, SIZE_T n, SIZE_T m, SIZE_T nclus, 
		    double **Mu, double alpha, SIZE_T deflate)
{
	SIZE_T i, j, k, vrbse = 0;
	ICAObject results;
	double **w_in, *xnorm, tollim = 1e-4;

	MAKE_MATRIX(w_in, nclus, nclus);

	for (i = 0; i < nclus; i++) 
		for (j = 0; j < nclus; j++)
			w_in[i][j] = rnorm(0.0, 1.0); /* initialize input unmixing matrix */
	
	if (nclus > m) 
		deflate = 1; /*symmetry method does not work here*/

	results = fastICA(x, w_in, n, m, nclus, alpha, 1, 1, logcosh, 
			  logcoshprime, 200, tollim, deflate, vrbse);
	
	FREE_MATRIX(w_in);

	MAKE_VECTOR(xnorm, n);
	for (i = 0; i < n; i++)
		xnorm[i] = L2norm(m, x[i]);
	
	for (i = 0; i < nclus; i++) {
		double dmax = fabs(results.S[0][i])/xnorm[0];
		k = 0;
		for (j = 1; j < n; j++) {
			double temp = fabs(results.S[j][i])/xnorm[j];
			if (temp > dmax) {
				k = j;
				dmax = temp;
			}
		}
		for (j = 0; j < m; j++)
			Mu[i][j] = x[k][j];
	}
	FREE_VECTOR(xnorm);
	
	FREE_VECTOR(results.Xpre);
	FREE_VECTOR(results.K);
	FREE_VECTOR(results.W);
	FREE_VECTOR(results.A);
	FREE_VECTOR(results.S);
}
	



