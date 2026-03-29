#include "fastICA.h"

SIZE_T right_eigens(double **X, SIZE_T n, SIZE_T p, SIZE_T rankmax,
		    double **v, double *d)
{
	double *ltxtx;
	int i, j;
	
	ltxtx = XprimeX(X, n, p);

	i = first_rank_eigens(p, rankmax, ltxtx, v, d);
	for (j = 0; j < i; j++)
		d[j] = sqrt(d[j]); /* true eigenvalues are square root of Ds*/
	FREE_VECTOR(ltxtx);
	return i; /* rank of matrix, no point in looking for eigenvector-values
		     beyond this */
}

/*
  calculates the means and the sample standard deviations for a n x p data 
  matrix X.
  @param X: input data matrix of dimension n x p
  @param n: number of observations n
  @param p: number of coordinates/attributes
  @param mean: output vector of p means
  @param sd: output vector of p standard deviations
 */

void meansd(double **X, SIZE_T n, SIZE_T p, double *mean, double *sd)
{
	SIZE_T i, j;
	double tmp;
	
	for (i = 0; i < p; i++) {
		mean[i] = 0;
		tmp = 0;
		for (j = 0; j < n;  j++) 
			mean[i] += X[j][i];
		mean[i] /= n;
		for (j = 0; j < n;  j++)
			tmp += SQ(X[j][i] - mean[i]);
		tmp /= (n - 1);
		sd[i] = sqrt(tmp);
	}
}


void normalize(SIZE_T p, double *x)
{
	SIZE_T i;
	double sum = 0;

	for (i = 0; i < p; i++)
		sum += SQ(x[i]);
	if (sum != 0) 
		for (i = 0; i < p; i++)
			x[i] /= sqrt(sum);
}

void update_ws(SIZE_T ncomp, SIZE_T i, double *w, double **W) {
	double *t, K;
	SIZE_T j, l;

	if (i) { /* ignore for the first component */
		MAKE_VECTOR(t, ncomp);
		for (j = 0; j < ncomp; j++) 
			t[j] = 0;
		for (j = 0; j < i; j++) {
			K = 0;
			for (l  = 0; l < ncomp; l++)
				K += w[l] * W[l][j];
			for (l = 0; l < ncomp; l++)
				t[l] += K * W[l][j];
		}
		for (j = 0; j < ncomp; j++)
			w[j] -= t[j];
		FREE_VECTOR(t);
	}
	normalize(ncomp, w);
}

double logcosh(double x, double alpha)
{
	return tanh(alpha * x);
}

double logcoshprime(double x, double alpha)
{
	return alpha * (1 - SQ(tanh(alpha * x)));
}

double exponential(double x, double alpha)
{
	return x * exp(-SQ(x)/2);
}

double exponentialprime(double x, double alpha)
{
	return (1 - SQ(x)) * exp(-SQ(x)/2);
}

double **ica_by_deflation(double **X, SIZE_T n, SIZE_T ncomp, 
			  double (*f)(double, double), 
			  double (*fprime)(double, double), 
			  double tol, double alpha, SIZE_T maxit, 
			  SIZE_T verbose, double **Winit)
{
	SIZE_T i, j, k, iter;
	double *w, *wx, **W, *w1, *v1;

	if (verbose)
		printf("Deflation fastICA\n"); 

	MAKE_MATRIX(W, ncomp, ncomp);

	for (i = 0; i < ncomp; i++)
		for (j = 0 ; j < ncomp; j++)
			W[i][j] = 0;
	
	MAKE_VECTOR(w, ncomp);
	MAKE_VECTOR(w1, ncomp);
	MAKE_VECTOR(v1, ncomp);
	MAKE_VECTOR(wx, n);
	for (i = 0; i < ncomp; i++) {
		double lim = INFINITY;
		if (verbose) 
			printf("Component %zu \n", i);
		for (j = 0; j < ncomp; j++)
			w[j] = Winit[j][i];

		update_ws(ncomp, i, w, W);
		
		for (iter = 0; iter < maxit && lim > tol; iter++) {
			double meanfprimedotx = 0, sum = 0;
			matXvec(X, n, ncomp, w, ncomp, wx);
			for (j = 0; j < ncomp; j++) {
				v1[j] = 0;
				for (k = 0; k < n; k++)
					v1[j] += X[k][j] * f(wx[k], alpha);
				v1[j] /= n;
			}
			for (k  = 0; k < n; k++)
				meanfprimedotx += fprime(wx[k], alpha);
			meanfprimedotx /= n;

			for (j = 0; j < ncomp; j++) 
				w1[j] = v1[j] - w[j] * meanfprimedotx;

			update_ws(ncomp, i, w1, W);

			for (j = 0; j < ncomp; j++) 
				sum += w[j] * w1[j];
			lim = fabs(fabs(sum) - 1);
			if (verbose)
				printf("Iteration %zu, tolerance = %e\n", iter, lim);
			for (j = 0; j < ncomp; j++)
				w[j] = w1[j];
		}
		for (j = 0; j < ncomp; j++) 
			W[j][i] = w[j];
	}
	FREE_VECTOR(wx);
	FREE_VECTOR(w1);
	FREE_VECTOR(v1);
	FREE_VECTOR(w);
	return W;
}

void orthogonalize_matrix(double **W, SIZE_T m, SIZE_T p, SIZE_T maxrank,
			  double **V)
{
	/* Given a matrix W with singular value decomposition W = UDV', this
	   function calculates its orthogonalized version UV'. We do this here
	   without calculating the SVD but rather only the right eigenvectors
	   of W, and the singular values given in the diagonal matrix D. Then 
	   UV' = WVD^{-1}V', this orthogonalized version is returned.*/
	double *d, **temp1, **Vt;
	SIZE_T rank, i, j;

	MAKE_MATRIX(temp1, p, maxrank);
	MAKE_VECTOR(d, maxrank);
	rank = right_eigens(W, m, p, maxrank, V, d);
	matxmat(W, m, p, V, p, rank, temp1); // temp1 = WV, m x ncomp matrix

	for (i = 0; i < p; i++)
		for (j = 0; j < rank; j++)
			temp1[i][j] /= d[j];  //WVD^{-1}, m x ncomp matrix
	FREE_VECTOR(d);
	
	MAKE_MATRIX(Vt, rank, p);
	matrpose(V, p, rank, Vt);

	matxmat(temp1, p, rank, Vt, rank, m, V); //WVD^{-1}V'
	
	FREE_MATRIX(Vt);
	FREE_MATRIX(temp1);
}

double materrdiff(double **W, double **V, SIZE_T p)
{
	double err = 0;
	SIZE_T i, j;
	
	for (i = 0; i < p; i++) {
		double e, d = 0;
		for (j = 0; j < p; j++)
			d += W[i][j] * V[j][i];
		e = fabs(1 - fabs(d));
		if (e > err)
			err = e;
	}
	return err;
}

double **ica_in_symmetry(double **X, SIZE_T n, SIZE_T ncomp, SIZE_T rnk,
			 double (*f)(double, double), 
			 double (*fprime)(double, double), 
			 double tol, double alpha, SIZE_T maxit, 
			 SIZE_T verbose, double **Win)
{
	SIZE_T i, j, iter;
	double **W, **Winit, **gwx, **wx, **V, lim = INFINITY;
	
	if (verbose)
		printf("Symmetric fastICA using fixed point methods\n");

	MAKE_MATRIX(Winit, ncomp, ncomp);
	matrpose(Win, ncomp, ncomp, Winit);
	MAKE_MATRIX(gwx, n, ncomp);
	MAKE_MATRIX(W, ncomp, ncomp);

	orthogonalize_matrix(Winit, ncomp, ncomp, rnk, W);
	
	MAKE_MATRIX(wx, n, ncomp);
	MAKE_MATRIX(V, ncomp, ncomp);
	
	for (iter = 0; iter < maxit && lim > tol; iter++) {
		double meanfprimedotx;
		matxmat(X, n, ncomp, W, ncomp, ncomp, wx);

		for (i = 0; i < n; i++)
			for (j = 0; j < ncomp; j++)
				gwx[i][j] = f(wx[i][j], alpha);
		AprimeB(X, gwx, n, ncomp, ncomp, V); //V = X'gwx;

		for (i = 0; i < ncomp; i++)
			for (j = 0; j < ncomp; j++)
				V[i][j] /= n; // V = X'gwx/n;

		for (i = 0; i < ncomp; i++) {
			meanfprimedotx = 0;
			for (j = 0; j < n; j++)
				meanfprimedotx += fprime(wx[j][i], alpha);
			meanfprimedotx /= n;
			for (j = 0; j < ncomp; j++)
				gwx[j][i] = W[j][i] * meanfprimedotx;  
		}

		for (i = 0; i < ncomp; i++)
			for (j = 0; j < ncomp; j++)
				V[i][j] -= gwx[i][j];
		orthogonalize_matrix(V, ncomp, ncomp, rnk, Winit); //reuse Winit
		lim = materrdiff(W, Winit, ncomp);
		cpy(Winit, ncomp, ncomp, W);
	}
	FREE_MATRIX(Winit);
	FREE_MATRIX(wx);
	FREE_MATRIX(V);
	FREE_VECTOR(gwx);
	return W;
}

/*
  standardize data to have zero row means and unit row variances.
  @param n number of observations (rows of data matrix)
  @param p number of attributed (columns of data matrix)
  @param data n x p data matrix
  @center 0 if no centering, nonnull if data need to be centered
  @scale 0 if no scaling needed, nonnull if data are to be scaled.
  @param output scaleddata standardized data
*/
void standardize(SIZE_T n, SIZE_T p, double **data, SIZE_T center, 
		 SIZE_T scale,  double **scaleddata)
{
	SIZE_T i, j;
	double *mean, *sd;
	
	for (i = 0; i < n; i++)
		for (j = 0; j < p; j++)
			scaleddata[i][j] = data[i][j];
	if ((center) || (scale)) {
		MAKE_VECTOR(mean, p);
		MAKE_VECTOR(sd, p);
		meansd(scaleddata, n, p, mean, sd);
		for (i = 0; i < n; i++)
			for (j = 0; j < p; j++) {
				if (center) 
					scaleddata[i][j] -= mean[j];
				if (scale) 
					scaleddata[i][j] /= sd[j];
			}
		FREE_VECTOR(sd);
		FREE_VECTOR(mean);
	}
}

/** 
 * Performs fastICA.
 * 
 * @param X data matrix, n x p
 * @param Winit initial unmixing matrix, ncomp x ncomp (default: N(0,1) elements) 
 * @param n number of columns of matrix
 * @param p number of rows of matrix
 * @param ncomp number of independent components
 * @param alpha: constant in range [1, 2] used in approximation to neg-entropy
                 when fun == "logcosh", default = 1
 * @param rowflag: (nonnull) 1 when rows need to be centered, 0 otherwise
 * @param colflag: (nonnull) 1 when columns need to be standardized, 0 otherwise
 * @param funflag: 1(log cosh), 2 (exponential) approximation to neg-entropy fn
 * @param maxit: maximum number of iterations
 * @param tol: tolerance limit for convergence, 1e-4
 * @param defflag: 0 if components calculated in symmetry, nonnull uses deflation
 * @param verbose: nonnull to indicate verbosity
 * @return ICAobject structure with elements
 *                numobs (number of observations, same as n)
                  numcomp (number of components, same as ncomp)
		  Xpre: pre-whitened data matrix of dimension n x p
		  K: pre-whitening matrix, ncomp x ncomp
		  W: final unmixing matrix, ncomp x ncomp
		  A: estimated mixing matrix, ncomp x ncomp
		  S: ICA components, one in each column in n x ncomp matrix
*/

ICAObject
fastICA (double **X, double **Winit, SIZE_T n, SIZE_T p, SIZE_T ncomp,	
	 double alpha, SIZE_T centerflag, SIZE_T scaleflag, 
	 double (*fun)(double, double), double (*funprime)(double, double),
	 SIZE_T maxit, double tol, SIZE_T deflate, SIZE_T verbose)
{
	SIZE_T i, j, k, rnk;
	double *d, *dinv, **Xstdzd, **W;
	ICAObject res;
	
	MAKE_MATRIX(Xstdzd, n, p);
	
	if (verbose)
		printf("Standardizing data if called for\n");
	
	standardize(n, p, X, centerflag, scaleflag, Xstdzd);
	
	if (verbose)
		printf("Whitening:\n");
	
	MAKE_MATRIX(res.K, p, MAX(ncomp, p));
	MAKE_VECTOR(d, p);
	
	rnk = right_eigens(Xstdzd, n, p, MIN(n, p), res.K, d);
	
	if (ncomp > rnk) {
//		 ncomp = MIN(ncomp, rnk);
		if (verbose)
			printf("Rank of data matrix is only %zu, calculating %zu independent components using pseudo-inverses, etc\n", rnk, ncomp); 
		deflate = 1; /*use only deflation method*/
	}

	for (i = 0; i < p; i++) {
		for (j = 0; j < MIN(ncomp, rnk); j++)
			res.K[i][j] /= d[j]/sqrt(n);  //VD^{-1}, p x ncomp matrix
		for (j = MIN(ncomp, rnk); j < ncomp; j++)
			res.K[i][j] = 0;
	}
	FREE_VECTOR(d);
	
	res.Xpre = multiply(Xstdzd, n, p, res.K, p, ncomp);//XVD^{-1}, same as UD^{1}
	if (deflate)
		res.W = ica_by_deflation(res.Xpre, n, ncomp, fun, 
					 funprime, tol, alpha, maxit, verbose, 
					 Winit);
	else 
		res.W = ica_in_symmetry(res.Xpre, n, ncomp, rnk, fun, funprime, 
					tol, alpha, maxit, verbose, Winit);
	
	res.numobs = n;
	res.numcomp = ncomp;
	
	W = multiply(res.K, p, ncomp, res.W, ncomp, ncomp); //W = KW

	res.S = multiply(Xstdzd, n, p, W, p, ncomp); //S = XW

	d = XprimeX(W, p, ncomp);
	MAKE_VECTOR(dinv, ncomp * (ncomp + 1)/2);
	if (ncomp > rnk)
		symm_Moore_Penrose_inverse(ncomp, rnk, d, dinv);
	else 
		cholesky_inverse(d, ncomp, dinv);
	FREE_VECTOR(d);
	MAKE_MATRIX(res.A, ncomp, p);
	for (i = 0; i < ncomp; i++)
		for (j = 0; j < p; j++) 
			res.A[i][j] = 0;
	for (i = 0; i < ncomp; i++)
		for (j = 0; j < p; j++)
			for (k = 0; k < ncomp; k++)
				res.A[i][j] += dinv[LTINDEX(i, k)] * W[j][k]; 
	FREE_VECTOR(dinv);
	FREE_MATRIX(W);
	FREE_MATRIX(Xstdzd);
	return res;
}

	 
