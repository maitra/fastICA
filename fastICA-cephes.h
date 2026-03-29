#include <stdlib.h>
#include "mat_vec.h"
#include "cephes_eigens.h"
#include "util.h"
#include "array.h"

typedef struct icaobject {
	SIZE_T numobs; //number of observatios in dataset
	SIZE_T numcomp; // number of independent components
	double **Xpre; // pre-whitened data of dimension n x ncomp
	double **K; // pre-whitening matrix of dimension ncomp x ncomp
	double **W; // final unmixing matrix of dimension ncomp x ncomp
	double **A; // estimated mixing matrix of dimension ncomp x ncomp
	double **S; // ICA components, of dimension n x ncomp
} ICAObject;


ICAObject
fastICA (double **X, double **Winit, SIZE_T n, SIZE_T p, SIZE_T ncomp,	
	 double alpha, SIZE_T centerflag, SIZE_T scaleflag, 
	 double (*fun)(double, double), double (*funprime)(double, double),
	 SIZE_T maxit, double tol, SIZE_T deflate, SIZE_T verbose);

double logcosh(double, double);

double logcoshprime(double, double);

double exponential(double, double);

double exponentialprime(double, double);
