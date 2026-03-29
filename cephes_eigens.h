#ifndef CEPHES_EIGENS_H
#define CEPHES_EIGENS_H

#include <stdlib.h>
#include "constants.h"
#include "util.h"

void eigens(double *A, double *EV, double *E, size_t N);
void ordered_eigens(double *a, double **evec, double *eval, size_t p);
SIZE_T top_few_eigens(SIZE_T p, SIZE_T rankmax, double *lta, double **evecs, 
		      double *evals);

#endif /* CEPHES_EIGENS_H */
