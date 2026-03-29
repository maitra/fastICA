#ifndef _UTIL_H
#define _UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include "constants.h"
#include "kmeans.h"

/**
 * Simple mathematical functions.
 */
#define SQ(x) ((x) * (x)) 
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQ(x) ((x) * (x))
#define LTINDEX(k, l) ((MAX(k, l) * (MAX(k, l) + 1))/2 + MIN(k, l))
#define LLTINDEX(k, l) (LTINDEX(k - 1, l))
enum {NO_ERROR,
	MEMORY_ALLOCATION_ERROR,	/*!< failure to allocate memory */
	INVALID_COMMAND_LINE,		/*!< user misuse of command line */
	FILE_OPEN_ERROR,		/*!< failure to open file */
	FILE_FORMAT_ERROR,		/*!< unexpected file format */
	FILE_WRITE_ERROR,		/*!< failure to write to file */
	INTERNAL_ERROR,			/*!< any other error */
	NUMBER_OF_ERRORS
};

#define NO_EXIT 0	/* do not let called function exit */
#define DO_EXIT 1	/* induce a called function to exit */

extern const char *error_string[];


/**
 * Convert upper triangle symmetric matrix indices to vector index.
 *
 * @param n size of matrix
 * @param i row index
 * @param j index, must have j>i
 * @return index of same position in row-order vector
 */
#define VINDEX(n, i, j) (                                     \
	(i) < (j)                                             \
		? (j) + (i) * (n) - ((i) + 1) * ((i) + 2) / 2 \
		: (i) + (j) * (n) - ((j) + 1) * ((j) + 2) / 2 \
)

/**
 * Generic matrix object.  When you need to operate on a matrix using generic
 * functions, like random samplers, that use call-backs to check user-required
 * conditions, can pass this object around as a void pointer.  It currently
 * just provides the data and the number of columns.
 */
struct data_obj {
	double **x;	/*!< data matrix */
	SIZE_T p;	/*!< number of columns in matrix */
};


/* distance calculations */
double sqdist(double *x, double *y, SIZE_T p);
double dist(double *x, double *y, SIZE_T p);
double sq_norm(double *, SIZE_T);

SIZE_T minindex(double *x, SIZE_T n);
int unique_points(SIZE_T *y, SIZE_T i, void *d);
int shuffle(SIZE_T *array, SIZE_T n);

/* k-means-related */
int check_wss(SIZE_T, SIZE_T, double *, double *, double *, SIZE_T *, SIZE_T *, double *, SIZE_T *, SIZE_T *);
void reset_wss(double *w, SIZE_T k);
SIZE_T null_cluster_downgrade(double **x, double **centers, SIZE_T n, SIZE_T p, SIZE_T k, double *dis, SIZE_T *sd_idx, int free);

double elapsed_seconds(struct timeval *x, struct timeval *y);

int error(int errnum, int do_exit, const char *filename, const char *fxn_name, int line_no, const char *message, ...);

char *string_copy(const char *);

/* io-related */
int check_newline(FILE *finp, const char *filename, const char *file, const char *fxn_name, int line_no);
char *make_full_filename(const char *path, const char *name);
const char *remove_path(const char *fullPath);	/* from optlist version 3 */


#endif
