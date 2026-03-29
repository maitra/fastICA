#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "util.h"



/**
 * Determine if a selected rows in a matrix is distinct from previously
 * selected rows.
 *
 * @param y indices of selected rows
 * @param i index of newest row
 * @param d data object with data matrix and number of columns
 * @return 1/0 for TRUE/FALSE unique
 */
int
unique_points(SIZE_T *y, SIZE_T i, void *d)
{
	struct data_obj *dobj = (struct data_obj *)d;
	SIZE_T j;

	for (j = 0; j < i; j++) {
		if (sqdist(dobj->x[y[i]], dobj->x[y[j]], dobj->p)
			< DBL_EPSILON)
			return 0;
	}

	return 1;
} /* unique_points */


/**
 * Generate random permutation of set of indices.  Note, will work with size_t
 * indices so long as they do not exceed RAND_MAX.  It issues an error if this
 * constraint is not satisfied.
 *
 * [TODO] Write non-rand() version that is not restricted to type int.
 *
 * @param indices indices to permute
 * @param n length of array
 * @return error condition
 */
int
shuffle(SIZE_T *indices, SIZE_T in_n) {
	SIZE_T l;
	int n = (int) in_n;
	int i;
	int j;
	int lim = RAND_MAX - RAND_MAX % n;

	/* check error conditions: cannot handle super large datasets */
	if (in_n > RAND_MAX)
		return error(INTERNAL_ERROR, NO_EXIT, __FILE__, __func__,
			__LINE__, "Permutation code only works for datasets "
			"with fewer than %u observations.",
			MIN(INT_MAX, RAND_MAX));


	for (i = n - 1; i > 0; i--) {

		/* generate integer ~ Unif(0, i) without modulus bias */
		lim = RAND_MAX - RAND_MAX % (i + 1);
		do {
			j = rand();
		} while (j >= lim);
		j = j % (i + 1);

		/* swap array elements */
		l = indices[j];
		indices[j] = indices[i];
		indices[i] = l;
	}

	return NO_ERROR;
} /* shuffle */




/**
 * Human description of errors.
 */
const char *error_string[NUMBER_OF_ERRORS] = {
	"Warning",
	"Memory allocation error",
	"Invalid command line",
	"File not found",
	"Invalid input file format",
	"File write error",
	"Internal error",
};


/**
 * Compute squared Euclidean distance between two p-vectors.
 *
 * @param x vector of length p
 * @param y vector of length p
 * @param p length of vector
 * @return squared distance between two vectors
 */
double
sqdist(double *x, double *y, SIZE_T p)
{
	double sum = 0;
	SIZE_T i;

	for (i = 0; i < p; i++)
		sum += SQ(x[i] - y[i]);
	return sum;
} /* sqdist */


/** Compute squared norm of a p-vector.
 *
 * @param vec array of length p
 * @param p length of vector
 * @return squared norm
 */
double
sq_norm(double *vec, SIZE_T p)
{
	SIZE_T i;
	double sum = 0;
	for (i = 0; i < p; i++)
		sum += vec[i] * vec[i];
	return sum;
} /* sq_norm */


/**
 * Compute Euclidian distance between two p-vectors.
 *
 * @param x vector of length p
 * @param y vector of length p
 * @param p length of vector
 * @return distance between two vectors
 */
double
dist(double *x, double *y, SIZE_T p)
{
	return sqrt(sqdist(x, y, p));
} /* dist */


/**
 * Find the index of the minimum element in an array.
 *
 * @param x array of doubles
 * @param length of array
 * @return index of minimum value
 */
SIZE_T
minindex(double *x, SIZE_T n)
{
	SIZE_T i, s1 = 0;
	double d1 = x[0];

	for (i = 1; i < n; i++)
		if (x[i] < d1) {
			d1 = x[i];
			s1 = i;
		}

	return s1;
} /* minindex */


/**
 * Propose a candidate and store a better K-means solution.
 *
 * This auxiliary function checks if the K-means solution in wss, iclass,
 * nclass is better than the one that yielded min_wss.  If so, udate min_wss
 * and return TRUE.  Also, if the solution is better than global_min_wss, update
 * global_min_wss and store solution in bestwss, besticlass, and bestnclass.
 * The point of the double check is to work within an iteration of an iterative
 * algorithm, like HK05.
 *
 * If any of bestwss, besticlass, or bestnclass is NULL, do not copy.
 *
 * @param n number of observations
 * @param k number of clusters
 * @param global_min_wss current global minimum within-sum-of-squares (may be updated)
 * @param min_wss current minimum within-sum-of-squares (may be updated)
 * @param wss candidate within-sum-of-squares per cluster
 * @param iclass candidate assignment of observations to clusters
 * @param nclass candidate cluster sizes
 * @param bestwss current best wss
 * @param besticlass current best iclass
 * @param bestnclass current best nclass
 * @param return 1 if a better solution is found
 */
int
check_wss(SIZE_T n, SIZE_T k, double *global_min_wss, double *min_wss,
	double *wss, SIZE_T *iclass, SIZE_T *nclass, double *bestwss,
	 SIZE_T *besticlass, SIZE_T *bestnclass)
{
	SIZE_T i;
	double total_wss = 0;

	/* compute proposed within-sum-of-squares */
	for (i = 0; i < k; i++)
		total_wss += wss[i];

	/* this is a better solution */
	if (global_min_wss && *global_min_wss > total_wss) {
		*global_min_wss = total_wss;
		if (min_wss)
			*min_wss = total_wss;
		/* copy solution */
		if (besticlass)
			for (i = 0; i < n; i++)
				besticlass[i] = iclass[i];
		for (i = 0; i < k; i++) {
			if (bestnclass)
				bestnclass[i] = nclass[i];
			if (bestwss)
				bestwss[i] = wss[i];
		}
		return 1;
	} else if (min_wss && *min_wss > total_wss) {
		*min_wss = total_wss;
		return 0;
	}
	return 0;
} /* check_wss */


/**
 * Error function, prints message and exits.
 *
 * @param errnum error number
 * @param message error message
 * @return error number
 */
int
error(int errnum, int do_exit, const char *filename, const char *fxn_name,
	int line_no, const char *message, ...)
{
	va_list args;

	errnum
		? fprintf(stderr, "\n[ERROR]")
		: fprintf(stderr, "[DEBUG]");
	fprintf(stderr, " at %s:%s():%d: ", filename, fxn_name, line_no);
	if (message) {
		va_start(args, message);
		vfprintf(stderr, message, args);
		va_end(args);
		if (errnum < NUMBER_OF_ERRORS)
		        fprintf(stderr, " (%s)\n", error_string[errnum]);
		else if (errno)
			fprintf(stderr, " (%s)\n", strerror(errno));
	} else if(errnum < NUMBER_OF_ERRORS)
	        fprintf(stderr, " %s.\n", error_string[errnum]);
	else if(errno)
		fprintf(stderr, " %s.\n", strerror(errno));
	else
		putc('\n', stderr);
	fflush(stderr);
	if (do_exit)
		exit(errnum);
	return errnum;
} /* error */


/**
 * Reset within-sums-of-squares.
 */
void
reset_wss(double *w, SIZE_T k)
{
	SIZE_T i;
	for (i = 0; i < k; i++)
		w[i] = DBL_MAX;
} /* reset_wss */


/**
 * Return filename with path.
 * Note, this function allocates memory for the string, and the caller
 * is responsible for freeing the memory.
 *
 * @param path path to file
 * @param name name of file
 * @return string containing full filename with path.
 */
char *
make_full_filename(const char *path, const char *name)
{
	const char *fxn_name = "make_full_filename";
	char *full_filename = malloc(((path ? strlen(path) : 0)
		+ strlen(name) + 2) * sizeof *full_filename);
	if (!full_filename) {
		free(full_filename);
		error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__, fxn_name,
			__LINE__, NULL);
		return NULL;
	}
	strcpy(full_filename, "");
	if (path) {
		strcpy(full_filename, path);
		strcat(full_filename, "/");
	}
	strcat(full_filename, name);
	return full_filename;
} /* make_full_filename */


/**
 * Check proposed seeds against a dataset for null clusters and remove null
 * cluster seeds.
 *
 * @param x dataset
 * @param seed proposed seeds, reordered to have non=null seeds first
 * @param n number of observations
 * @param p dimension of dataset
 * @param K number of seeds
 * @param dis distance matrix (upper triangle in vector) [optional]
 * @param sd_idx index of seeds in data [optional]
 * @param key (if 1, then free dis and sd_idx; if 2, then only check first seed)
 * @return number of remaining seeds
 */
SIZE_T
null_cluster_downgrade(double **x, double **seed, SIZE_T n, SIZE_T p, SIZE_T K,
	double *dis, SIZE_T *sd_idx, int key)
{
	double d, *dptr;
	SIZE_T i, k, s=0;
	double dmin;		/* minimum distance to seed */
	SIZE_T newK;		/* number of non-null clusters */
	SIZE_T *nc = NULL;	/* number of observations per cluster */
	int err = NO_ERROR;	/* error status */

	newK = K;

	/* we only need to know that the first cluster is non-null */
	if (key == 2)
		for (i = 0; i < n; i++) {
			dmin = INFINITY;
			for (k = 0; k < K; k++) {
				if (dis && sd_idx)
					d = dis[VINDEX(n, i, sd_idx[k])];
				else
					d = sqdist(x[i], seed[k], p);
				if (d < dmin) {
					dmin = d;
					s = k;
				}
			}
			if (!s)
				goto RETURN;
		}

	CMAKE_VECTOR(nc, K);
	if (!nc) {
		err = MEMORY_ALLOCATION_ERROR;
		goto RETURN;
	}

	/* count number in each cluster */
	for (i = 0; i < n; i++) {
		dmin = INFINITY;
		for (k = 0; k < K; k++) {
			if (dis && sd_idx)
				d = dis[VINDEX(n, i, sd_idx[k])];
			else
				d = sqdist(x[i], seed[k], p);

			if (d < dmin) {
				dmin = d;
				s = k;
			}
		}
		nc[s]++;
	}

	/* check for null clusters */
	for (k = 0; k < newK; k++) {
		if (!nc[k]) {
			do {
				newK--;
			} while (k < newK && !nc[newK]);
			if (newK == k)
				break;
			s = nc[newK];
			nc[newK] = nc[k];
			nc[k] = s;
			dptr = seed[newK];
			seed[newK] = seed[k];
			seed[k] = dptr;
		}
	}

RETURN:
	FREE_VECTOR(nc);
	if (key == 1) {
		if (dis)
			FREE_VECTOR(dis);
		if (sd_idx)
			FREE_VECTOR(sd_idx);
	}
	return err ? 0 : newK;
} /* null_cluster_downgrade */


/**
 * Compute the seconds between two processor time values.
 * From Gnu.org.
 */
double
elapsed_seconds(struct timeval *x, struct timeval *y)
{
	long int esec, eusec;
	if (y->tv_usec < x->tv_usec) {
		int nsec = (x->tv_usec - x->tv_usec) / 1000000 + 1;
		x->tv_usec -= 1000000 * nsec;
		x->tv_sec += nsec;
	}
	if (y->tv_usec - x->tv_usec > 1000000) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}
	esec = y->tv_sec - x->tv_sec;
	eusec = y->tv_usec - x->tv_usec;
	return esec + (double)eusec/1000000;
} /* elapsed_seconds */ 


/**
 * Verify there are no non-blank characters before next newline.
 * It is a reasonable way to check for the expected file format.
 *
 * @param finp open file for reading
 * @param filename name of file opened for reading
 * @param file name of coding file calling this function
 * @param fxn_name name of function calling this funciton
 * @param line_no line number calling this function
 * @return 1 if there is an error
 */
int
check_newline(FILE *finp, const char *filename, const char *file,
	 const char *fxn_name, int line_no)
{
	char c;
	do {
		if (!fscanf(finp, "%c", &c)) {
			error(FILE_FORMAT_ERROR, 0, file, fxn_name, line_no,
				"Failed to read data file '%s'", filename);
			return 1;
		}
		if (c != ' ' && c != '\n') {
			error(FILE_FORMAT_ERROR, 0, file, fxn_name, line_no,
				"Invalid format in data file '%s' (%c).  Are your dimensions correct?",
				filename, c);
			return 1;
		}
	} while (c != '\n');
	return 0;
} /* check_newline */


/**
 * Remove path from a filename including full path.
 *
 * This function is from the optlist library, available on the
 * <a href="http://michael.dipperstein.com/optlist/">web</a>.
 *
 * @param full_path name of file including full path
 * @return name of file without full path
 */
const char *
remove_path(const char *full_path)
{
	int i;
	const char *start, *tmp;		/* start of file name */
	const char delim[3] = {'\\', '/', ':'};	/* path deliminators */

	start = full_path;

	/* find the first character after all file path delimiters */
	for (i = 0; i < 3; i++) {
		tmp = strrchr(start, delim[i]);

		if (tmp != NULL)
			start = tmp + 1;
	}

	return start;
} /* remove_path */


/**
 * Copy a string, allocating space.  The caller is responsible for 
 * freeing the memory.
 *
 * @param str the string to copy
 * @return newly allocated string
 */
char *
string_copy(const char *str)
{
	char *str_cpy = malloc((strlen(str)+1) * sizeof *str_cpy);
	strcpy(str_cpy, str);
	return str_cpy;
} /* string_copy */
