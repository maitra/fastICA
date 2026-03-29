#ifndef KMEANS_H
#define KMEANS_H

#include "constants.h"
#include "srswor.h"
#include "util.h"

/**
 * Include macros for allocating vectors and matrices.  Also define base
 * allocators to handle errors according to the program's style.  There are two
 * versions, one that immediately returns with the MEMORY_ALLOCATION_ERROR and
 * another that that simply issues an error.  The caller is responsible for
 * checking failed allocation.  By default, the base macros MAKE_1ARRAY and
 * CMAKE_1ARRAY are set to the return options.  Redefine wherever you do not
 * wish immediate return, for example when you want to perform extra clean-up
 * before leaving a function that has failed.
 */
#include "array.h"
#define MAKE_1ARRAY_RETURN(a,n) do {                                           \
	(a) = malloc((n) * sizeof *(a));                                       \
	if ((a) == NULL)                                                       \
		return error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__,       \
			__func__, __LINE__, NULL);                             \
} while (0)
#define MAKE_1ARRAY_ERROR(a,n) do {                                            \
	(a) = malloc((n) * sizeof *(a));                                       \
	if ((a) == NULL)                                                       \
		error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__, __func__,    \
			__LINE__, NULL);                                       \
} while (0)
#define CMAKE_1ARRAY_RETURN(a,n) do {                                          \
	(a) = calloc((n), sizeof *(a));                                        \
	if ((a) == NULL)                                                       \
		return error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__,       \
			__func__, __LINE__, NULL);                             \
} while (0)
#define CMAKE_1ARRAY_ERROR(a,n) do {                                           \
	(a) = calloc((n), sizeof *(a));                                        \
	if ((a) == NULL)                                                       \
		error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__, __func__,    \
			__LINE__, NULL);                                       \
} while (0)

/* avoid gcc warnings by undef'ing */
#undef MAKE_1ARRAY
#undef CMAKE_1ARRAY
/* before reseting to rewritten macros */
#define MAKE_1ARRAY MAKE_1ARRAY_ERROR
#define CMAKE_1ARRAY CMAKE_1ARRAY_ERROR

/**
 * Default maximum number of K-means iterations.
 * The user can specify any number of maximum iterations when calling 
 * #kmeans(), but for consistency across the code, this define can be 
 * used.
 */
#define MAX_KMEANS_ITERATIONS	1000

/**
 * Default buffer size for buffered reads.
 */
#define DEFAULT_BUFFER_SIZE	1000

/**
 * Verbosity levels.
 */
enum {
	SILENT,		/*!< no output to screen */
	TERSE,		/*!< minimal output */
	VERBOSE,	/*!< most verbose */
	NUM_VERBOSITY_LEVELS
};

/**
 * Errors produced by kmeans().
 */
enum {
	KMEANS_NO_ERROR,		/*!< error-free condition */
	KMEANS_NULL_CLUSTER_ERROR,	/*!< K-means produced null cluster */
	KMEANS_CALLER_INPUT_ERROR,	/*!< invalid caller input */
	KMEANS_UNUSED_ERROR,		/*!< not used */
	KMEANS_NO_SEEDS,		/*!< initialization gives 0 seeds */
	KMEANS_NUMBER_ERRORS,		/*!< number of K-means errors */
	KMEANS_EXCEED_ITER_WARNING,	/*!< K-means exceeded max iterations */
	KMEANS_NUMBER_FAULTS		/*!< total number of fault codes */
};

/** 
 * Human-friendly character strings describing each K-means error.
 */
extern const char *KMEANS_ERROR_STRING[];

typedef struct _options options;
typedef struct _data data;
typedef struct _model model;

/**
 * Options structure.
 */
struct _options {
	/* data i/o */
	const char *path;		/*!< path to files */
	const char *xname;		/*!< data filename */
	const char *idname;		/*!< classification output filename */
	const char *wssname;		/*!< wss output filename */
	const char *muname;		/*!< filename with mu */
	const char *sfile;		/*!< seed filename */
	const char *tfile;		/*!< time filename */
	const char *efile;		/*!< error filename */
	char *full_filename;		/*!< data file name */
	FILE *finp;			/*!< FILE pointer to data file */
	char *mu_full_filename;		/*!< full path mu filename */
	FILE *fmup;			/*!< FILE pointer to true clusters */

	/* run settings */
	SIZE_T K;			/*!< number of clusters [TODO: more] */
	SIZE_T n_datasets;		/*!< number of datasets in data file */
	int conserve_memory;		/*!< computationally intensive, low
					 *   memory version */
	SIZE_T n_max_iter;		/*!< maximum number of iterations */
	double run_seconds;		/*!< no. secs to run stoch. methods */
	int seed_from_file;		/*!< indicate if seed(s) from file */
	unsigned int seed1;		/*!< random number seeds */
	unsigned int seed2;		/*!< random number seeds */
	int verbosity;			/*!< verbosity level */

	/* initialization method */
	int init_method;		/*!< initialization method; see enums */
	const char *method_name;	/*!< human-friendly name of method */

	/* initialization method parameters */
	double as70_subsample_prop;	/*!< subsample size; prop. of data::n */
	double as70_sphere_radius;	/*!< compute density in spheres */
	double as70_separation_radius;	/*!< pick seeds separated by this */
	int ha75a_dim;			/*!< pick single dim for partition */
	int hh93a_dim;			/*!< pick single dim for splitting */
	int hh93a_cluster_choice;	/*!< how to pick cluster to split */
	int bh67_steinley;		/*!< BH67 with Steinley's recommend. */
	double lw67_subsample_prop;	/*!< subsample size; prop. of data::n */
	double use_tg74;		/*!< use TG74 variant of BW67 */
	double bf98_subsample_prop;	/*!< subsample size; prop. of data::n */
	SIZE_T ltty08_dim;		/*!< precision of integer transform */
	double hk05_alpha;		/*!< init. prob. to move observations */
	double hk05_beta;		/*!< prop. decrease in move prob. */
	SIZE_T hk05_Bs;			/*!< stop if Bs iters. w/ no change */
	SIZE_T hk05_Bl;			/*!< stop after Bl iterations */
	double ORSS06b_rho;		/*!< sqrt(epsilon) parameter of ORSS06b
					 * method, measures support for cluster
					 * k over cluster k+1 solution. */
};

/**
 * Data structure.
 */
struct _data {
	SIZE_T n;		/*!< number of observations */
	SIZE_T p;		/*!< number of dimensions */
	double **x;		/*!< data matrix */
};

/**
 * Model structure.  Variables defining model to fit and variables for k-means
 * algorithm.
 */
struct _model {
	SIZE_T K;		/*!< number of clusters */
	SIZE_T real_K;		/*!< realized K */
	int ifault;		/*!< error status of k-means */

	double **means;		/*!< initialization seeds */

	SIZE_T *iclass;		/*!< inferred class assignments */
	SIZE_T *nc;		/*!< inferred class sizes */
	double *w;		/*!< cluster within-sums-of-squares */
	double wsum;		/*!< within-sum-of-squares */
	double time;		/*!< actual time used */
};


void kmeans(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	SIZE_T *nc, SIZE_T iter, double *wss, int *ifault); 

const char *kmeans_error(int err);

#endif /* KMEANS_H */
