#include "fastICA.h"
#include<string.h>
#include "util.h"
#include "optlist.h"

/**
 * Remove path from a filename including full path.
 *
 * This function is from the optlist library, available on the
 * <a href="http://michael.dipperstein.com/optlist/">web</a>.
 *
 * @param full_path name of file including full path
 * @return name of file without full path
 */

/**
 * Print usage information.
 *
 * @param err Error message, if any.
 * @param
 */
void
print_usage(int status, const char *cmd)
{
	const char *cmd_name = remove_path(cmd);
	printf("\
NAME\n\									\
\t%s - Run the basic fastICA algorithm.\n\
SYNOPSIS\n\
\t%s -n <int> -p <int> -# <int> -D <dir> -i <file> -v -h\n",
	       cmd_name, cmd_name);
	if (status)
		printf("\nTry option -h for more information\n");
	else {
		printf("\
OPTIONS\n\
\t-n <int>   number of observations\n\
\t-p <int>   number of dimensions\n\
\t-# <int>   desired number of independent components\n");
		printf("\
\t-D <dir>   working directory (default: OUTPUT)\n\
\t-i <file>  file containing the dataset\n\
\t-v verbose output, default is no verbosity\n"
			);
	}
	exit(status);
} /* print_usage */

int main(int argc, char *argv[]) {
	SIZE_T i, j, vrbse = 0, deflate = 1; 
	unsigned int seed1, seed2;
	FILE *finp, *fout;
	ICAObject results;

	double **data_matrix, **w_in, alpha = 1, tollim = 1e-4;
	SIZE_T nn, pp, ee;

	option_t *ol = NULL, *copt;	/* command-line processing */
	const char *path = NULL;	/* path to files */
	const char *filename = NULL;	/* data filename */
	const char *fxn_name = "main";
	char *full_filename;		/* full path data filename */
	char *fname;

	ol = GetOptList(argc, argv, "n:p:#:D:i:v:h?");
	copt = ol;

	while (copt != NULL) {
		if (copt->option == 'n') {
			nn = strtoul(copt->argument, NULL, 0);
			if (!nn) {
				error(INVALID_COMMAND_LINE, 0, __FILE__,
				      fxn_name, __LINE__, "Invalid use of -n");
				print_usage(1, argv[0]);
			}
		} else if (copt->option == 'p') {
			pp = strtoul(copt->argument, NULL, 0);
			if (!pp) {
				error(INVALID_COMMAND_LINE, 0, __FILE__,
				      fxn_name, __LINE__, "Invalid use of -p");
				print_usage(1, argv[0]);
			}
		} else if (copt->option == '#') {
			ee = strtoul(copt->argument, NULL, 0);
			if (!ee) {
				error(INVALID_COMMAND_LINE, 0, __FILE__,
				      fxn_name, __LINE__, "Invalid use of -#");
				print_usage(1, argv[0]);
			}
		} else 
			if (copt->option == 'D') 
				path = copt->argument;
			else 
				if (copt->option == 'i') 
					filename = copt->argument;
				else
					if (copt->option == 'v') 
						vrbse = 1;
					else 
						if (copt->option == 'h')
							print_usage(0, argv[0]);
						else 
						{
							error(INVALID_COMMAND_LINE, 0, __FILE__, fxn_name,
							      __LINE__, "Unrecognized option (-%c)", copt->option);
							print_usage(1, argv[0]);
						}
		copt = copt->next;
	}	
	FreeOptList(ol);

	printf("%d %d %d\n", nn, pp, ee);
	
	if (!path)
		path = "OUTPUT";

	MAKE_MATRIX(data_matrix, nn, pp);
	MAKE_MATRIX(w_in, ee, ee);

	if (vrbse)
		printf("Start reading data\n");

	finp = fopen(filename, "r");
	for (i = 0; i < nn; i++)
		for (j = 0; j < pp; j++)
			fscanf(finp, "%lf ", &data_matrix[i][j]);
	fclose(finp);

	if (vrbse)
		printf("Done reading data\n");


	finp = fopen("random.seed", "r");
	fscanf(finp, "%u %u", &seed1, &seed2);
	fclose(finp);

	if (vrbse)
		printf("Done reading in seeds\n");
		
	set_seed(seed1, seed2);

	for (i = 0; i < ee; i++) 
		for (j = 0; j < ee; j++)
			w_in[i][j] = rnorm(0.0, 1.0); /* initialize input unmixing matrix */

/*	fout = fopen("in-W.out", "w");
	for (i = 0; i < ee; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", w_in[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);
*/	
	if (vrbse)
		printf("Done initializing W\n");

	results = fastICA (data_matrix, w_in, nn, pp, ee, alpha, 1, 1, 
			   logcosh, logcoshprime, 200, tollim, deflate, vrbse);

	full_filename = make_full_filename(path, "ICA-S.out");
	fout = fopen(full_filename, "w");
	for (i = 0; i < nn; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", results.S[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);

	full_filename = make_full_filename(path, "ICA-W.out");
	fout = fopen(full_filename, "w");
	for (i = 0; i < ee; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", results.W[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);
	
	full_filename = make_full_filename(path, "ICA-A.out");
	fout = fopen(full_filename, "w");
	for (i = 0; i < ee; i++) {
		for (j = 0; j < pp; j++)
			fprintf(fout, " %e", results.A[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);

	full_filename = make_full_filename(path, "ICA-K.out");
	fout = fopen(full_filename, "w");
	for (i = 0; i < pp; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", results.K[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);

	full_filename = make_full_filename(path, "ICA-Xpre.out");
	fout = fopen(full_filename, "w");
	for (i = 0; i < nn; i++) {
		for (j = 0; j < pp; j++)
			fprintf(fout, " %e", results.Xpre[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);



	FREE_VECTOR(data_matrix);
	FREE_VECTOR(w_in);
	FREE_VECTOR(results.Xpre);
	FREE_VECTOR(results.K);
	FREE_VECTOR(results.W);
	FREE_VECTOR(results.A);
	FREE_VECTOR(results.S);

	get_seed(&seed1,&seed2);

	fout=fopen("random.seed","w");
	fprintf(fout,"%d %d\n",seed1,seed2);
	fclose(fout); 

	return EXIT_SUCCESS;
}
	



