#include "fastICA.h"

int main(void) {
	SIZE_T i, j, vrbse = 0, deflate = 1; 
	unsigned int seed1, seed2;
	FILE *finp, *fout;
	ICAObject results;

	double **data_matrix, **w_in, alpha = 1, tollim = 1e-4;
	SIZE_T nn = 262144, pp = 3, ee = 64;

	MAKE_MATRIX(data_matrix, nn, pp);
	MAKE_MATRIX(w_in, ee, ee);

	if (vrbse)
		printf("Start reading data\n");

	finp = fopen("/home/maitra/C.libs/functions/kmeans/Lena.dat", "r");
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

	fout = fopen("in-W.out", "w");
	for (i = 0; i < ee; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", w_in[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);
	
	if (vrbse)
		printf("Done initializing W\n");

	results = fastICA (data_matrix, w_in, nn, pp, ee, alpha, 1, 1, 
			   logcosh, logcoshprime, 200, tollim, deflate, vrbse);
	
	fout = fopen("ICA-S.out", "w");
	for (i = 0; i < nn; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", results.S[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);

	fout = fopen("ICA-W.out", "w");
	for (i = 0; i < ee; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", results.W[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);
	
	fout = fopen("ICA-A.out", "w");
	for (i = 0; i < ee; i++) {
		for (j = 0; j < pp; j++)
			fprintf(fout, " %e", results.A[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);

	fout = fopen("ICA-K.out", "w");
	for (i = 0; i < pp; i++) {
		for (j = 0; j < ee; j++)
			fprintf(fout, " %e", results.K[i][j]);
		fprintf(fout, "\n");
	}
	fclose(fout);

	fout = fopen("ICA-Xpre.out", "w");
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
	



