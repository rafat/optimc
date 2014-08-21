

#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {
	int N, j, i, retval;
	double fsval, delta,maxstep;
	double *xi, *xf, *dx;

	custom_function brown_min = { brown, 0 };
	custom_gradient browngrad_min = { browngrad, 0 };

	N = 2;

	xi = (double*)malloc(sizeof(double)* N);
	xf = (double*)malloc(sizeof(double)* N);
	dx = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		xi[i] = 1;
		//dx[i] = 1;
	}
	fsval = 1;// Default Value
	delta = -1.0;//Default Value
	dx[0] = 1.0e05; dx[1] = 1.0e-05;
	maxstep = 1000.0;//Default Value
	printf("\n\n%-25s%-20s%-20s \n", "Method", "Return Value", "Objective Function");



	for (j = 1; j < 4; ++j) {
		retval = fminnewt(&brown_min, &browngrad_min, N, xi, delta, dx, fsval,maxstep, j, xf);
		printf("\n\n%-25d%-20d%-20g \n", j, retval, brown(xf, N,0));
		printf("Function Minimized At : ");
		for (i = 0; i < N; ++i) {
			printf("%g ", xf[i]);
		}

		for (i = 0; i < N; ++i) {
			xi[i] = 1;
		}

	}

	free(xi);
	free(xf);
	free(dx);

	return 0;
}
