#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {
	int N, j, i, retval;
	double *xi, *xf;
	double maxstep;

	custom_function myvalue_min = { myvalue, 0 };
	custom_gradient myvaluegrad_min = { myvaluegrad, 0 };

	N = 4;

	xi = (double*)malloc(sizeof(double)* N);
	xf = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		xi[i] = 1;
	}
	printf("\n\n%-25s%-20s%-20s \n", "Method", "Return Value", "Objective Function");
	maxstep = 1000.0; // Set the value accordingly 



	for (j = 0; j < 7; ++j) {
		retval = fminunc(&myvalue_min, &myvaluegrad_min, N, xi,maxstep, j, xf);
		printf("\n\n%-25d%-20d%-20g \n", j, retval, myvalue(xf, N,0));
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

	return 0;
}
