/*
ex1_fminsearch : demonstrates usage of fminsearch. fminsearch is a simple implementation of
Nelder-Mead method and it requires only the knowledge of function values and an initial
point.

retval = fminsearch(custom_function *funcpt,int N,double *xi,double *xf);

funcpt - function to be minimized. See below and check testfunctions.c for function examples.
N - Problem size. Number of Variables. fminsearch works best for N <= 15
xi - Initial point (N-dimensional)
xf - Minimized point (N-dimensional)
retval - Return Value ( 1 for success. 0 and 4 for failure. 0 - Input error. 4 - Maximum
Iterations exceeded)

*/

#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {
	int N, i, retval;
	double *xi, *xf;

	custom_function quartic_min = {quartic, 0};
	custom_function rosenbrock_min = {rosenbrock, 0};

	N = 4;

	xi = (double*)malloc(sizeof(double)* N);
	xf = (double*)malloc(sizeof(double)* N);

	xi[0] = 3; xi[1] = -1; xi[2] = 0; xi[3] = -1;


	retval = fminsearch(&quartic_min, N, xi, xf);

	printf("Powell's Quartic Function \n");
	printf("Return Value %d \nObjective Function %g \n", retval, quartic(xf, N,0));

	printf("Function minimized at : ");
	for (i = 0; i < N; ++i) {
		printf("%g ", xf[i]);
	}
	printf("\n \n");
	free(xi);
	free(xf);
	printf("\nRosenbrock Function \n");
	N = 2;
	xi = (double*)malloc(sizeof(double)* N);
	xf = (double*)malloc(sizeof(double)* N);

	xi[0] = -1.2; xi[1] = 1;
	retval = fminsearch(&rosenbrock_min, N, xi, xf);

	printf("Return Value %d \nObjective Function %g \n", retval, rosenbrock(xf, N,0));

	printf("Function minimized at : ");
	for (i = 0; i < N; ++i) {
		printf("%g ", xf[i]);
	}


	return 0;
}
