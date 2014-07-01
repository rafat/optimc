/*

 */

#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {
	int N,j,i,retval;
	double *xi,*xf;
	opt_object optim;

	N = 2;

	xi = (double*) malloc(sizeof(double) * N);
	xf = (double*) malloc(sizeof(double) * N);

	for (i = 0; i < N;++i) {
		xi[i] = 2;
	}
	xi[0] = 0.5;

	optim = opt_init(N);

	for(j = 0; j < 7;++j) {
		optimize(optim,froth,NULL,N,xi,j);
		summary(optim);
		for (i = 0; i < N;++i) {
			xi[i] = 2;
		}
		xi[0] = 0.5;
	}

	free(xi);
	free(xf);
	free_opt(optim);

	return 0;
}
