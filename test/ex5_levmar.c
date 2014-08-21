#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {
	int M, N, i, ret;
	custom_funcmult fpowell_min = { fpowell, 0 };

	M = N = 3;

	double xi[3] = { -1, 0, 0 };
	double xf[3] = { 0, 0, 0 };

	//levmar simply returns the end point and the termination code
	// For more options use the object-oriented approach employed by nls and nls_scale
	ret = levmar(&fpowell_min, NULL, xi, M, N, xf);

	printf("Return Value : %d \n", ret);
	printf("\n");
	printf("Termination At : ");
	for (i = 0; i < N; ++i) {
		printf("%g ", xf[i]);
	}

	return 0;
}
