#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {

	int M, N, i;

	nls_object obj;
	custom_funcmult fk10_min = { fk10, 0 };
	custom_jacobian fk10jac_min = { fk10jac, 0 };

	//Matlab Example http://www.mathworks.in/help/optim/ug/lsqnonlin.html See Examples

	M = 10; N = 2;
	obj = nls_init(M, N);// M >= N
	double xi[2] = { 0.3, 0.4 };
	double diag[2] = { 1.0, 1.0 };

	// Using Custom Scaling nls_scale function that allows diag values to be set by the user.
	// The program will only accept strictly positive diagonal values. diag has the size N

	nls_scale(obj, &fk10_min, &fk10jac_min, diag, xi);
	for (i = 0; i < obj->N; ++i) {
		printf("%g ", obj->xopt[i]);
	}
	nlssummary(obj);


	free_nls(obj);

	// Using Internal Scaling. Function nls.

	obj = nls_init(M, N);

	nls(obj, &fk10_min, &fk10jac_min, xi);
	for (i = 0; i < obj->N; ++i) {
		printf("%g ", obj->xopt[i]);
	}
	nlssummary(obj);

	free_nls(obj);

	return 0;
}
