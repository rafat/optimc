/*
fminbnd - minimizes a single variable function within a bound given by a and b.

 */

#include <stdio.h>
#include <stdlib.h>
#include "../header/optimc.h"
#include "testfunctions.h"

int main(void) {
	double a,b,oup;
	
	a = 0.3;
	b = 1;

	oup = fminbnd(humps,a,b);
	printf("OUP %g \n",oup);

	return 0;
}
