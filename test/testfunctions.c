/*
 * testfunctions.c
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "testfunctions.h"

//levmar test

void fk10(double *x, int M, int N, double *f, void *params) {
	int i;

	for (i = 1; i <= M; ++i) {
		f[i - 1] = 2 + 2*i - exp(i*x[0]) - exp(i*x[1]);
	}
}

void fk10jac(double *x, int M, int N, double *f, void *params) {
	int i;

	for(i = 1; i <= M;++i) {
		f[(i-1)*N] = -i * exp(i*x[0]);
		f[(i-1)*N+1] = -i * exp(i*x[1]);
	}

}

void fpowell(double *x, int M, int N, double *f, void *params) {
	double pi;

	pi = 3.14159265359;

	if (x[0] > 0.0) {
		f[0] = 10 * (x[2] - 10 * (atan(x[1]/x[0])/(2*pi)));
	} else if (x[0] < 0.0) {
		f[0] = 10 * (x[2] - 10 * (atan(x[1]/x[0])/(2*pi)) - 5);
	}

	f[1] = 10 * (sqrt(x[0]*x[0] + x[1]*x[1]) - 1);
	f[2] = x[2];
}

double myvalue
(
    double   *x,
	int       n, 
	void *params
)
{
    double f ;
    int i ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        //t = i+1 ;
        //t = sqrt (t) ;
        //f += exp (x [i]) - t*x [i] ;
    	f += x[i] * x[i] * x[i] * x[i];
    }
    return (f) ;
}

void myvaluegrad(
		double *x,
		int n,
		double *jac, 
		void *params
)
{
	int i;
	for(i = 0; i < n; ++i) {
		jac[i] = 4 * x[i] * x[i] * x[i];
	}
}

double quartic(double *x,int N,void *params) {
	double f;
	// Powell's Quartic Function
	f = (x[0] + 10 * x[1] ) * (x[0] + 10 * x[1] ) + 5 * (x[2] - x[3]) * (x[2] - x[3])
	+ (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2])
	+ 10 * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) ;
	return f;
}

void quarticgrad(double *x, int N, double *g, void *params) {
	double t1,t2,t3,t4;
	// Gradient of Powell's Quartic Function
	t1 = (x[0] + 10 * x[1]);
	t2 = (x[2] - x[3]);
	t3 = (x[1] - 2 * x[2]);
	t4 = (x[0] - x[3]);

	g[0] = 2*t1 + 40*t4*t4*t4;
	g[1] = 20*t1 + 4*t3*t3*t3;
	g[2] = 10*t2 - 8*t3*t3*t3;
	g[3] = -10*t2 - 40*t4*t4*t4;

}

double rosenbrock(double *x,int N,void *params) {
	double f,alpha,alpha2;
	alpha = 1;
	alpha2 = alpha*alpha;
	f = 100 * (x[0]*x[0]*alpha2 - x[1]/alpha)* (x[0]*x[0]*alpha2 - x[1]/alpha) + (1 - x[0]*alpha) * (1 - x[0]*alpha);

	return f;
}

void rosenbrockgrad(double *x, int N, double *g, void *params) {
	double t,alpha,alpha2;
	alpha = 1;
	alpha2 = alpha * alpha;

	t = (x[0]*x[0]*alpha2 - x[1]/alpha);

	g[0] = 400*t*alpha2*x[0] - 2 * alpha * (1 - x[0]*alpha);
	g[1] = -200*t/alpha;
}

double brown(double *x, int N, void *params) {
	double f;
	f = (x[0] - 1e06)*(x[0] - 1e06) + (x[1] - 2*1e-06)*(x[1] - 2*1e-06) + (x[0]*x[1] - 2)*(x[0]*x[1] - 2);

	return f;
}

void browngrad(double *x, int N, double *g, void *params) {
	double t1,t2,t3;

	t1 = x[0] - 1e06;
	t2 = x[1] - 2*1e-06;
	t3 = x[0]*x[1] - 2;
	//printf("gradient %g %g %g \n",t1,t2,t3);

	g[0] = 2 * t1 + 2 * t3 * x[1];
	g[1] = 2 * t2 + 2 * t3 * x[0];
	//printf("gradient %g %g \n",2 * t2,2 * t3 * x[0]);
}

double powell(double *x, int N, void *params) {
	double f;
	f = (10000*x[0]*x[1] - 1) * (10000*x[0]*x[1] - 1) + (exp(-x[0]) + exp(-x[1]) - 1.0001) * (exp(-x[0]) + exp(-x[1]) - 1.0001);
	return f;
}

void powellgrad(double *x, int N, double *g, void *params) {
	double t1,t2;

	t1 = (10000*x[0]*x[1] - 1);
	t2 = (exp(-x[0]) + exp(-x[1]) - 1.0001);

	g[0] = 2 * t1* 10000 * x[1] - 2 * exp(-x[0]) * t2;
	g[1] = 2 * t1* 10000 * x[0] - 2 * exp(-x[1]) * t2;
}

double beale(double *x, int N, void *params) {
	double f;
	f = (1.5 - x[0] *(1 - x[1])) * (1.5 - x[0] *(1 - x[1])) + (2.25 - x[0] *(1 - x[1] * x[1])) * (2.25 - x[0] *(1 - x[1] * x[1])) +
			 	 (2.625 - x[0] *(1 - x[1] * x[1] * x[1])) * (2.625 - x[0] *(1 - x[1] * x[1] * x[1]));
	return f;
}

double froth(double *x, int N, void *params) {
	double f;
	f = (-13 + x[0] + ((5 - x[1]) * x[1] - 2.0) *x[1]) * (-13 + x[0] + ((5 - x[1]) * x[1] - 2.0) *x[1]) +
			(-29 + x[0] + ((x[1] + 1) * x[1] - 14.0) *x[1]) * (-29 + x[0] + ((x[1] + 1) * x[1] - 14.0) *x[1]);
	return f;
}

double humps(double x,void *params) {
	double f;

	f = 1./((x - 0.3)*(x - 0.3) + 0.01) + 1./((x - 0.9)*(x - 0.9) + 0.04) - 6;

	return f;
}

double func4(double *x, int N, void *params) {
	double f;
	f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);


	return f;
}

double func1(double *x, int N, void *params)
{
	double f;

	//f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);
	//f = 100 * (x[0]*x[0]*alpha2 - x[1]/alpha)* (x[0]*x[0]*alpha2 - x[1]/alpha) + (1 - x[0]*alpha) * (1 - x[0]*alpha);
	//f = (x[0] + 2 * x[1] - 7) * (x[0] + 2 * x[1] - 7) + (2*x[0] + x[1] - 5) * (2*x[0] + x[1] - 5);
/*
	f = (x[0] + 10 * x[1] ) * (x[0] + 10 * x[1] ) + 5 * (x[2] - x[3]) * (x[2] - x[3])
	+ (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2])
	+ 10 * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) ;
*/

	//f = (1.0 /(1.0 + (x[0] - x[1]) * (x[0] - x[1]))) + sin(pi *x[1] * x[2] / 2.0)
	//+ exp(-(((x[0]+x[2])/x[1]) - 2) * (((x[0] +x[2])/x[1]) - 2));
	//f = -1.0 * f;


	//f = (10000*x[0]*x[1] - 1) * (10000*x[0]*x[1] - 1) + (exp(-x[0]) + exp(-x[1]) - 1.0001) * (exp(-x[0]) + exp(-x[1]) - 1.0001);
	//f = log(f);
	//printf("fval %g \n",f);

	//f = (10000*x[0]*x[1] - 1) + (exp(-x[0]) + exp(-x[1]) - 1.0001) ;

	//f = (x[0] - 1e06)*(x[0] - 1e06) + (x[1] - 2*1e-06)*(x[1] - 2*1e-06) + (x[0]*x[1] - 2)*(x[0]*x[1] - 2);
	f = 100 * (x[1]-x[0]*x[0]) * (x[1]-x[0]*x[0]) + ( 1.0 - x[0] ) * ( 1.0 - x[0] )  + 90 *(x[3]-x[2]*x[2])*(x[3]-x[2]*x[2])
	+ ( 1.0 - x[2])*( 1.0 - x[2]) + 10 * (x[1] + x[3] - 2)*(x[1] + x[3] - 2) + 0.1 * (x[1] - x[3]) *(x[1] - x[3]);
	return f;
}

double tf6(double *x,int N,void *params) {
	double f;
	double *p = (double*) params;

	f = p[0] * (x[1]-x[0]*x[0]) * (x[1]-x[0]*x[0]) + ( 1.0 - x[0] ) * ( 1.0 - x[0] )  + p[1] *(x[3]-x[2]*x[2])*(x[3]-x[2]*x[2])
	+ ( 1.0 - x[2])*( 1.0 - x[2]) + 10 * (x[1] + x[3] - 2)*(x[1] + x[3] - 2) + p[2] * (x[1] - x[3]) *(x[1] - x[3]);
	return f;
}
