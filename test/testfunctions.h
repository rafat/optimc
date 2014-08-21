/*
 * testfunctions.h
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#ifndef TESTFUNCTIONS_H_
#define TESTFUNCTIONS_H_


#ifdef __cplusplus
extern "C" {
#endif

	void fk10(double *x, int M, int N, double *f, void *params);

	void fk10jac(double *x, int M, int N, double *f, void *params);

	void fpowell(double *x, int M, int N, double *f, void *params);

double myvalue
(
    double   *x,
	int       n, 
	void *params
);

void myvaluegrad(
		double *x,
		int n,
		double *jac,
		void *params
);

double quartic(double *x,int N,void *params);

void quarticgrad(double *x, int N, double *jac, void *params);

double func4(double *x, int N, void *params);

double rosenbrock(double *x,int N,void *params);

void rosenbrockgrad(double *x, int N, double *g, void *params);

double brown(double *x, int N, void *params);

void browngrad(double *x, int N, double *g, void *params);

double powell(double *x, int N, void *params);

void powellgrad(double *x, int N, double *g, void *params);

double beale(double *x, int N, void *params);

double froth(double *x, int N, void *params);

double humps(double x,void *params);

double func1(double *x, int N, void *params);

double tf6(double *x,int N,void *params);

#ifdef __cplusplus
}
#endif

#endif /* TESTFUNCTIONS_H_ */
