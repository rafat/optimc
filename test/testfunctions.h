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

void fk10(double *x, int M, int N, double *f);

void fk10jac(double *x, int M, int N, double *f);

void fpowell(double *x,int M,int N,double *f);

double myvalue
(
    double   *x,
    int       n
);

void myvaluegrad(
		double *x,
		int n,
		double *jac
);

double quartic(double *x,int N);

void quarticgrad(double *x,int N,double *jac);

double func4(double *x,int N);

double rosenbrock(double *x,int N);

void rosenbrockgrad(double *x,int N,double *g);

double brown(double *x,int N);

void browngrad(double *x,int N,double *g);

double powell(double *x,int N);

void powellgrad(double *x,int N,double *g);

double beale(double *x,int N);

double froth(double *x,int N);

double humps(double x);

double func1(double *x,int N);

#ifdef __cplusplus
}
#endif

#endif /* TESTFUNCTIONS_H_ */
