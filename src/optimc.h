/*
 * optimc.h
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#ifndef OPTIMC_H_
#define OPTIMC_H_

#include "nls.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct opt_set* opt_object;

opt_object opt_init(int N);

struct opt_set{
	int N;
	double objfunc;
	double eps;
	double gtol;
	double stol;
	double ftol;
	double xtol;
	int MaxIter;
	int Iter;
	int Method;
	int retval;
	char MethodName[50];
	double xopt[1];
};

typedef struct nls_set* nls_object;

nls_object nls_init(int M,int N);

struct nls_set{
	int M;
	int N;
	double eps;
	double epsfcn;
	double factor;
	double gtol;
	double ftol;
	double xtol;
	int MaxIter;
	int Maxfev;
	int Iter;
	int nfev;
	int njev;
	int ldfjac;
	int mode;
	int retval;
	double xopt[1];
};

void setnlsTOL(nls_object obj,double gtol,double ftol,double xtol);

void summary(opt_object obj);

void setMaxIter(opt_object obj,int MaxIter);

void setTOL(opt_object obj,double gtol,double stol,double ftol,double xtol);

int fminsearch(double (*funcpt)(double *,int),int N,double *xi,double *xf);

double fminbnd(double (*funcpt)(double),double a, double b);

int fminunc(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,int method,double *xf);

int fminnewt(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
	     double delta,double *dx,double fsval,int method,double *xf);

double brent_local_min(double (*funcpt)(double ),double a, double b, double t, double eps, double *x);

void optimize(opt_object obj,double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		int method);

void free_opt(opt_object object);

int levmar(void (*funcmult)(double *,int,int,double *),void(*jacobian)(double *, int,int,double *),
		double *xi,int M, int N,double *xf);

void nls(nls_object obj,void (*funcmult)(double *,int,int,double *),void(*jacobian)(double *, int,int,double *),
		double *xi);

void nls_scale(nls_object obj,void (*funcmult)(double *,int,int,double *),void(*jacobian)(double *, int,int,double *),
		double *diag,double *xi);

void nlssummary(nls_object obj);

void free_nls(nls_object object);

#ifdef __cplusplus
}
#endif

#endif /* OPTIMC_H_ */
