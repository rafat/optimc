#ifndef SECANT_H_
#define SECANT_H_

#include "conjgrad.h"


#ifdef __cplusplus
extern "C" {
#endif

void bfgs_naive(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min_naive(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,int N,double *dx,double fsval,int MAXITER,
		double eps, double *xf);

void bfgs_factored(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,int N,double *dx,double fsval,int MAXITER,int *niter,
		double eps,double gtol,double stol,double *xf);

void inithess_l(double *H, int N, int k, double *tsk, double *tyk, double *dx);

void bfgs_rec(double *H, int N, int iter, int m, double *jac, double *sk, double *yk, double *r);

int bfgs_l_min(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,int N,int m,double *dx,double fsval,int MAXITER,int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf);

#ifdef __cplusplus
}
#endif

#endif /* SECANT_H_ */
