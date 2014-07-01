#ifndef NEWTONMIN_H_
#define NEWTONMIN_H_

#include "lnsrchmp.h"


#ifdef __cplusplus
extern "C" {
#endif

double signx(double x);

double cholmod(double *A, int N, double *L, double eps,double maxinp);

double modelhess(double *A,int N,double *dx,double eps,double *L);

void linsolve_lower(double *L,int N,double *b,double *x);

int hessian_fd(double (*funcpt)(double *,int),double *x,int N,double *dx,double eps,double *f);

int hessian_fd2(double (*funcpt)(double *,int),double *x,int N,double *dx,double eps,double *f);

void fdjac(void (*funcgrad)(double *,int,double *),double *x,int N,double *jac,double *dx,double eps2,double *J);

void hessian_fdg(void (*funcgrad)(double *,int,double *),double *x,int N,double *jac,double *dx,double eps2,double *H);

int hessian_opt(double (*funcpt)(double *,int),void (*funcgrad)(double *,int,double *),double *x,int N,double *jac,
		double *dx,double eps,double eps2,double *H);

int lnsrch(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double * dx,double maxstep,double stol,double *x);

int lnsrchmod(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,double *jac,double *p,int N,double * dx,double maxstep,
	double eps2,double stol,double *x,double *jacf);
	
int lnsrchcg(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,double *jac,double *p,int N,double * dx,double maxstep,
	double eps2,double stol,double *x,double *jacf);
	
int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval);

int newton_min_func(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,int N,double *dx,double fsval,int MAXITER,
		int *niter,double eps,double gtol,double stol,double *xf);

int trsrch(double (*funcpt)(double *,int),double *xi,double *jac,double *sN,int N,double * dx,double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double eps,double *x);
		
void trstep(double *jac,double *sN,int N,double * dx,double *L,double *hess,double nlen,double *ioval,double eps,
		double *step);

int trupdate(double (*funcpt)(double *,int),double *xi,double *jac,double *step,int N,double * dx,double maxstep,
		int retcode,double *L,double *hess,double stol,int method,double *ioval,double *xprev,double *funcprev,double *x);		

void trstep_ddl(double *jac,double *sN,int N,double * dx,double maxstep,double *L,double *hess,double nlen,double *ioval,
		double *ssd,double *v,double *step);		
		
int trsrch_ddl(double (*funcpt)(double *,int),double *xi,double *jac,double *sN,int N,double * dx,double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double *x);		

int newton_min_trust(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),double *xi,int N,double *dx,double fsval,double delta,
		int method,int MAXITER,int *niter,double eps,double gtol,double stol,double *xf);

#ifdef __cplusplus
}
#endif

#endif /* NEWTONMIN_H_ */
