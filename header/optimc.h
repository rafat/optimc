/*
 * optimc.h
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#ifndef OPTIMC_H_
#define OPTIMC_H_

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
	double maxstep;
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

typedef struct custom_function_set custom_function;

struct custom_function_set{
	double(*funcpt) (double *x,int N,void *params);// Function in N variables
	void *params;
};

typedef struct custom_gradient_set custom_gradient;

struct custom_gradient_set{
	void(*funcgrad) (double *x, int N,double *g, void *params);
	void *params;
};

typedef struct custom_funcuni_set custom_funcuni;

struct custom_funcuni_set{
	double(*funcuni) (double x, void *params);// Function in one variable
	void *params;
};

typedef struct custom_funcmult_set custom_funcmult;

struct custom_funcmult_set{
	void(*funcmult) (double *x,int M, int N, double *f, void *params);// M Functions in N variables
	void *params;
};

typedef struct custom_jacobian_set custom_jacobian;

struct custom_jacobian_set{
	void(*jacobian) (double *x, int M, int N, double *jac, void *params);
	void *params;
};

#define FUNCPT_EVAL(F,x,N) (*((F)->funcpt))(x,N,(F)->params)

#define FUNCGRAD_EVAL(F,x,N,g) (*((F)->funcgrad))(x,N,(g),(F)->params)

#define FUNCUNI_EVAL(F,x) (*((F)->funcuni))(x,(F)->params)

#define FUNCMULT_EVAL(F,x,M,N,f) (*((F)->funcmult))(x,M,N,(f),(F)->params)

#define JACOBIAN_EVAL(F,x,M,N,jac) (*((F)->jacobian))(x,M,N,(jac),(F)->params)

void setnlsTOL(nls_object obj,double gtol,double ftol,double xtol);

void summary(opt_object obj);

void setMaxIter(opt_object obj,int MaxIter);

void setMaxStep(opt_object obj, double maxstep);

void setTOL(opt_object obj,double gtol,double stol,double ftol,double xtol);

int fminsearch(custom_function *funcpt,int N,double *xi,double *xf);

double fminbnd(custom_funcuni *funcuni,double a, double b);

int fminunc(custom_function *funcpt,custom_gradient *funcgrad,int N,double *xi,double maxstep, int method,double *xf);

int fminnewt(custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,
	     double delta,double *dx,double fsval,double maxstep,int method,double *xf);

double brent_local_min(custom_funcuni *funcuni,double a, double b, double t, double eps, double *x);

void optimize(opt_object obj, custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,
		int method);

void free_opt(opt_object object);

int levmar(custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *xi,int M, int N,double *xf);

void nls(nls_object obj, custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *xi);

void nls_scale(nls_object obj, custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *diag,double *xi);

void nlssummary(nls_object obj);

void free_nls(nls_object object);

#ifdef __cplusplus
}
#endif

#endif /* OPTIMC_H_ */
