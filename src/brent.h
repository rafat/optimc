#ifndef BRENT_H_
#define BRENT_H_

#include "secant.h"


#ifdef __cplusplus
extern "C" {
#endif


double brent_zero(double (*funcpt)(double ),double a, double b, double tol, double eps);

double brent_local_min(double (*funcpt)(double ),double a, double b, double t, double eps, double *x);

#ifdef __cplusplus
}
#endif

#endif /* SECANT_H_ */
