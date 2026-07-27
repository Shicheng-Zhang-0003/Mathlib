#ifndef MATHLIB_ML_OPTIMIZATION_H
#define MATHLIB_ML_OPTIMIZATION_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_numerical.h"

typedef double (*ml_opt_func_t)(double x);

ML_API double ml_optimize_golden(ml_opt_func_t f, double a, double b, double tol, int max_iter);
ML_API double ml_optimize_gradient_descent(ml_opt_func_t f, double start, double lr, double tol, int max_iter);

#endif /* MATHLIB_ML_OPTIMIZATION_H */
