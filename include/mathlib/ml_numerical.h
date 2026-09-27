#ifndef MATHLIB_ML_NUMERICAL_H
#define MATHLIB_ML_NUMERICAL_H

#include "ml_compiler.h"
#include "ml_core.h"

typedef double (*ml_func_t)(double);

ML_API double ml_newton_raphson(ml_func_t f, ml_func_t df, double x0, double epsilon, int max_iter);
ML_API double ml_bisection(ml_func_t f, double a, double b, double epsilon, int max_iter);
ML_API double ml_brent(ml_func_t f, double a, double b, double tol, int max_iter);
ML_API double ml_kepler(double M, double e, double tol);
ML_API double ml_derivative(ml_func_t f, double x, double h);
ML_API double ml_second_derivative(ml_func_t f, double x, double h);
ML_API double ml_integral_simpson(ml_func_t f, double a, double b, int n);
ML_API double ml_integral_adaptive(ml_func_t f, double a, double b, double tol, int max_depth);
ML_API double ml_integral_tanhsinh(ml_func_t f, double a, double b, double tol);
ML_API double ml_integral_trapezoid(ml_func_t f, double a, double b, int n);
ML_API double ml_arith_nth(double a1, double d, long long n);
ML_API double ml_arith_sum(double a1, double d, long long n);
ML_API double ml_geom_nth(double a1, double r, long long n);
ML_API double ml_geom_sum(double a1, double r, long long n);
/* Exactness guards (integer precision, not range caps): the closed forms
 * are returned only where they are exact in double; beyond n=2^32-1
 * (sum_k) / n=2^21-1 (sum_k2) they return NaN. sum_k3 inherits sum_k. */
ML_API double ml_sum_k(long long n);
ML_API double ml_sum_k2(long long n);
ML_API double ml_sum_k3(long long n);
ML_API int ml_nim_win(const uint64_t *piles, int n);
ML_API int ml_majorizes(const double *a, const double *b, int n, double tol);

#endif /* MATHLIB_ML_NUMERICAL_H */
