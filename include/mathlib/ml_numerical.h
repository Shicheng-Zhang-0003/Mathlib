#ifndef MATHLIB_ML_NUMERICAL_H
#define MATHLIB_ML_NUMERICAL_H

#include "ml_compiler.h"
#include "ml_core.h"

typedef double (*ml_func_t)(double);

ML_API double ml_newton_raphson(ml_func_t f, ml_func_t df, double x0, double epsilon, int max_iter);
ML_API double ml_bisection(ml_func_t f, double a, double b, double epsilon, int max_iter);
ML_API double ml_derivative(ml_func_t f, double x, double h);
ML_API double ml_second_derivative(ml_func_t f, double x, double h);
ML_API double ml_integral_simpson(ml_func_t f, double a, double b, int n);

#endif /* MATHLIB_ML_NUMERICAL_H */
