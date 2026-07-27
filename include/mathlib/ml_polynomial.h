#ifndef MATHLIB_ML_POLYNOMIAL_H
#define MATHLIB_ML_POLYNOMIAL_H

#include "ml_compiler.h"
#include "ml_core.h"

/* FIX: Added const to non-mutating coefficient arrays */
ML_API double ml_polynomial_eval(const double *coeffs, int degree, double x);
ML_API void ml_polynomial_derivative(const double *coeffs, int degree, double *out);
ML_API double ml_polynomial_newton(const double *coeffs, int degree, double x0, double epsilon, int max_iter);

#endif /* MATHLIB_ML_POLYNOMIAL_H */
