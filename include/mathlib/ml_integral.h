#ifndef MATHLIB_ML_INTEGRAL_H
#define MATHLIB_ML_INTEGRAL_H
#include "ml_core.h"

#include "ml_compiler.h"
#include "ml_exp_log.h"

#ifndef MATHLIB_PI
#define MATHLIB_PI ML_PI
#endif
#ifndef MATHLIB_E
#define MATHLIB_E ML_E
#endif

ML_API double ml_factorial_float(double x);
ML_API double ml_integral_traditional(double a, double b, double exponent, double additive, double d);
ML_API double ml_gamma_new(double x);

#endif /* MATHLIB_ML_INTEGRAL_H */
