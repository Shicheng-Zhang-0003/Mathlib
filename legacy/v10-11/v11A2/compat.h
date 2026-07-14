#ifndef MATHLIB_COMPAT_H
#define MATHLIB_COMPAT_H

#warning "MathLib V1.1: Using deprecated V1.0 API. Migrate to ml_ prefixed functions."

#define sine(x) ml_sin(x)
#define cosine(x) ml_cos(x)
#define tangent(x) ml_tan(x)
#define exponential(x) ml_exp(x)
#define logarithm(x) ml_log(x)
#define power(x, y) ml_pow(x, y)
#define arctangent(x) ml_atan(x)
#define arcsine(x) ml_asin(x)
#define arccosine(x) ml_acos(x)
#define arccotangent(x) ml_acot(x)
#define arctangent2(y, x) ml_atan2(y, x)
#define logarithm_base(x, b) ml_logb(x, b)
#define hyperbolic_sine(x) ml_sinh(x)
#define hyperbolic_cosine(x) ml_cosh(x)
#define hyperbolic_tangent(x) ml_tanh(x)
#define inverse_hyperbolic_sine(x) ml_asinh(x)
#define inverse_hyperbolic_cosine(x) ml_acosh(x)
#define inverse_hyperbolic_tangent(x) ml_atanh(x)

#endif
