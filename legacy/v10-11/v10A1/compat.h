#ifndef MATHLIB_COMPAT_H
#define MATHLIB_COMPAT_H
// Legacy API Shims (Deprecated)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#define sine(x) ml_sin(x)
#define cosine(x) ml_cos(x)
#define exponential(x) ml_exp(x)
#define logarithm(x) ml_log(x)
#pragma GCC diagnostic pop
#endif
