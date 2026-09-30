#ifndef MATHLIB_ML_EXP_LOG_H
#define MATHLIB_ML_EXP_LOG_H

#include "ml_compiler.h"
#include "ml_core.h"
#ifdef __cplusplus
extern "C" {
#endif

ML_API double ml_exp(double x);
ML_API double ml_log(double x);
ML_API void ml_log_split(double x, double *log_hi, double *log_lo);
ML_API double ml_pow(double x, double y);
ML_API double ml_logb(double x, double b);
ML_API double ml_expm1(double x);
ML_API double ml_log1p(double x);
ML_API double ml_exp2(double x);
ML_API double ml_log2(double x);
ML_API double ml_log10(double x);
ML_API double ml_cbrt(double x);
ML_API double ml_erf(double x);
ML_API double ml_erfc(double x);
ML_API double ml_erfinv(double p);
ML_API double ml_exp10(double x);
ML_API double ml_sech(double x);
ML_API double ml_csch(double x);
ML_API double ml_coth(double x);
ML_API double ml_asech(double x);
ML_API double ml_acsch(double x);
ML_API double ml_acoth(double x);
ML_API double ml_lambert_w0(double x);
ML_API double ml_lambert_wm1(double x);
ML_API double ml_sinh(double x);
ML_API double ml_cosh(double x);
ML_API double ml_tanh(double x);
ML_API double ml_asinh(double x);
ML_API double ml_acosh(double x);
ML_API double ml_atanh(double x);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_EXP_LOG_H */
