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
ML_API double ml_lgamma(double x);
ML_API double ml_digamma(double x);
ML_API double ml_gamma_p(double a, double x);
ML_API double ml_gamma_q(double a, double x);
ML_API double ml_beta(double a, double b);
ML_API double ml_bessel_j0(double x);
ML_API double ml_bessel_j1(double x);
ML_API double ml_bessel_y0(double x);
ML_API double ml_bessel_y1(double x);
ML_API double ml_bessel_i0(double x);
ML_API double ml_bessel_i1(double x);
ML_API double ml_bessel_k0(double x);
ML_API double ml_bessel_k1(double x);
ML_API double ml_airy_ai(double x);
ML_API double ml_carlson_rf(double x, double y, double z);
ML_API double ml_carlson_rc(double x, double y);
ML_API double ml_carlson_rj(double x, double y, double z, double p);
ML_API double ml_carlson_rd(double x, double y, double z);
ML_API double ml_ellip_k(double k);
ML_API double ml_ellip_e(double k);
ML_API double ml_ellip_f(double phi, double m);
ML_API double ml_ellip_e_inc(double phi, double m);
/* Gauss 2F1(a,b;c;z): Taylor series valid for |z|<1 (plus Gauss sum at
 * z=1 when c-a-b>0). NO analytic continuation: |z|>=1 (z!=1) returns NaN
 * by design rather than a wrong value. c must not be a non-positive integer. */
ML_API double ml_hyp2f1(double a, double b, double c, double z);
ML_API double ml_hyp1f1(double a, double b, double z);
ML_API double ml_zeta(double s);

#endif /* MATHLIB_ML_INTEGRAL_H */
