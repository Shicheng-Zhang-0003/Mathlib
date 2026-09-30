#ifndef MATHLIB_ML_TRIG_H
#define MATHLIB_ML_TRIG_H

#include "ml_compiler.h"
#include "ml_core.h"
#ifdef __cplusplus
extern "C" {
#endif


ML_API double ml_sin(double x);
ML_API double ml_cos(double x);
ML_API double ml_tan(double x);
ML_API double ml_sinpi(double x);
ML_API double ml_cospi(double x);
ML_API double ml_sec(double x);
ML_API double ml_csc(double x);
ML_API double ml_cot(double x);
ML_API double ml_sinc(double x);
ML_API double ml_r_form(double a, double b, double *R, double *alpha);
ML_API double ml_triangle_area_sas(double a, double b, double C);
ML_API double ml_law_cos_side(double a, double b, double C);
ML_API double ml_law_sin_side(double a, double A, double B);
ML_API void ml_polar_to_cart(double r, double theta, double *x, double *y);
ML_API void ml_cart_to_polar(double x, double y, double *r, double *theta);
ML_API double ml_heron(double a, double b, double c);
ML_API double ml_stewart(double a, double b, double c, double m, double n);
ML_API int ml_ceva(double a1, double a2, double b1, double b2, double c1, double c2, double tol);
ML_API double ml_atan(double x);
ML_API double ml_asin(double x);
ML_API double ml_acos(double x);
ML_API double ml_acot(double x);
ML_API double ml_atan2(double y, double x);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_TRIG_H */
