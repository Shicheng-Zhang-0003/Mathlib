#ifndef MATHLIB_ML_HARMONIC_H
#define MATHLIB_ML_HARMONIC_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif
ML_API void ml_haar_fwt(const double *x, double *avg, double *det, int n);
ML_API void ml_haar_iwt(const double *avg, const double *det, double *x, int n);
ML_API double ml_morlet_cwt(const double *x, int n, double dt, double s, int k);
ML_API ml_status_t ml_fft_real(const double *x, double *re, double *im, int n);
ML_API ml_status_t ml_fft2d_pow2(const double *re_in, const double *im_in, double *re_out, double *im_out, int n);
#ifdef __cplusplus
}
#endif
#endif
