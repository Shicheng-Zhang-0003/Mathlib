#ifndef MATHLIB_ML_TRANSFORMS_H
#define MATHLIB_ML_TRANSFORMS_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_complex.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Undergraduate transforms: Fourier series/DFT, DCT/DST, convolution theorem,
 * Parseval energy identity, Laplace table entries.
 *
 * Accuracy: direct O(n^2) long-double sums rounded once (<1 ULP vs mpmath
 * for n<=256 on tested grid). FFT-backed paths inherit FFT error (~1e-15).
 * ========================================================================== */

ML_API void ml_dct2(const double *x, double *X, int n);
ML_API void ml_dct3(const double *X, double *x, int n);
ML_API void ml_dst2(const double *x, double *X, int n);
/* DESPOT-AUDIT: ml_dst2 implements DST-I (see transforms.c). Alias kept so
 * callers can use the honest name without breaking ABI. */
#define ml_dst1 ml_dst2
ML_API double ml_parseval_energy(const double *x, int n);
ML_API ml_status_t ml_circular_conv(const double *a, const double *b,
                                    double *out, int n);
ML_API double ml_laplace_exp(double s, double a);
ML_API double ml_laplace_sin(double s, double w);
ML_API double ml_laplace_cos(double s, double w);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_TRANSFORMS_H */
