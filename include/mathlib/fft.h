#ifndef LIBMATHC_FFT_H
#define LIBMATHC_FFT_H

#include "ml_compiler.h"
#include "ml_complex.h"
#include "ml_core.h"
#ifdef __cplusplus
extern "C" {
#endif

/*
 * v11S FFT contract (R6: is_supported now enforces the same cap):
 * ml_fft_execute() and ml_ifft_execute() support power-of-two lengths
 * up to 2^24 only. Unsupported lengths (including NULL, n <= 0, or
 * n > 2^24) are safe no-ops. Use ml_fft_is_supported() to query support
 * before calling.
 */
ML_API int ml_fft_is_supported(int n);

ML_API void ml_fft_execute(cplx *x, int n);
ML_API void ml_ifft_execute(cplx *x, int n);

#ifdef __cplusplus
}
#endif
#endif /* LIBMATHC_FFT_H */
