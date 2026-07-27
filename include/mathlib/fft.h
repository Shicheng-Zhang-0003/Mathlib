#ifndef LIBMATHC_FFT_H
#define LIBMATHC_FFT_H

#include "ml_compiler.h"
#include "ml_complex.h"
#include "ml_core.h"

/*
 * v11S FFT contract:
 * ml_fft_execute() and ml_ifft_execute() support power-of-two lengths only.
 * Unsupported lengths (including NULL or n <= 0) are safe no-ops.
 * Use ml_fft_is_supported() to query support before calling.
 */
ML_API int ml_fft_is_supported(int n);

ML_API void ml_fft_execute(cplx *x, int n);
ML_API void ml_ifft_execute(cplx *x, int n);

#endif /* LIBMATHC_FFT_H */
