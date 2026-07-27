#include "ml_compiler.h"
#include "fft.h"
#include <stddef.h>

/* v11S CLOSURE IP-11: FFT ergonomics + safe size_t loops */

static inline int ml_fft_is_power_of_two_internal(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

ML_API int ml_fft_is_supported(int n) {
    return ml_fft_is_power_of_two_internal(n);
}

ML_API void ml_fft_execute(cplx *x, int n) {
    if (ML_UNLIKELY(x == NULL)) return;
    if (!ml_fft_is_power_of_two_internal(n)) return;

    size_t ns = (size_t)n;

    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < ns; i++) {
        size_t bit = ns >> 1;

        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }

        j ^= bit;

        if (i < j) {
            cplx temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    /* Cooley-Tukey iterative radix-2 FFT */
    for (size_t len = 2; len <= ns; len <<= 1) {
        double ang = -2.0 * ML_PI / (double)len;
        cplx wlen = {ml_cos(ang), ml_sin(ang)};
        size_t half = len >> 1;

        for (size_t i = 0; i < ns; i += len) {
            cplx w = {1.0, 0.0};

            for (size_t j = 0; j < half; j++) {
                if ((j & 63) == 0) {
                    double theta = ang * (double)j;
                    w = (cplx){ml_cos(theta), ml_sin(theta)};
                }

                cplx u = x[i + j];
                cplx v = ml_cplx_mul(x[i + j + half], w);

                x[i + j] = ml_cplx_add(u, v);
                x[i + j + half] = ml_cplx_sub(u, v);

                w = ml_cplx_mul(w, wlen);
            }
        }

        /* Prevent unnecessary final shift beyond ns */
        if (len > ns / 2) break;
    }
}

ML_API void ml_ifft_execute(cplx *x, int n) {
    if (ML_UNLIKELY(x == NULL)) return;
    if (!ml_fft_is_power_of_two_internal(n)) return;

    size_t ns = (size_t)n;

    for (size_t i = 0; i < ns; i++) {
        x[i].imag = -x[i].imag;
    }

    ml_fft_execute(x, n);

    double inv = 1.0 / (double)ns;

    for (size_t i = 0; i < ns; i++) {
        x[i].imag = -x[i].imag;
        x[i].real *= inv;
        x[i].imag *= inv;
    }
}
