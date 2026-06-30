#ifndef LIBMATHC_FFT_H
#define LIBMATHC_FFT_H

#include "complex.h"
#include "ml_core.h"

static inline int is_power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

static inline void fft_execute(cplx *x, int n) {
    if (!is_power_of_two(n)) return;

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
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

    // Cooley-Tukey iterative FFT
    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * math_pi / len;
        cplx wlen = {cosine(ang), sine(ang)};
        for (int i = 0; i < n; i += len) {
            cplx w = {1.0, 0.0};
            for (int j = 0; j < len / 2; j++) {
                cplx u = x[i + j];
                cplx v = cplx_mul(x[i + j + len / 2], w);
                x[i + j] = cplx_add(u, v);
                x[i + j + len / 2] = cplx_sub(u, v);
                w = cplx_mul(w, wlen);
            }
        }
    }
}

static inline void ifft_execute(cplx *x, int n) {
    if (!is_power_of_two(n)) return;

    for (int i = 0; i < n; i++) {
        x[i].imag = -x[i].imag;
    }

    fft_execute(x, n);

    for (int i = 0; i < n; i++) {
        x[i].imag = -x[i].imag;
        x[i].real /= n;
        x[i].imag /= n;
    }
}
#endif
