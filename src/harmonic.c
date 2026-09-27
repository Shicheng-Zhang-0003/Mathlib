#include "ml_compiler.h"
#include "ml_harmonic.h"
#include "fft.h"
ML_API void ml_haar_fwt(const double *x, double *avg, double *det, int n) {
    if (!x || !avg || !det || n <= 0 || (n & (n - 1)) != 0) return;
    /* Stack-local scratch (thread-safe; no static storage). */
    double buf[4096];
    if (n > 4096) return;
    for (int i = 0; i < n; i++) { if (!ml_isfinite(x[i])) { avg[0] = ml_make_nan(); return; } buf[i] = x[i]; }
    int m = n;
    int lvl = 0;
    double tmp[4096];
    (void)tmp;
    int off = 0;
    while (m > 1) {
        int h = m / 2;
        for (int i = 0; i < h; i++) {
            double a = (buf[2*i] + buf[2*i+1]) * 0.70710678118654752440;
            double d = (buf[2*i] - buf[2*i+1]) * 0.70710678118654752440;
            buf[i] = a;
            det[off + i] = d;
        }
        /* compact avgs to front */
        for (int i = 0; i < h; i++) tmp[i] = buf[i];
        for (int i = 0; i < h; i++) buf[i] = tmp[i];
        off += h;
        m = h;
        lvl++;
        if (off > 4096) return;
    }
    avg[0] = buf[0];
    (void)lvl;
}
ML_API void ml_haar_iwt(const double *avg, const double *det, double *x, int n) {
    if (!avg || !det || !x || n <= 0 || (n & (n - 1)) != 0) return;
    /* Stack-local scratch (thread-safe; no static storage). */
    double buf[4096];
    if (n > 4096) return;
    buf[0] = avg[0];
    int m = 1, off = 0;
    /* det layout from fwt above is level-major; recompute offsets */
    int total = 0; { int mm = n; while (mm > 1) { total += mm / 2; mm /= 2; } }
    /* invert: walk levels coarse->fine */
    int lens[16]; int nl = 0; { int mm = n; while (mm > 1) { lens[nl++] = mm / 2; mm /= 2; } }
    /* det stored fine-first in our fwt? we stored coarse-last; rebuild by reversing */
    int roff = total;
    for (int li = nl - 1; li >= 0; li--) {
        int h = lens[li];
        roff -= h;
        for (int i = h - 1; i >= 0; i--) {
            double a = buf[i], d = det[roff + i];
            buf[2*i] = (a + d) * 0.70710678118654752440;
            buf[2*i+1] = (a - d) * 0.70710678118654752440;
        }
        m = h * 2;
        (void)m;
    }
    for (int i = 0; i < n; i++) x[i] = buf[i];
    (void)off;
}
ML_API double ml_morlet_cwt(const double *x, int n, double dt, double s, int k) {
    if (!x || n <= 0 || n > 8192) return ml_make_nan();
    if (!(dt > 0.0) || !(s > 0.0) || !ml_isfinite(dt) || !ml_isfinite(s)) return ml_make_nan();
    if (k < 0 || k >= n) return ml_make_nan();
    long double ssum = 0.0L;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(x[i])) return ml_make_nan();
        long double tt = ((long double)i - (long double)k) * (long double)dt / (long double)s;
        long double env = __builtin_expl(-0.5L * tt * tt) / __builtin_sqrtl((long double)s);
        long double w = env * __builtin_cosl(5.0L * tt);
        ssum += (long double)x[i] * w;
    }
    double r = (double)(ssum * (long double)dt * 0.7511255444649425L);
    return ml_isfinite(r) ? r : ml_make_nan();
}
ML_API ml_status_t ml_fft_real(const double *x, double *re, double *im, int n) {
    if (!x || !re || !im || n <= 0 || (n & (n - 1)) != 0 || n > 4096) return ML_ERR_INVALID_ARG;
    /* Stack-local scratch (thread-safe; no static storage). */
    cplx buf[4096];
    for (int i = 0; i < n; i++) { if (!ml_isfinite(x[i])) return ML_ERR_NAN_INPUT; buf[i].real = x[i]; buf[i].imag = 0.0; }
    ml_fft_execute(buf, n);
    for (int i = 0; i < n; i++) { re[i] = buf[i].real; im[i] = buf[i].imag; }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_fft2d_pow2(const double *re_in, const double *im_in, double *re_out, double *im_out, int n) {
    if (!re_in || !im_in || !re_out || !im_out || n <= 0 || (n & (n - 1)) != 0 || n > 128) return ML_ERR_INVALID_ARG;
    /* Stack-local scratch (thread-safe; no static storage). */
    cplx rows[128][128], cols[128];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        if (!ml_isfinite(re_in[i*n+j]) || !ml_isfinite(im_in[i*n+j])) return ML_ERR_NAN_INPUT;
        rows[i][j].real = re_in[i*n+j]; rows[i][j].imag = im_in[i*n+j];
    }
    for (int i = 0; i < n; i++) ml_fft_execute(rows[i], n);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) cols[i] = rows[i][j];
        ml_fft_execute(cols, n);
        for (int i = 0; i < n; i++) rows[i][j] = cols[i];
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { re_out[i*n+j] = rows[i][j].real; im_out[i*n+j] = rows[i][j].imag; }
    return ML_SUCCESS;
}
