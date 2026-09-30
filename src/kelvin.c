#include "ml_kelvin.h"

#define ML_KELVIN_XMAX 20.0
#define ML_KELVIN_GAMMA 0.57721566490153286061

static double ml_kelvin_ber_series(double x) {
    long double y = (long double)x * 0.5L;
    long double y4 = y * y * y * y;
    long double t = 1.0L, s = 1.0L, c = 0.0L;
    int k;
    for (k = 1; k < 300; k++) {
        long double w, tt, d;
        d = (long double)(2 * k - 1) * (long double)(2 * k);
        t *= -y4 / (d * d);
        w = t - c;
        tt = s + w;
        c = (tt - s) - w;
        s = tt;
        if (__builtin_fabsl(t) < 1e-22L * __builtin_fabsl(s)) break;
    }
    return (double)s;
}

static double ml_kelvin_bei_series(double x) {
    long double y = (long double)x * 0.5L;
    long double y4 = y * y * y * y;
    long double t = y * y, s = y * y, c = 0.0L;
    int k;
    for (k = 1; k < 300; k++) {
        long double w, tt, d;
        d = (long double)(2 * k) * (long double)(2 * k + 1);
        t *= -y4 / (d * d);
        w = t - c;
        tt = s + w;
        c = (tt - s) - w;
        s = tt;
        if (__builtin_fabsl(t) < 1e-22L * __builtin_fabsl(s)) break;
    }
    return (double)s;
}

static void ml_kelvin_kerkei_series(double x, double *ker, double *kei) {
    long double y = (long double)x * 0.5L;
    long double y2 = y * y;
    long double ber = (long double)ml_kelvin_ber_series(x);
    long double bei = (long double)ml_kelvin_bei_series(x);
    long double log_term = (ml_log(x * 0.5) + ML_KELVIN_GAMMA);
    long double sr = 0.0L, si = 0.0L, cr = 0.0L, ci = 0.0L, hn = 0.0L;
    long double t = 1.0L;
    int k;
    for (k = 1; k < 300; k++) {
        long double wr, wi;
        hn += 1.0L / (long double)k;
        t *= y2 / ((long double)k * (long double)k);
        if (k % 4 == 0)       { wr = t * hn; wi = 0.0L; }
        else if (k % 4 == 1)  { wr = 0.0L; wi = t * hn; }
        else if (k % 4 == 2)  { wr = -t * hn; wi = 0.0L; }
        else                  { wr = 0.0L; wi = -t * hn; }
        {
            long double w = wr - cr, tt = sr + w;
            cr = (tt - sr) - w; sr = tt;
        }
        {
            long double w = wi - ci, tt = si + w;
            ci = (tt - si) - w; si = tt;
        }
        if (t * hn < 1e-22L) break;
    }
    *ker = (double)(-log_term * ber + 0.78539816339744830962L * bei + sr);
    *kei = (double)(-log_term * bei - 0.78539816339744830962L * ber + si);
}

static void ml_kelvin_asym(double x, double *ber, double *bei, double *ker, double *kei) {
    long double ax = (long double)x;
    long double e = __builtin_expl(ax * 0.70710678118654752440L);
    long double f = __builtin_expl(-ax * 0.70710678118654752440L);
    long double sqb = 0.39894228040143267794L / __builtin_sqrtl(ax);
    long double sqk = 1.25331413731550025121L / __builtin_sqrtl(ax);
    long double ph = ax * 0.70710678118654752440L + 0.39269908169872415481L;
    long double c = __builtin_cosl(ph), sn = __builtin_sinl(ph);
    long double invx = 1.0L / ax;
    long double invx2 = invx * invx;
    long double P = 1.0L + (9.0L / 128.0L) * invx2;
    long double Q = (1.0L / 8.0L) * invx + (225.0L / 3072.0L) * invx2 * invx;
    *ber = (double)(sqb * e * (c * P + sn * Q));
    *bei = (double)(sqb * e * (sn * P - c * Q));
    *ker = (double)(sqk * f * (c * P + sn * Q));
    *kei = (double)(sqk * f * (c * Q - sn * P));
}

ML_API double ml_kelvin_ber(double x) {
    if (ml_isnan(x)) return x;
    /* DESPOT-AUDIT: ber grows-oscillates ~e^{x/sqrt2}/sqrt(x), unbounded.
     * Previous return 0.0 at +Inf was wrong (that is ker/kei behavior).
     * No limit exists: return NaN fail-loud. */
    if (ml_isinf(x)) return ml_make_nan();
    if (x == 0.0) return 1.0;
    if (x < 0.0) return ml_kelvin_ber(-x);
    if (x < ML_KELVIN_XMAX) return ml_kelvin_ber_series(x);
    double b, bi, k, ki;
    ml_kelvin_asym(x, &b, &bi, &k, &ki);
    return b;
}

ML_API double ml_kelvin_bei(double x) {
    if (ml_isnan(x)) return x;
    /* DESPOT-AUDIT: same as ber — bei is unbounded oscillatory, NaN at Inf. */
    if (ml_isinf(x)) return ml_make_nan();
    if (x == 0.0) return 0.0;
    if (x < 0.0) return -ml_kelvin_bei(-x);
    if (x < ML_KELVIN_XMAX) return ml_kelvin_bei_series(x);
    double b, bi, k, ki;
    ml_kelvin_asym(x, &b, &bi, &k, &ki);
    return bi;
}

ML_API double ml_kelvin_ker(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) return ml_make_nan();
    if (ml_isinf(x)) return 0.0;
    if (x < ML_KELVIN_XMAX) {
        double k, ki;
        ml_kelvin_kerkei_series(x, &k, &ki);
        return k;
    }
    double b, bi, k, ki;
    ml_kelvin_asym(x, &b, &bi, &k, &ki);
    return k;
}

ML_API double ml_kelvin_kei(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) return ml_make_nan();
    if (ml_isinf(x)) return 0.0;
    if (x < ML_KELVIN_XMAX) {
        double k, ki;
        ml_kelvin_kerkei_series(x, &k, &ki);
        return ki;
    }
    double b, bi, k, ki;
    ml_kelvin_asym(x, &b, &bi, &k, &ki);
    return ki;
}
