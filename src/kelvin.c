#include "ml_kelvin.h"

#define ML_KELVIN_XMAX 20.0
#define ML_KELVIN_KXMAX 12.0
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
    /* DLMF 10.67.3-4 full sums (verified vs mpmath 80-dps):
     * ber/bei: pref_b * sum a_k/x^k {cos,sin}(phb+3k pi/4),
     * ker/kei: pref_k * sum a_k/x^k {cos,-sin}(phk+k pi/4),
     * a_k(0)=prod_{j<=k}(0-(2j-1)^2)/(k! 8^k), least-term truncation.
     * Round-3: replaces 2-term P/Q (0.5% at 20) with 1e-12 at 20. */
    long double ax = (long double)x;
    long double e = __builtin_expl(ax * 0.70710678118654752440L);
    long double f = __builtin_expl(-ax * 0.70710678118654752440L);
    long double sqb = 0.39894228040143267794L / __builtin_sqrtl(ax);
    long double sqk = 1.25331413731550025121L / __builtin_sqrtl(ax);
    long double phb = ax * 0.70710678118654752440L - 0.39269908169872415481L;
    long double phk = ax * 0.70710678118654752440L + 0.39269908169872415481L;
    long double invx = 1.0L / ax;
    /* a_k recurrence: a_0=1, a_k=a_{k-1}*(-(2k-1)^2)/(8k). */
    long double sb_cos = 0.0L, sb_sin = 0.0L, sk_cos = 0.0L, sk_sin = 0.0L;
    long double ak = 1.0L, pw = 1.0L;
    long double prev_b = 1e300L, prev_k = 1e300L;
    for (int k = 0; k <= 12; k++) {
        long double thb = phb + (long double)k * 2.35619449019234492885L;
        long double thk = phk + (long double)k * 0.78539816339744830962L;
        long double tb = ak * pw;
        long double cb = __builtin_cosl(thb) * tb, sb = __builtin_sinl(thb) * tb;
        long double ck = __builtin_cosl(thk) * tb, sk = __builtin_sinl(thk) * tb;
        long double magb = __builtin_fabsl(tb), magk = __builtin_fabsl(tb);
        if (k > 0 && (magb > prev_b || magk > prev_k)) break;
        sb_cos += cb; sb_sin += sb; sk_cos += ck; sk_sin += sk;
        prev_b = magb; prev_k = magk;
        long double ff = (long double)(2 * k + 1);
        ak = ak * (-(ff * ff)) / ((long double)(8 * (k + 1)));
        pw *= invx;
        if (__builtin_fabsl(tb) < 1e-22L) break;
    }
    *ber = (double)(sqb * e * sb_cos);
    *bei = (double)(sqb * e * sb_sin);
    *ker = (double)(sqk * f * sk_cos);
    *kei = (double)(sqk * f * (-sk_sin));
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
    if (x < ML_KELVIN_KXMAX) {
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
    if (x < ML_KELVIN_KXMAX) {
        double k, ki;
        ml_kelvin_kerkei_series(x, &k, &ki);
        return ki;
    }
    double b, bi, k, ki;
    ml_kelvin_asym(x, &b, &bi, &k, &ki);
    return ki;
}
