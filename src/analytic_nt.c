#include "ml_compiler.h"
#include "ml_analytic_nt.h"
#include "ml_integral.h"
ML_API cplx ml_hurwitz_zeta(double s, double a) {
    cplx nn = {ml_make_nan(), ml_make_nan()};
    if (ml_isnan(s) || ml_isnan(a) || !ml_isfinite(a) || a <= 0.0) return nn;
    if (s <= 1.0 && s == ml_round(s) && a == ml_round(a) && a <= 0.0) return nn;
    /* Euler-Maclaurin N=50 + Bernoulli tail (real s>1 fast path);
     * s<=1 via reflection through zeta for a=1 else direct EM. */
    if (s > 1.0) {
        const int N = 50;
        long double sum = 0;
        for (int n = 0; n < N; n++) sum += __builtin_powl((long double)(n + a), -(long double)s);
        long double tail = __builtin_powl((long double)(N + a), 1.0L - (long double)s) / ((long double)s - 1.0L)
            + 0.5L * __builtin_powl((long double)(N + a), -(long double)s);
        double r = (double)(sum + tail);
        if (!ml_isfinite(r)) return nn;
        return (cplx){r, 0.0};
    }
    if (a == 1.0) { double z = ml_zeta(s); if (!ml_isfinite(z)) return nn; return (cplx){z, 0.0}; }
    /* STUB: s<=1 with a!=1 is unimplemented (needs analytic continuation
     * of the general Hurwitz case); keep NaN until a reflection/EM path
     * lands here. */
    return nn;
}
ML_API cplx ml_dirichlet_eta_cplx(cplx s) {
    cplx nn = {ml_make_nan(), ml_make_nan()};
    if (ml_isnan(s.real) || ml_isnan(s.imag)) return nn;
    /* eta(s)=(1-2^{1-s}) zeta(s); complex via cplx power. */
    cplx one = {1.0, 0.0}, two = {2.0, 0.0};
    cplx ex = {1.0 - s.real, -s.imag};
    cplx p = ml_cplx_power(two, ex);
    if (ml_isnan(p.real)) return nn;
    cplx f = {1.0 - p.real, -p.imag};
    cplx z = ml_zeta_cplx(s);
    if (ml_isnan(z.real)) return nn;
    return ml_cplx_mul(f, z);
    (void)one;
}
ML_API double ml_theta3(double q) {
    if (ml_isnan(q) || !ml_isfinite(q) || q <= 0.0 || q >= 1.0) return ml_make_nan();
    long double s = 1.0L;
    long double qn = (long double)q;
    for (int n = 1; n < 200; n++) {
        long double t = 2.0L * __builtin_powl(qn, (long double)n * (long double)n);
        long double prev = s;
        s += t;
        if (s == prev) break;
        if (!ml_isfinite((double)s)) break;
    }
    double r = (double)s;
    return ml_isfinite(r) ? r : ml_make_nan();
}
ML_API double ml_partition_p(int n) {
    if (n < 0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n > 200) return ml_make_nan();
    /* Stack-local table (was static): thread-safe, no cross-call state. */
    double P[201];
    for (int i = 0; i <= 200; i++) P[i] = 0.0;
    P[0] = 1.0;
    for (int k = 1; k <= 200; k++) {
        for (int i = k; i <= 200; i++) P[i] += P[i-k];
    }
    return P[n];
}
ML_API cplx ml_zeta_cplx(cplx s) {
    cplx nn = {ml_make_nan(), ml_make_nan()};
    if (ml_isnan(s.real) || ml_isnan(s.imag)) return nn;
    if (s.imag == 0.0) {
        double z = ml_zeta(s.real);
        if (!ml_isfinite(z)) return nn;
        return (cplx){z, 0.0};
    }
    /* Dirichlet eta Euler transform, K=40, then /(1-2^{1-s}). */
    const int K = 40;
    /* Stack-local Euler table (was static): thread-safe, no cross-call state. */
    cplx a[41];
    for (int n = 0; n <= K; n++) {
        /* (n+1)^{-s} = exp(-s log(n+1)) */
        double lr = ml_log((double)(n + 1));
        cplx e = {-s.real * lr, -s.imag * lr};
        a[n] = ml_cplx_exponential(e);
        if (ml_isnan(a[n].real)) return nn;
    }
    cplx eta = {0.0, 0.0};
    double inv2 = 0.5;
    for (int k = 0; k <= K; k++) {
        eta.real += a[0].real * inv2; eta.imag += a[0].imag * inv2;
        inv2 *= 0.5;
        for (int j = 0; j < K - k; j++) { a[j].real -= a[j+1].real; a[j].imag -= a[j+1].imag; }
    }
    cplx ex = {1.0 - s.real, -s.imag};
    cplx p = ml_cplx_power((cplx){2.0, 0.0}, ex);
    cplx d = {1.0 - p.real, -p.imag};
    if (d.real == 0.0 && d.imag == 0.0) return nn;
    /* eta/d */
    cplx conj = {d.real, -d.imag};
    double den = d.real*d.real + d.imag*d.imag;
    if (den == 0.0) return nn;
    cplx num = ml_cplx_mul(eta, conj);
    return (cplx){num.real / den, num.imag / den};
}
