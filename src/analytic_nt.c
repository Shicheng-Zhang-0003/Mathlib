#include "ml_compiler.h"
#include "ml_analytic_nt.h"
#include "ml_integral.h"
/* Euler-Maclaurin N=200 + Bernoulli tail through B12, all long double
 * with Kahan summation.
 * Round-4: was N=50 + B2 only (~1e-10: 515k ULP at s=1.1, 1.5M at
 * s=1.5 — the B4 term alone is 9e-10 there). N=200 with B2..B12
 * leaves truncation ~1e-20; LD rounding dominates (~1e-19).
 * BB[1..6] = B2,B4,B6,B8,B10,B12. */
static cplx ml_hurwitz_em_ld(double s, double a) {
    cplx nn = {ml_make_nan(), ml_make_nan()};
    const int N = 200;
    static const long double BB[7] = {0.0L, 1.0L/6.0L, -1.0L/30.0L,
        1.0L/42.0L, -1.0L/30.0L, 5.0L/66.0L, -691.0L/2730.0L};
    long double sum = 0.0L, comp = 0.0L;
    for (int n = 0; n < N; n++) {
        /* (long double)n + (long double)a (was (long double)(n+a):
         * double addition first rounds 203.7 (57 significant bits) by
         * 0.5 ULP = 1.4e-14, amplified by the pole term into 1.1e-14
         * absolute at (0.01,3.7). LD holds 64 bits: exact here. */
        long double t = __builtin_powl((long double)n + (long double)a, -(long double)s);
        long double w = t - comp;
        long double tt = sum + w;
        comp = (tt - sum) - w;
        sum = tt;
    }
    /* B2..B12 tail: sum B_{2k}(s)_{2k-1}/(2k)! Na^{-s-2k+1}. */
    {
        long double Na = (long double)N + (long double)a;
        long double tail = __builtin_powl(Na, 1.0L - (long double)s) / ((long double)s - 1.0L)
            + 0.5L * __builtin_powl(Na, -(long double)s);
        long double pw = __builtin_powl(Na, -(long double)s - 1.0L);
        for (int k = 1; k <= 6; k++) {
            long double rising = 1.0L;
            for (int j = 0; j < 2 * k - 1; j++)
                rising *= (long double)s + (long double)j;
            long double fact = 1.0L;
            for (int j = 1; j <= 2 * k; j++) fact *= (long double)j;
            tail += BB[k] * rising / fact * pw;
            pw /= Na * Na;
        }
        {
            long double r = sum + tail;
            double o = (double)r;
            if (!ml_isfinite(o)) return nn;
            return (cplx){o, 0.0};
        }
    }
}
ML_API cplx ml_hurwitz_zeta(double s, double a) {
    cplx nn = {ml_make_nan(), ml_make_nan()};
    if (ml_isnan(s) || ml_isnan(a) || !ml_isfinite(a) || a <= 0.0) return nn;
    if (s <= 1.0 && s == ml_round(s) && a == ml_round(a) && a <= 0.0) return nn;
    /* Real s>1 fast path: shared LD Euler-Maclaurin helper above. */
    if (s > 1.0) {
        return ml_hurwitz_em_ld(s, a);
    }
    if (a == 1.0) { double z = ml_zeta(s); if (!ml_isfinite(z)) return nn; return (cplx){z, 0.0}; }
    /* Round-4: 0<s<=1, a!=1 via the same LD Euler-Maclaurin as above
     * (was a STUB returning NaN for valid inputs like hurwitz(0.9,2)).
     * The pole term is the analytic continuation; pole subtraction
     * cancels at most ~130x (s=0.01), harmless in LD. s<=0, a!=1 stays
     * a fail-loud NaN stub (needs the functional equation; documented). */
    if (s > 0.0) {
        return ml_hurwitz_em_ld(s, a);
    }
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
