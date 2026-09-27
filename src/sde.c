#include "ml_compiler.h"
#include "ml_sde.h"
static uint64_t ml_sde_rng(uint64_t *s) {
    uint64_t x = *s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    *s = x ? x : 0x9E3779B97F4A7C15ULL;
    return *s;
}
static double ml_sde_normal(uint64_t *s) {
    /* Box-Muller with splitmix stream; deterministic in seed. */
    double u1 = 0, u2 = 0;
    do { u1 = (double)(ml_sde_rng(s) >> 11) * (1.0 / 9007199254740992.0); } while (u1 <= 0.0);
    u2 = (double)(ml_sde_rng(s) >> 11) * (1.0 / 9007199254740992.0);
    return __builtin_sqrt(-2.0 * __builtin_log(u1)) * __builtin_cos(6.28318530717958647692 * u2);
}
ML_API ml_status_t ml_sde_euler_maruyama(ml_sde_drift_t a, ml_sde_diff_t b, double t0, double x0,
                                         double t1, int steps, uint64_t seed, double *x1) {
    if (!a || !b || !x1 || steps <= 0 || steps > 1000000) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(t0) || !ml_isfinite(t1) || !ml_isfinite(x0)) return ML_ERR_NAN_INPUT;
    double dt = (t1 - t0) / steps;
    if (!ml_isfinite(dt) || dt == 0.0) return ML_ERR_INVALID_ARG;
    uint64_t s = seed ? seed : 0x123456789ABCDEFULL;
    double t = t0, x = x0;
    for (int i = 0; i < steps; i++) {
        double ai = a(t, x, 0), bi = b(t, x, 0);
        if (!ml_isfinite(ai) || !ml_isfinite(bi)) return ML_ERR_SINGULAR;
        double dW = __builtin_sqrt(ml_fabs(dt)) * ml_sde_normal(&s);
        x += ai * dt + bi * dW;
        t += dt;
        if (!ml_isfinite(x) || !ml_isfinite(t)) return ML_ERR_SINGULAR;
    }
    *x1 = x;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_sde_milstein(ml_sde_drift_t a, ml_sde_diff_t b,
                                   double (*dbdx)(double t, double x, void *ctx),
                                   double t0, double x0, double t1, int steps, uint64_t seed, double *x1) {
    if (!a || !b || !dbdx || !x1 || steps <= 0 || steps > 1000000) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(t0) || !ml_isfinite(t1) || !ml_isfinite(x0)) return ML_ERR_NAN_INPUT;
    double dt = (t1 - t0) / steps;
    if (!ml_isfinite(dt) || dt == 0.0) return ML_ERR_INVALID_ARG;
    uint64_t s = seed ? seed : 0x123456789ABCDEFULL;
    double t = t0, x = x0;
    for (int i = 0; i < steps; i++) {
        double ai = a(t, x, 0), bi = b(t, x, 0), di = dbdx(t, x, 0);
        if (!ml_isfinite(ai) || !ml_isfinite(bi) || !ml_isfinite(di)) return ML_ERR_SINGULAR;
        double dW = __builtin_sqrt(ml_fabs(dt)) * ml_sde_normal(&s);
        x += ai * dt + bi * dW + 0.5 * bi * di * (dW * dW - dt);
        t += dt;
        if (!ml_isfinite(x)) return ML_ERR_SINGULAR;
    }
    *x1 = x;
    return ML_SUCCESS;
}
ML_API double ml_brownian_bridge(double t, double T, double a, double b, double wT) {
    if (ml_isnan(t) || ml_isnan(T) || ml_isnan(a) || ml_isnan(b) || ml_isnan(wT)) return ml_make_nan();
    if (!(T > 0.0) || t < 0.0 || t > T || !ml_isfinite(T)) return ml_make_nan();
    /* MEAN-ONLY pinned bridge E[W_t | W_0=a, W_T=wT]: a+(wT-a)*t/T.
     * b is the free endpoint mean (kept for API symmetry, must equal wT
     * for a pinned bridge, else linear interpolation). No variance term:
     * for a random draw use ml_brownian_bridge_sample. */
    (void)b;
    return a + (wT - a) * (t / T);
}
ML_API double ml_brownian_bridge_mean(double t, double T, double a, double b, double wT) {
    /* Honest name for the mean-only bridge above. */
    return ml_brownian_bridge(t, T, a, b, wT);
}
ML_API double ml_brownian_bridge_sample(double a, double b, double T, double t, uint64_t seed) {
    /* TRUE Brownian bridge sampler: W_t | W_0=a, W_T=b ~
     * N(a+(b-a)*t/T, t*(T-t)/T). Gaussian via the module's
     * xorshift+Box-Muller stream, deterministic in seed. */
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(T) || ml_isnan(t)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b)) return ml_make_nan();
    if (!(T > 0.0) || !ml_isfinite(T)) return ml_make_nan();
    if (t < 0.0 || t > T || !ml_isfinite(t)) return ml_make_nan();
    if (t == 0.0) return a;
    if (t == T) return b;
    {
        double frac = t / T;
        double mean = a + (b - a) * frac;
        double v = t * (T - t) / T;
        uint64_t s;
        double z, r;
        if (!(v >= 0.0) || !ml_isfinite(v) || !ml_isfinite(mean)) return ml_make_nan();
        if (v == 0.0) return mean;
        s = seed ? seed : 0x123456789ABCDEFULL;
        /* Warm up the xorshift stream: tiny seeds (1,2,...) leave the
         * top bits unmixed on the first outputs, which would skew the
         * first Box-Muller draw. 16 discarded rounds fully mix any seed. */
        for (int w = 0; w < 16; w++) ml_sde_rng(&s);
        z = ml_sde_normal(&s);
        if (!ml_isfinite(z)) return ml_make_nan();
        r = mean + __builtin_sqrt(v) * z;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}
ML_API double ml_ou_exact(double x0, double theta, double mu, double sigma, double t, double dw) {
    if (ml_isnan(x0) || ml_isnan(theta) || ml_isnan(mu) || ml_isnan(sigma) || ml_isnan(t) || ml_isnan(dw)) return ml_make_nan();
    if (!(theta > 0.0) || !(sigma >= 0.0) || t < 0.0 || !ml_isfinite(t)) return ml_make_nan();
    /* X_t = mu+(x0-mu)e^{-th t}+sigma*sqrt((1-e^{-2th t})/2th)*N(0,1); dw=N here */
    long double e1 = __builtin_expl(-(long double)theta * (long double)t);
    long double var = (1.0L - e1 * e1) / (2.0L * (long double)theta);
    if (var < 0) var = 0;
    double r = (double)((long double)mu + ((long double)x0 - (long double)mu) * e1 + (long double)sigma * __builtin_sqrtl(var) * (long double)dw);
    return ml_isfinite(r) ? r : ml_make_nan();
}
