#include "ml_compiler.h"
#include "ml_mcmc.h"
static uint64_t ml_mcmc_rng(uint64_t *s) {
    uint64_t x = *s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    *s = x ? x : 0x9E3779B97F4A7C15ULL;
    return *s;
}
static double ml_mcmc_normal(uint64_t *s) {
    double u1 = 0, u2 = 0;
    do { u1 = (double)(ml_mcmc_rng(s) >> 11) * (1.0 / 9007199254740992.0); } while (u1 <= 0.0);
    u2 = (double)(ml_mcmc_rng(s) >> 11) * (1.0 / 9007199254740992.0);
    return __builtin_sqrt(-2.0 * __builtin_log(u1)) * __builtin_cos(6.28318530717958647692 * u2);
}
ML_API ml_status_t ml_mh_sample_ctx(ml_logpdf_t logp, const double *x0, int n, uint64_t seed,
                                    double proposal_std, int steps, double *mean_out, double *acc_rate,
                                    void *ctx) {
    if (!logp || !x0 || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (!(proposal_std > 0.0) || !ml_isfinite(proposal_std)) return ML_ERR_INVALID_ARG;
    if (steps <= 0 || steps > 1000000) return ML_ERR_INVALID_ARG;
    double cur[16], prop[16];
    for (int j = 0; j < n; j++) { if (!ml_isfinite(x0[j])) return ML_ERR_NAN_INPUT; cur[j] = x0[j]; }
    double cur_lp = logp(cur, n, ctx);
    if (!ml_isfinite(cur_lp)) return ML_ERR_NAN_INPUT;
    uint64_t s = seed ? seed : 0x243F6A8885A308D3ULL;
    long double acc = 0;
    long double sum[16];
    for (int j = 0; j < n; j++) sum[j] = 0;
    for (int t = 0; t < steps; t++) {
        for (int j = 0; j < n; j++) prop[j] = cur[j] + proposal_std * ml_mcmc_normal(&s);
        double plp = logp(prop, n, ctx);
        if (!ml_isfinite(plp)) continue;
        double alpha = plp - cur_lp;
        double u = (double)(ml_mcmc_rng(&s) >> 11) * (1.0 / 9007199254740992.0);
        double lu = __builtin_log(u + 1e-300);
        if (lu < alpha) { for (int j = 0; j < n; j++) cur[j] = prop[j]; cur_lp = plp; acc += 1.0L; }
        int burn = steps / 2;
        if (t >= burn) for (int j = 0; j < n; j++) sum[j] += cur[j];
    }
    if (mean_out) for (int j = 0; j < n; j++) mean_out[j] = (double)(sum[j] / (steps - steps / 2));
    if (acc_rate) *acc_rate = (double)(acc / steps);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_mh_sample(ml_logpdf_t logp, const double *x0, int n, uint64_t seed,
                                double proposal_std, int steps, double *mean_out, double *acc_rate) {
    /* Legacy ABI wrapper: no ctx slot, passes NULL ctx to logp.
     * Prefer ml_mh_sample_ctx when the target needs user data. */
    return ml_mh_sample_ctx(logp, x0, n, seed, proposal_std, steps, mean_out, acc_rate, 0);
}
ML_API double ml_kde_gaussian(const double *data, int n, double x, double h) {
    if (!data || n <= 0 || n > 100000) return ml_make_nan();
    if (ml_isnan(x) || ml_isinf(x) || !(h > 0.0) || !ml_isfinite(h)) return ml_make_nan();
    long double s = 0.0L;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(data[i])) return ml_make_nan();
        long double z = ((long double)x - (long double)data[i]) / (long double)h;
        s += __builtin_expl(-0.5L * z * z);
    }
    long double r = s / ((long double)n * (long double)h * 2.50662827463100050242L);
    double o = (double)r;
    return ml_isfinite(o) ? o : ml_make_nan();
}
ML_API double ml_ess(const double *xs, int n) {
    if (!xs || n < 4) return ml_make_nan();
    long double m = 0; for (int i = 0; i < n; i++) { if (!ml_isfinite(xs[i])) return ml_make_nan(); m += xs[i]; }
    m /= n;
    long double v = 0; for (int i = 0; i < n; i++) { long double d = xs[i] - m; v += d * d; }
    if (!(v > 0)) return ml_make_nan();
    v /= n;
    long double rho_sum = 0;
    for (int lag = 1; lag < n / 2; lag++) {
        long double c = 0;
        for (int i = 0; i + lag < n; i++) c += (xs[i] - m) * (xs[i+lag] - m);
        c /= n * v;
        if (c <= 0) break;
        rho_sum += c;
    }
    double r = (double)((long double)n / (1.0L + 2.0L * rho_sum));
    return ml_isfinite(r) ? r : ml_make_nan();
}
ML_API double ml_gelman_rubin(const double *chain1, const double *chain2, int n) {
    if (!chain1 || !chain2 || n < 2) return ml_make_nan();
    long double m1 = 0, m2 = 0;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(chain1[i]) || !ml_isfinite(chain2[i])) return ml_make_nan();
        m1 += chain1[i]; m2 += chain2[i];
    }
    m1 /= n; m2 /= n;
    long double v1 = 0, v2 = 0;
    for (int i = 0; i < n; i++) { long double d1 = chain1[i]-m1, d2 = chain2[i]-m2; v1 += d1*d1; v2 += d2*d2; }
    v1 /= (n - 1); v2 /= (n - 1);
    long double W = (v1 + v2) * 0.5L;
    if (!(W > 0)) return ml_make_nan();
    long double B = (long double)n * (m1 - m2) * (m1 - m2) * 0.5L;
    long double Vh = ((n - 1.0L) / n) * W + B / n;
    double r = (double)__builtin_sqrtl(Vh / W);
    return ml_isfinite(r) ? r : ml_make_nan();
}
