#include "ml_compiler.h"
#include "ml_stats_inv.h"
#include "ml_integral.h"
#include "ml_statistics.h"
ML_API double ml_gamma_inv(double p, double a, double b) {
    if (ml_isnan(p) || ml_isnan(a) || ml_isnan(b)) return ml_make_nan();
    if (!(p > 0.0) || !(p < 1.0)) { if (p == 0.0) return 0.0; if (p == 1.0) return ml_make_inf(0); return ml_make_nan(); }
    if (!(a > 0.0) || !(b > 0.0) || !ml_isfinite(a) || !ml_isfinite(b)) return ml_make_nan();
    double lo = 0.0, hi = a / b;
    if (hi <= 0) hi = 1.0;
    while (ml_gamma_cdf(hi, a, b) < p) { hi *= 2.0; if (!ml_isfinite(hi)) return ml_make_nan(); if (hi > 1e308) break; }
    for (int i = 0; i < 200; i++) {
        double mid = lo * 0.5 + hi * 0.5;
        if (mid == lo || mid == hi) break;
        double c = ml_gamma_cdf(mid, a, b);
        if (!ml_isfinite(c)) return ml_make_nan();
        if (c < p) lo = mid; else hi = mid;
        if (hi - lo <= 1e-14 * (1.0 + hi)) break;
    }
    return lo * 0.5 + hi * 0.5;
}
ML_API double ml_beta_inv(double p, double a, double b) {
    if (ml_isnan(p) || ml_isnan(a) || ml_isnan(b)) return ml_make_nan();
    if (!(p > 0.0) || !(p < 1.0)) { if (p == 0.0) return 0.0; if (p == 1.0) return 1.0; return ml_make_nan(); }
    if (!(a > 0.0) || !(b > 0.0) || !ml_isfinite(a) || !ml_isfinite(b)) return ml_make_nan();
    double lo = 0.0, hi = 1.0;
    for (int i = 0; i < 200; i++) {
        double mid = lo * 0.5 + hi * 0.5;
        if (mid == lo || mid == hi) break;
        double c = ml_beta_cdf(mid, a, b);
        if (!ml_isfinite(c)) return ml_make_nan();
        if (c < p) lo = mid; else hi = mid;
        if (hi - lo < 1e-15) break;
    }
    return lo * 0.5 + hi * 0.5;
}
ML_API double ml_chi2_inv(double p, int k) {
    if (!(p > 0.0) || !(p < 1.0) || k <= 0) { if (p == 0.0) return 0.0; if (p == 1.0) return ml_make_inf(0); return ml_make_nan(); }
    return ml_gamma_inv(p, (double)k * 0.5, 0.5);
}
ML_API double ml_student_t_inv(double p, int nu) {
    /* Starters: symmetric bracket around the median 0 (CDF(0)=0.5), then
     * doubling expansion of the open side until CDF brackets p -- the same
     * scheme as ml_gamma_inv/ml_f_inv. Fixed [-1e6,1e6] cannot hold heavy
     * tails (e.g. Cauchy p=1-1e-12 needs ~3e11). Total CDF budget <= 300
     * evals: <=50 per expansion side, 200 for bisection. */
    if (ml_isnan(p) || nu <= 0) return ml_make_nan();
    if (!(p > 0.0) || !(p < 1.0)) { if (p == 0.0) return -ml_make_inf(0); if (p == 1.0) return ml_make_inf(0); return ml_make_nan(); }
    if (p == 0.5) return 0.0;
    {
        double lo = -1.0, hi = 1.0;
        if (p < 0.5) hi = 0.0; else lo = 0.0;
        if (p < 0.5) {
            for (int i = 0; i < 50; i++) {
                double c = ml_student_t_cdf(lo, nu);
                if (!ml_isfinite(c)) return ml_make_nan();
                if (c <= p) break;
                lo *= 2.0;
                if (!ml_isfinite(lo) || lo < -1e308) return ml_make_nan();
            }
        } else {
            for (int i = 0; i < 50; i++) {
                double c = ml_student_t_cdf(hi, nu);
                if (!ml_isfinite(c)) return ml_make_nan();
                if (c >= p) break;
                hi *= 2.0;
                if (!ml_isfinite(hi) || hi > 1e308) return ml_make_nan();
            }
        }
        for (int i = 0; i < 200; i++) {
            double mid = lo * 0.5 + hi * 0.5;
            if (mid == lo || mid == hi) break;
            double c = ml_student_t_cdf(mid, nu);
            if (!ml_isfinite(c)) return ml_make_nan();
            if (c < p) lo = mid; else hi = mid;
            if (hi - lo < 1e-13 * (1.0 + ml_fabs(mid))) break;
        }
        return lo * 0.5 + hi * 0.5;
    }
}
ML_API double ml_f_inv(double p, int d1, int d2) {
    if (ml_isnan(p) || d1 <= 0 || d2 <= 0) return ml_make_nan();
    if (!(p > 0.0) || !(p < 1.0)) { if (p == 0.0) return 0.0; if (p == 1.0) return ml_make_inf(0); return ml_make_nan(); }
    double lo = 0.0, hi = 1.0;
    while (ml_f_cdf(hi, d1, d2) < p) { hi *= 2.0; if (!ml_isfinite(hi) || hi > 1e308) return ml_make_nan(); }
    for (int i = 0; i < 300; i++) {
        double mid = lo * 0.5 + hi * 0.5;
        if (mid == lo || mid == hi) break;
        double c = ml_f_cdf(mid, d1, d2);
        if (!ml_isfinite(c)) return ml_make_nan();
        if (c < p) lo = mid; else hi = mid;
        if (hi - lo < 1e-13 * (1.0 + hi)) break;
    }
    return lo * 0.5 + hi * 0.5;
}
ML_API double ml_normal_logcdf(double x, double mu, double sigma) {
    if (ml_isnan(x) || ml_isnan(mu) || ml_isnan(sigma) || !(sigma > 0.0) || !ml_isfinite(sigma)) return ml_make_nan();
    if (!ml_isfinite(mu)) return ml_make_nan();
    if (ml_isinf(x)) return (x > 0.0) ? 0.0 : -ml_make_inf(0);
    {
        double z = (x - mu) / sigma;
        if (z >= 0.0) {
            double c = ml_normal_cdf(x, mu, sigma);
            if (c <= 0.0) return -ml_make_inf(0);
            return ml_log(c);
        }
        /* lower tail via erfcx-style: log(0.5*erfc(-z/sqrt2)) with expm1 for accuracy */
        double w = -z / 1.41421356237309504880;
        /* erfc(w) for w>0 large: use gamma_q tail then log */
        double q = ml_gamma_q(0.5, w * w);
        if (q <= 0.0) return -ml_make_inf(0);
        /* 0.5*erfc includes 0.5 factor; erfc(w)=2*Q... use: cdf=0.5*erfc(w) */
        return ml_log(0.5) + ml_log(q > 1.0 ? 1.0 : q);
    }
}
