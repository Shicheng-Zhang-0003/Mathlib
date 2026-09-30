#include "ml_compiler.h"
#include "ml_info.h"
#include "ml_exp_log.h"
ML_API double ml_entropy(const double *p, int n) {
    if (!p || n <= 0 || n > 100000) return ml_make_nan();
    long double s = 0;
    for (int i = 0; i < n; i++) {
        if (ml_isnan(p[i]) || p[i] < 0.0 || !ml_isfinite(p[i])) return ml_make_nan();
        if (p[i] == 0.0) continue;
        s -= (long double)p[i] * __builtin_logl((long double)p[i]);
    }
    double r = (double)s;
    return ml_isfinite(r) ? r : ml_make_nan();
}
ML_API double ml_kl_div(const double *p, const double *q, int n) {
    if (!p || !q || n <= 0 || n > 100000) return ml_make_nan();
    long double s = 0;
    for (int i = 0; i < n; i++) {
        if (ml_isnan(p[i]) || ml_isnan(q[i]) || p[i] < 0.0 || q[i] < 0.0) return ml_make_nan();
        if (!ml_isfinite(p[i]) || !ml_isfinite(q[i])) return ml_make_nan();
        if (p[i] == 0.0) continue;
        if (q[i] == 0.0) return ml_make_inf(0);
        s += (long double)p[i] * (__builtin_logl((long double)p[i]) - __builtin_logl((long double)q[i]));
    }
    double r = (double)s;
    return ml_isfinite(r) ? r : ml_make_nan();
}
ML_API double ml_cross_entropy(const double *p, const double *q, int n) {
    if (!p || !q || n <= 0) return ml_make_nan();
    /* DESPOT-AUDIT: H(p)+KL(p||q) with support mismatch is +Inf, not NaN.
     * Previous code mapped any non-finite h/kl to NaN, hiding the legitimate
     * +Inf case (p>0 where q==0). Propagate Inf; NaN only for bad input. */
    double h = ml_entropy(p, n), kl = ml_kl_div(p, q, n);
    if (ml_isnan(h) || ml_isnan(kl)) return ml_make_nan();
    if (!ml_isfinite(h) || !ml_isfinite(kl)) {
        if (ml_isinf(h) || ml_isinf(kl)) {
            /* +Inf dominates: support mismatch. -Inf cannot occur here. */
            return ml_make_inf(0);
        }
        return ml_make_nan();
    }
    return h + kl;
}
ML_API double ml_mi_discrete(const double *joint, int nr, int nc) {
    if (!joint || nr <= 0 || nc <= 0 || nr > 1024 || nc > 1024) return ml_make_nan();
    /* Stack-local marginals (was static): thread-safe, no cross-call state. */
    long double pr[1024], pc[1024];
    for (int i = 0; i < nr; i++) pr[i] = 0;
    for (int j = 0; j < nc; j++) pc[j] = 0;
    for (int i = 0; i < nr; i++) for (int j = 0; j < nc; j++) {
        double v = joint[i*nc+j];
        if (ml_isnan(v) || v < 0.0 || !ml_isfinite(v)) return ml_make_nan();
        pr[i] += v; pc[j] += v;
    }
    long double mi = 0;
    for (int i = 0; i < nr; i++) for (int j = 0; j < nc; j++) {
        long double v = joint[i*nc+j];
        if (v == 0) continue;
        long double d = pr[i] * pc[j];
        if (!(d > 0)) return ml_make_inf(0);
        mi += v * (__builtin_logl(v) - __builtin_logl(d));
    }
    double r = (double)mi;
    return ml_isfinite(r) ? r : ml_make_nan();
}
ML_API double ml_logistic(double x) {
    if (ml_isnan(x)) return x;
    if (x >= 0) { double e = ml_exp(-x); return 1.0 / (1.0 + e); }
    double e = ml_exp(x);
    return e / (1.0 + e);
}
ML_API double ml_softplus(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? x : 0.0;
    if (x > 20.0) return x;
    if (x < -20.0) return ml_exp(x);
    return ml_log1p(ml_exp(x));
}
