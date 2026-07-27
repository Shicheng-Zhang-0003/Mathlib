#include "ml_compiler.h"
#include "ml_statistics.h"

/* v11S CLOSURE IP-13: statistics invalid-argument hardening */

ML_API double ml_mean(const double *data, int n) {
    if (ML_UNLIKELY(data == NULL || n <= 0)) {
        return ml_make_nan();
    }

    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        if (ML_UNLIKELY(!ml_isfinite(data[i]))) {
            return ml_make_nan();
        }

        sum += data[i];
    }

    return sum / (double)n;
}

ML_API double ml_variance(const double *data, int n) {
    if (ML_UNLIKELY(data == NULL || n <= 0)) {
        return ml_make_nan();
    }

    double m = ml_mean(data, n);

    if (ml_isnan(m)) {
        return m;
    }

    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double diff = data[i] - m;

        if (ML_UNLIKELY(ml_isnan(diff))) {
            return ml_make_nan();
        }

        sum += diff * diff;
    }

    double var = sum / (double)n;

    /*
     * Variance is mathematically non-negative.
     * A tiny negative value can appear only from floating-point rounding.
     */
    if (var < 0.0) {
        var = 0.0;
    }

    return var;
}

ML_API double ml_stddev(const double *data, int n) {
    double var = ml_variance(data, n);

    if (var < 0.0) {
        var = 0.0;
    }

    return ml_sqrt(var);
}

ML_API double ml_binomial_pmf(int n, int k, double p) {
    if (ML_UNLIKELY(n < 0 || k < 0 || k > n)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(p) || p < 0.0 || p > 1.0)) {
        return ml_make_nan();
    }

    /*
     * Exact edge cases.
     *
     * These must be handled before taking logarithms, because log(0)
     * would otherwise introduce -Inf into the evaluation.
     */
    if (p == 0.0) {
        return (k == 0) ? 1.0 : 0.0;
    }

    if (p == 1.0) {
        return (k == n) ? 1.0 : 0.0;
    }

    /*
     * Use the smaller tail to reduce work:
     * C(n, k) == C(n, n-k)
     */
    int r = (k < n - k) ? k : (n - k);

    /*
     * v11S closure safety ceiling.
     *
     * The log-space loop is exact in structure, but a gigantic r can
     * become a denial-of-vector. If extremely large binomial support is
     * required later, v12 can add an asymptotic approximation path.
     */
    const int max_terms = 1000000;

    if (ML_UNLIKELY(r > max_terms)) {
        return ml_make_nan();
    }

    double log_coeff = 0.0;

    for (int i = 1; i <= r; i++) {
        log_coeff += ml_log((double)(n - r + i)) - ml_log((double)i);
    }

    double log_pmf = log_coeff;

    if (k > 0) {
        log_pmf += (double)k * ml_log(p);
    }

    if (n - k > 0) {
        log_pmf += (double)(n - k) * ml_log(1.0 - p);
    }

    double pmf = ml_exp(log_pmf);

    /*
     * Rounding can push a theoretically <= 1 probability barely above 1.
     */
    if (pmf > 1.0) {
        pmf = 1.0;
    }

    return pmf;
}

ML_API double ml_normal_pdf(double x, double mu, double sigma) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(mu) || ml_isnan(sigma) || sigma <= 0.0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isinf(sigma))) {
        return 0.0;
    }

    double z = (x - mu) / sigma;
    double exponent = -0.5 * z * z;

    /*
     * If x is far enough from mu that the exponent underflows to -Inf,
     * the PDF is effectively zero. Handle this before denominator logic.
     */
    if (ml_isinf(exponent) && exponent < 0.0) {
        return 0.0;
    }

    double denom = ml_sqrt(2.0 * ML_PI) * sigma;

    if (denom == 0.0) {
        return ml_make_inf(0);
    }

    return ml_exp(exponent) / denom;
}

ML_API void ml_linear_regression(const double *x, const double *y, int n, double *out_m, double *out_b) {
    if (ML_UNLIKELY(!x || !y || !out_m || !out_b || n <= 0)) {
        if (out_m) *out_m = ml_make_nan();
        if (out_b) *out_b = ml_make_nan();
        return;
    }

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xy = 0.0;
    double sum_x2 = 0.0;

    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];

        if (ML_UNLIKELY(!ml_isfinite(xi) || !ml_isfinite(yi))) {
            *out_m = ml_make_nan();
            *out_b = ml_make_nan();
            return;
        }

        sum_x += xi;
        sum_y += yi;
        sum_xy += xi * yi;
        sum_x2 += xi * xi;
    }

    double denom = (double)n * sum_x2 - sum_x * sum_x;

    if (ML_UNLIKELY(denom == 0.0 || !ml_isfinite(denom))) {
        *out_m = ml_make_nan();
        *out_b = ml_make_nan();
        return;
    }

    double m = ((double)n * sum_xy - sum_x * sum_y) / denom;
    double b = (sum_y * sum_x2 - sum_x * sum_xy) / denom;

    if (ML_UNLIKELY(!ml_isfinite(m) || !ml_isfinite(b))) {
        *out_m = ml_make_nan();
        *out_b = ml_make_nan();
        return;
    }

    *out_m = m;
    *out_b = b;
}
