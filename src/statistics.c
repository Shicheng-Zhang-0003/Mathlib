#include "ml_compiler.h"
#include "ml_statistics.h"
#include "ml_integral.h"

/* v11S CLOSURE IP-13: statistics invalid-argument hardening */

ML_API double ml_mean(const double *data, int n) {
    if (ML_UNLIKELY(data == NULL || n <= 0)) {
        return ml_make_nan();
    }

    /* Welford incremental mean: avoids overflow of naive sum
     * (e.g. mean([1e308,1e308]) must be 1e308, not Inf). */
    double mean = 0.0;
    for (int i = 0; i < n; i++) {
        double x = data[i];
        if (ML_UNLIKELY(!ml_isfinite(x))) {
            return ml_make_nan();
        }
        mean += (x - mean) / (double)(i + 1);
    }
    return mean;
}

ML_API double ml_variance(const double *data, int n) {
    /* Population variance (divide by n). For the sample variant /(n-1)
     * see ml_variance_s; ml_covariance below is sample (/(n-1)). */
    if (ML_UNLIKELY(data == NULL || n <= 0)) {
        return ml_make_nan();
    }

    /* Welford's online algorithm for numerical stability */
    double mean = 0.0;
    double m2 = 0.0;

    for (int i = 0; i < n; i++) {
        double x = data[i];
        if (ML_UNLIKELY(!ml_isfinite(x))) {
            return ml_make_nan();
        }
        double delta = x - mean;
        mean += delta / (i + 1.0);
        double delta2 = x - mean;
        m2 += delta * delta2;
    }

    double var = m2 / (double)n;

    /*
     * Variance is mathematically non-negative.
     * A tiny negative value can appear only from floating-point rounding.
     */
    if (var < 0.0) {
        var = 0.0;
    }

    return var;
}

ML_API double ml_variance_s(const double *data, int n) {
    /* Sample variance (divide by n-1, Bessel's correction). n<=1 has no
     * degrees of freedom left, so it returns NaN. */
    if (ML_UNLIKELY(data == NULL || n <= 1)) {
        return ml_make_nan();
    }

    /* Welford's online algorithm for numerical stability */
    double mean = 0.0;
    double m2 = 0.0;

    for (int i = 0; i < n; i++) {
        double x = data[i];
        if (ML_UNLIKELY(!ml_isfinite(x))) {
            return ml_make_nan();
        }
        double delta = x - mean;
        mean += delta / (i + 1.0);
        double delta2 = x - mean;
        m2 += delta * delta2;
    }

    double var = m2 / (double)(n - 1);

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

    /* Use lgamma for log(C(n, r)) = lgamma(n+1) - lgamma(r+1) - lgamma(n-r+1) */
    double log_coeff = ml_lgamma((double)n + 1.0) - ml_lgamma((double)r + 1.0) - ml_lgamma((double)(n - r) + 1.0);

    double log_pmf = log_coeff;

    if (k > 0) {
        log_pmf += (double)k * ml_log(p);
    }

    if (n - k > 0) {
        /* ml_log1p(-p) stays accurate for tiny p where 1-p rounds to 1. */
        log_pmf += (double)(n - k) * ml_log1p(-p);
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

    /* Infinite mean is invalid input (x=Inf correctly yields 0 tail below). */
    if (ML_UNLIKELY(!ml_isfinite(mu))) {
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

    /* Centered two-pass OLS: textbook raw sums overflow (sum_x2) and
     * cancel catastrophically when |mean(x)|>>std(x) (e.g. 1e8±1).
     * Pass 1: stable means via Welford. Pass 2: centered covariances. */
    double mx = 0.0, my = 0.0;
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];
        if (ML_UNLIKELY(!ml_isfinite(xi) || !ml_isfinite(yi))) {
            *out_m = ml_make_nan();
            *out_b = ml_make_nan();
            return;
        }
        mx += (xi - mx) / (double)(i + 1);
        my += (yi - my) / (double)(i + 1);
    }

    double sxx = 0.0, sxy = 0.0;
    double cxx = 0.0, cxy = 0.0; /* Kahan compensations */
    for (int i = 0; i < n; i++) {
        double dx = x[i] - mx;
        double dy = y[i] - my;
        {
            double p = dx * dx - cxx;
            double t = sxx + p;
            cxx = (t - sxx) - p;
            sxx = t;
        }
        {
            double p = dx * dy - cxy;
            double t = sxy + p;
            cxy = (t - sxy) - p;
            sxy = t;
        }
    }

    if (ML_UNLIKELY(sxx == 0.0 || !ml_isfinite(sxx) || !ml_isfinite(sxy))) {
        *out_m = ml_make_nan();
        *out_b = ml_make_nan();
        return;
    }

    double m = sxy / sxx;
    double b = my - m * mx;

    if (ML_UNLIKELY(!ml_isfinite(m) || !ml_isfinite(b))) {
        *out_m = ml_make_nan();
        *out_b = ml_make_nan();
        return;
    }

    *out_m = m;
    *out_b = b;
}

ML_API double ml_normal_cdf(double x, double mu, double sigma) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(mu) || ml_isnan(sigma) || sigma <= 0.0)) {
        return ml_make_nan();
    }
    if (ML_UNLIKELY(!ml_isfinite(mu))) return ml_make_nan();
    if (ML_UNLIKELY(ml_isinf(sigma))) return 0.5;
    if (ml_isinf(x)) return (x > 0.0) ? 1.0 : 0.0;
    {
        double z = (x - mu) / (sigma * 1.41421356237309504880);
        double c = ml_erfc(-z);
        return 0.5 * c;
    }
}

ML_API double ml_normal_inv(double p, double mu, double sigma) {
    /* Acklam rational starter + Newton refinements on the CDF. ~1e-12. */
    if (ML_UNLIKELY(ml_isnan(p) || ml_isnan(mu) || ml_isnan(sigma))) return ml_make_nan();
    if (ML_UNLIKELY(!(sigma > 0.0) || !ml_isfinite(sigma))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(mu))) return ml_make_nan();
    if (!(p > 0.0) || !(p < 1.0)) {
        if (p == 0.0) return -ml_make_inf(0);
        if (p == 1.0) return ml_make_inf(0);
        return ml_make_nan();
    }
    {
        double z;
        if (p < 0.02425 && p > 0.0) {
            double q = ml_sqrt(-2.0 * ml_log(p));
            z = (((((-7.784894002430293e-03 * q - 0.3223964580411365) * q - 2.400758277161838) * q - 2.549732539343734) * q + 4.374664141464968) * q + 2.938163982698783) / ((((7.784695709041462e-03 * q + 0.3224671290700398) * q + 2.445134137142996) * q + 3.754408661907416) * q + 1.0);
        } else if (p > 0.97575 && p < 1.0) {
            double q = ml_sqrt(-2.0 * ml_log(1.0 - p));
            z = -(((((-7.784894002430293e-03 * q - 0.3223964580411365) * q - 2.400758277161838) * q - 2.549732539343734) * q + 4.374664141464968) * q + 2.938163982698783) / ((((7.784695709041462e-03 * q + 0.3224671290700398) * q + 2.445134137142996) * q + 3.754408661907416) * q + 1.0);
        } else {
            double q = p - 0.5;
            double r = q * q;
            z = (((((-39.69683028665376 * r + 220.9460984245205) * r - 275.9285104469687) * r + 138.3577518672690) * r - 30.66479806614716) * r + 2.506628277459239) * q / (((((-54.47609879822406 * r + 161.5858368580409) * r - 155.6989798598866) * r + 66.80131188771972) * r - 13.28068155288572) * r + 1.0);
        }
        for (int i = 0; i < 10; i++) {
            double e = ml_normal_cdf(mu + sigma * z, mu, sigma) - p;
            double pdf = ml_exp(-0.5 * z * z) / 2.50662827463100050242;
            if (pdf == 0.0 || !ml_isfinite(pdf)) break;
            {
                double d = e / (sigma * pdf);
                z -= d;
                if (ml_fabs(d) < 1e-15) break;
            }
        }
        return mu + sigma * z;
    }
}

ML_API double ml_exponential_pdf(double x, double lambda) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(lambda) || !(lambda > 0.0) || !ml_isfinite(lambda))) {
        return ml_make_nan();
    }
    if (x < 0.0) return 0.0;
    if (x == 0.0) return lambda;
    if (ml_isinf(x)) return 0.0;
    {
        double e = ml_exp(-lambda * x);
        return lambda * e;
    }
}

ML_API double ml_poisson_pmf(int k, double lambda) {
    if (ML_UNLIKELY(k < 0 || !(lambda > 0.0) || !ml_isfinite(lambda))) {
        return ml_make_nan();
    }
    /* log pmf = k*log(l) - l - lgamma(k+1), via DD-stable lgamma. */
    {
        double lp = (k == 0) ? 0.0 : (double)k * ml_log(lambda);
        double lpmf = lp - lambda - ml_lgamma((double)k + 1.0);
        double p = ml_exp(lpmf);
        if (p > 1.0) p = 1.0;
        if (!ml_isfinite(p) || p < 0.0) return ml_make_nan();
        return p;
    }
}

ML_API double ml_median(const double *data, int n, double *tmp) {
    if (ML_UNLIKELY(!data || n <= 0 || !tmp)) return ml_make_nan();
    for (int i = 0; i < n; i++) {
        if (ML_UNLIKELY(!ml_isfinite(data[i]))) return ml_make_nan();
        tmp[i] = data[i];
    }
    for (int i = 1; i < n; i++) {
        double k = tmp[i];
        int j = i - 1;
        while (j >= 0 && tmp[j] > k) { tmp[j + 1] = tmp[j]; j--; }
        tmp[j + 1] = k;
    }
    if ((n % 2) == 1) return tmp[n / 2];
    return tmp[n / 2 - 1] * 0.5 + tmp[n / 2] * 0.5;
}

ML_API double ml_covariance(const double *x, const double *y, int n) {
    /* Sample covariance (divide by n-1), matching ml_variance_s.
     * Contrast ml_variance, which is population (divide by n).
     * DESPOT-AUDIT: n<=1 has no degrees of freedom -> NaN (was 0.0 for
     * n==1, inconsistent with variance_s). */
    if (ML_UNLIKELY(!x || !y || n <= 1)) return ml_make_nan();
    {
        double mx = 0.0, my = 0.0;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(x[i]) || !ml_isfinite(y[i]))) return ml_make_nan();
            mx += (x[i] - mx) / (double)(i + 1);
            my += (y[i] - my) / (double)(i + 1);
        }
        double s = 0.0, c = 0.0;
        for (int i = 0; i < n; i++) {
            double p = (x[i] - mx) * (y[i] - my) - c;
            double t = s + p;
            c = (t - s) - p;
            s = t;
        }
        return s / (double)(n - 1);
    }
}

ML_API double ml_correlation(const double *x, const double *y, int n) {
    if (ML_UNLIKELY(!x || !y || n <= 0)) return ml_make_nan();
    {
        double mx = ml_mean(x, n), my = ml_mean(y, n);
        if (!ml_isfinite(mx) || !ml_isfinite(my)) return ml_make_nan();
        /* Long-double centered accumulators: same 1e8+-1 hardening as
         * linear_regression plus immunity to sxx*syy overflow
         * (double product Inf -> r=0 instead of ~1). */
        long double sxx = 0.0L, syy = 0.0L, sxy = 0.0L;
        for (int i = 0; i < n; i++) {
            long double dx = (long double)x[i] - (long double)mx;
            long double dy = (long double)y[i] - (long double)my;
            sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
        }
        if (!(sxx > 0.0L) || !(syy > 0.0L)) return ml_make_nan();
        if (!ml_isfinite((double)sxx) && !ml_isfinite((double)(sxx / 1e300L))) {
            /* sxx/syy far beyond double: normalize by scale first. */
            long double sc = sxx > syy ? sxx : syy;
            sxx /= sc; syy /= sc; sxy /= sc;
        }
        {
            long double den = __builtin_sqrtl(sxx * syy);
            if (!(den > 0.0L)) return ml_make_nan();
            long double r = sxy / den;
            if (r > 1.0L) r = 1.0L;
            if (r < -1.0L) r = -1.0L;
            return (double)r;
        }
    }
}

ML_API double ml_hmean(const double *data, int n) {
    /* Reciprocal sum can overflow (inputs ~1e-308 -> 1/x ~1e308 each):
     * kept as-is; IEEE then yields 0/NaN rather than a wrong mean. */
    if (ML_UNLIKELY(!data || n <= 0)) return ml_make_nan();
    {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(data[i]) || !(data[i] > 0.0))) return ml_make_nan();
            s += 1.0 / data[i];
        }
        return (double)n / s;
    }
}

ML_API double ml_gmean(const double *data, int n) {
    /* Log-sum then exp: exp overflows to +Inf only when the true geometric
     * mean itself exceeds ~1e308; kept as-is (no spurious overflow, since
     * the product is never formed directly). */
    if (ML_UNLIKELY(!data || n <= 0)) return ml_make_nan();
    {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(data[i]) || !(data[i] > 0.0))) return ml_make_nan();
            s += ml_log(data[i]);
        }
        return ml_exp(s / (double)n);
    }
}

ML_API double ml_geometric_pmf(int k, double p) {
    if (ML_UNLIKELY(k < 1 || ml_isnan(p) || p <= 0.0 || p > 1.0)) return ml_make_nan();
    if (p == 1.0) return (k == 1) ? 1.0 : 0.0;
    /* Log-space exp((k-1)*log1p(-p)+log(p)): pow(1-p,k-1) underflows to 0
     * for far tails and rounds 1-p to 1 for tiny p; log1p keeps both. */
    {
        double lpmf = (double)(k - 1) * ml_log1p(-p) + ml_log(p);
        double r = ml_exp(lpmf);
        if (!ml_isfinite(r)) return 0.0;
        if (r > 1.0) r = 1.0;
        return r;
    }
}

ML_API double ml_uniform_pdf(double x, double a, double b) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(a) || ml_isnan(b))) return ml_make_nan();
    if (ML_UNLIKELY(!(b > a) || !ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    if (x < a || x > b) return 0.0;
    {
        double w = b - a;
        if (ml_isinf(w)) {
            /* b-a overflowed (e.g. [-1e308,1e308]): halve first. */
            double wh = b * 0.5 - a * 0.5;
            if (!(wh > 0.0) || !ml_isfinite(wh)) return ml_make_nan();
            return 0.5 / wh;
        }
        return 1.0 / w;
    }
}

ML_API double ml_binom_mean(int n, double p) {
    if (ML_UNLIKELY(n < 0 || ml_isnan(p) || p < 0.0 || p > 1.0)) return ml_make_nan();
    return (double)n * p;
}

ML_API double ml_binom_var(int n, double p) {
    if (ML_UNLIKELY(n < 0 || ml_isnan(p) || p < 0.0 || p > 1.0)) return ml_make_nan();
    return (double)n * p * (1.0 - p);
}

ML_API double ml_bayes(double p_b_given_a, double p_a, double p_b) {
    if (ML_UNLIKELY(ml_isnan(p_b_given_a) || ml_isnan(p_a) || ml_isnan(p_b))) return ml_make_nan();
    if (ML_UNLIKELY(p_b == 0.0 || !ml_isfinite(p_b))) return ml_make_nan();
    if (ML_UNLIKELY(p_b_given_a < 0.0 || p_b_given_a > 1.0 || p_a < 0.0 || p_a > 1.0 || p_b <= 0.0 || p_b > 1.0)) {
        return ml_make_nan();
    }
    {
        double r = p_b_given_a * p_a / p_b;
        /* Inconsistent inputs give r>1 (e.g. p(b|a)p(a)>p(b)); clamping
         * hid them. Allow 1e-12 rounding slack, else NaN. */
        if (r < 0.0) {
            if (r > -1e-12) r = 0.0;
            else return ml_make_nan();
        }
        if (r > 1.0) {
            if (r <= 1.0 + 1e-12) r = 1.0;
            else return ml_make_nan();
        }
        return r;
    }
}

/* Undergraduate distributions: exact transforms of gamma_p/q. ULP
 * inherits from gamma kernels (<=2 ULP central, <1 ULP on tested grid
 * after long-double argument scaling). */

ML_API double ml_gamma_cdf(double x, double a, double b) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(a) || ml_isnan(b))) return ml_make_nan();
    if (ML_UNLIKELY(!(a > 0.0) || !(b > 0.0) || !ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    if (x <= 0.0) return 0.0;
    if (ml_isinf(x)) return 1.0;
    {
        long double bx = (long double)b * (long double)x;
        if (bx > 1e308L) return 1.0;
        return ml_gamma_p(a, (double)bx);
    }
}

ML_API double ml_chi2_cdf(double x, int k) {
    if (ML_UNLIKELY(ml_isnan(x) || k <= 0)) return ml_make_nan();
    if (x <= 0.0) return 0.0;
    if (ml_isinf(x)) return 1.0;
    return ml_gamma_p((double)k * 0.5, x * 0.5);
}

ML_API double ml_exponential_cdf(double x, double lambda) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(lambda) || !(lambda > 0.0) || !ml_isfinite(lambda))) return ml_make_nan();
    if (x <= 0.0) return 0.0;
    if (ml_isinf(x)) return 1.0;
    return -ml_expm1(-lambda * x);
}

static double ml_beta_cf(double a, double b, double x) {
    /* Lentz for incomplete beta I_x(a,b), |x|<1. Standard continued
     * fraction (Abramowitz-Stegun 26.5.8). */
    static const double FPMIN = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (ml_fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m < 500; m++) {
        int m2 = 2 * m;
        double aa = (double)m * (b - (double)m) * x / ((qam + (double)m2) * (a + (double)m2));
        d = 1.0 + aa * d;
        if (ml_fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (ml_fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -((a + (double)m) * (qab + (double)m) * x / ((a + (double)m2) * (qap + (double)m2)));
        d = 1.0 + aa * d;
        if (ml_fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (ml_fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        {
            double del = d * c;
            h *= del;
            if (ml_fabs(del - 1.0) < 1e-15) break;
        }
    }
    return h;
}

ML_API double ml_beta_cdf(double x, double a, double b) {
    if (ML_UNLIKELY(ml_isnan(x) || ml_isnan(a) || ml_isnan(b))) return ml_make_nan();
    if (ML_UNLIKELY(!(a > 0.0) || !(b > 0.0) || !ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    {
        /* Front factor x^a(1-x)^b/(a B(a,b)) in log space, then CF.
         * Use symmetry I_x(a,b)=1-I_{1-x}(b,a) for x above the mean. */
        double lb = ml_lgamma(a + b) - ml_lgamma(a) - ml_lgamma(b);
        if (x < (a + 1.0) / (a + b + 2.0)) {
            double lf = a * ml_log(x) + b * ml_log1p(-x) - ml_log(a) + lb;
            double front = ml_exp(lf);
            if (!ml_isfinite(front)) return ml_make_nan();
            return front * ml_beta_cf(a, b, x);
        }
        {
            double y = 1.0 - x;
            double lf = b * ml_log(y) + a * ml_log(x) - ml_log(b) + lb;
            double front = ml_exp(lf);
            if (!ml_isfinite(front)) return ml_make_nan();
            return 1.0 - front * ml_beta_cf(b, a, y);
        }
    }
}

ML_API double ml_student_t_cdf(double t, int nu) {
    if (ML_UNLIKELY(ml_isnan(t) || nu <= 0)) return ml_make_nan();
    if (ml_isinf(t)) return (t > 0.0) ? 1.0 : 0.0;
    if (t == 0.0) return 0.5;
    {
        /* t_nu relation to I_x: x = nu/(nu+t^2). */
        long double tt = (long double)t;
        long double x = (long double)nu / ((long double)nu + tt * tt);
        double ib = ml_beta_cdf((double)x, (double)nu * 0.5, 0.5);
        if (!ml_isfinite(ib)) return ml_make_nan();
        if (t > 0.0) return 1.0 - 0.5 * ib;
        return 0.5 * ib;
    }
}

ML_API double ml_f_cdf(double x, int d1, int d2) {
    if (ML_UNLIKELY(ml_isnan(x) || d1 <= 0 || d2 <= 0)) return ml_make_nan();
    if (x <= 0.0) return 0.0;
    if (ml_isinf(x)) return 1.0;
    {
        long double z = ((long double)d1 * (long double)x) / ((long double)d1 * (long double)x + (long double)d2);
        if (z >= 1.0L) return 1.0;
        if (z <= 0.0L) return 0.0;
        return ml_beta_cdf((double)z, (double)d1 * 0.5, (double)d2 * 0.5);
    }
}
