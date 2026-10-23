#include "ml_compiler.h"
#include "ml_stats_extra.h"
#include "ml_statistics.h"
#include "ml_stats_inv.h"
#include "ml_exp_log.h"
#include "ml_integral.h"
#include "ml_trig.h"

/* ============================================================================
 * v12R2 statistics depth pack.
 *
 * Everything here is self-contained (only depends on the public cdf/inverse
 * and gamma/erf machinery) and respects the no-internal-malloc contract:
 * scratch space is supplied by the caller or lives on the stack.
 * ========================================================================== */

#define SE_PI 3.14159265358979323846
/* ------------------------------------------------------------- helpers --- */

static void se_sort (double *a, int n) {
    int i, j, gap;
    for (gap = n / 2; gap > 0; gap /= 2) {
        for (i = gap; i < n; i++) {
            double v = a [i];
            for (j = i; j >= gap && a [j - gap] > v; j -= gap) a [j] = a [j - gap];
            a [j] = v;
        }
    }
}
static double se_log_choose (int n, int k) {
    if (k < 0 || k > n) return -HUGE_VAL;
    return ml_lgamma ((double) n + 1.0) - ml_lgamma ((double) k + 1.0) - ml_lgamma ((double) (n - k) + 1.0);
}
/* Gaussian elimination with partial pivoting; A row-major n x n, b length n.
 * On success writes solution to x. Returns 0 or ML_ERR_SINGULAR. */
static ml_status_t se_solve (double *A, double *b, double *x, int n) {
    int i, j, k;
    for (k = 0; k < n; k++) {
        int piv = k;
        double mx = ml_fabs (A [(size_t) k * (size_t) n + (size_t) k]);
        for (i = k + 1; i < n; i++) {
            double v = ml_fabs (A [(size_t) i * (size_t) n + (size_t) k]);
            if (v > mx) {
                mx = v;
                piv = i;
            }
        }
        if (mx == 0.0) return ML_ERR_SINGULAR;
        if (piv != k) {
            for (j = 0; j < n; j++) {
                double t = A [(size_t) k * (size_t) n + (size_t) j];
                A [(size_t) k * (size_t) n + (size_t) j] = A [(size_t) piv * (size_t) n + (size_t) j];
                A [(size_t) piv * (size_t) n + (size_t) j] = t;
            }
            {
                double t = b [k];
                b [k] = b [piv];
                b [piv] = t;
            }
        }
        {
            double d = A [(size_t) k * (size_t) n + (size_t) k];
            for (i = k + 1; i < n; i++) {
                double f = A [(size_t) i * (size_t) n + (size_t) k] / d;
                b [i] -= f * b [k];
                for (j = k; j < n; j++)
                    A [(size_t) i * (size_t) n + (size_t) j] -=
                        f * A [(size_t) k * (size_t) n + (size_t) j];
            }
        }
    }
    for (i = n - 1; i >= 0; i--) {
        double s = b [i];
        for (j = i + 1; j < n; j++)
            s -= A [(size_t) i * (size_t) n + (size_t) j] * x [j];
        x [i] = s / A [(size_t) i * (size_t) n + (size_t) i];
    }
    return ML_SUCCESS;
}
static double se_twosided_t (double t, int df) {
    return 2.0 * (1.0 - ml_student_t_cdf (ml_fabs (t), df));
}
static double se_twosided_normal (double z) {
    return 2.0 * (1.0 - ml_normal_cdf (ml_fabs (z), 0.0, 1.0));
}
/* ------------------------------------------------------- descriptive ----- */

ML_API double ml_stat_skewness (const double *d, int n) {
    double m, s, s3 = 0.0;
    int i;
    if (ML_UNLIKELY (!d || n < 2)) return ml_make_nan ();
    m = ml_mean (d, n);
    s = ml_stddev (d, n);
    if (!(s > 0.0)) return 0.0;
    for (i = 0; i < n; i++) {
        double z = (d [i] - m) / s;
        s3 += z * z * z;
    }
    return s3 / (double) n;
}
ML_API double ml_stat_kurtosis (const double *d, int n) {
    double m, s, s4 = 0.0;
    int i;
    if (ML_UNLIKELY (!d || n < 2)) return ml_make_nan ();
    m = ml_mean (d, n);
    s = ml_stddev (d, n);
    if (!(s > 0.0)) return 0.0;
    for (i = 0; i < n; i++) {
        double z = (d [i] - m) / s;
        s4 += z * z * z * z;
    }
    return s4 / (double) n - 3.0;
}
ML_API double ml_quantile (const double *d, int n, double q, double *tmp) {
    double pos, frac;
    int lo, hi;
    if (ML_UNLIKELY (!d || !tmp || n <= 0)) return ml_make_nan ();
    if (q <= 0.0) q = 0.0;
    if (q >= 1.0) q = 1.0;
    {
        int i;
        for (i = 0; i < n; i++) tmp [i] = d [i];
    }
    se_sort (tmp, n);
    pos = q * (double) (n - 1);
    lo = (int) pos;
    hi = lo + 1;
    frac = pos - (double) lo;
    if (hi >= n) return tmp [n - 1];
    return tmp [lo] + frac * (tmp [hi] - tmp [lo]);
}
ML_API double ml_percentile (const double *d, int n, double p, double *tmp) {
    return ml_quantile (d, n, p / 100.0, tmp);
}
ML_API double ml_iqr (const double *d, int n, double *tmp) {
    double q1, q3;
    if (ML_UNLIKELY (!d || !tmp || n <= 0)) return ml_make_nan ();
    q1 = ml_quantile (d, n, 0.25, tmp);
    q3 = ml_quantile (d, n, 0.75, tmp);
    return q3 - q1;
}
ML_API double ml_cv (const double *d, int n) {
    double m, s;
    if (ML_UNLIKELY (!d || n <= 0)) return ml_make_nan ();
    m = ml_mean (d, n);
    s = ml_stddev (d, n);
    if (m == 0.0) return ml_make_nan ();
    return s / m;
}
ML_API double ml_range (const double *d, int n) {
    double mn, mx;
    int i;
    if (ML_UNLIKELY (!d || n <= 0)) return ml_make_nan ();
    mn = mx = d [0];
    for (i = 1; i < n; i++) {
        if (d [i] < mn) mn = d [i];
        if (d [i] > mx) mx = d [i];
    }
    return mx - mn;
}
ML_API double ml_mad (const double *d, int n, double *tmp) {
    double med, s = 0.0;
    int i;
    if (ML_UNLIKELY (!d || !tmp || n <= 0)) return ml_make_nan ();
    med = ml_median (d, n, tmp);
    for (i = 0; i < n; i++) tmp [i] = ml_fabs (d [i] - med);
    s = ml_median (tmp, n, tmp);
    return s;
}
ML_API double ml_trimmed_mean (const double *d, int n, double alpha, double *tmp) {
    int cut, i, m;
    double s = 0.0;
    if (ML_UNLIKELY (!d || !tmp || n <= 0)) return ml_make_nan ();
    if (alpha < 0.0) alpha = 0.0;
    if (alpha >= 0.5) alpha = 0.5 - 1e-9;
    for (i = 0; i < n; i++) tmp [i] = d [i];
    se_sort (tmp, n);
    cut = (int) (alpha * (double) n);
    m = n - 2 * cut;
    if (m <= 0) return ml_make_nan ();
    for (i = 0; i < m; i++) s += tmp [cut + i];
    return s / (double) m;
}
ML_API ml_status_t ml_cov_matrix (const double *const *cols, int n, int p,
                                  double *cov) {
    int i, j, k;
    if (ML_UNLIKELY (!cols || !cov || n < 2 || p <= 0)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < p; i++) {
        for (j = i; j < p; j++) {
            double mi = ml_mean (cols [i], n), mj = ml_mean (cols [j], n), s = 0.0;
            for (k = 0; k < n; k++)
                s += (cols [i][k] - mi) * (cols [j][k] - mj);
            s /= (double) (n - 1);
            cov [(size_t) i * (size_t) p + (size_t) j] = s;
            cov [(size_t) j * (size_t) p + (size_t) i] = s;
        }
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_corr_matrix (const double *const *cols, int n, int p,
                                   double *corr) {
    int i, j, k;
    if (ML_UNLIKELY (!cols || !corr || n < 2 || p <= 0)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < p; i++) {
        for (j = i; j < p; j++) {
            double mi = ml_mean (cols [i], n), mj = ml_mean (cols [j], n);
            double si = 0.0, sj = 0.0, sij = 0.0, r;
            for (k = 0; k < n; k++) {
                double a = cols [i][k] - mi, b = cols [j][k] - mj;
                si += a * a;
                sj += b * b;
                sij += a * b;
            }
            r = (si > 0.0 && sj > 0.0) ? sij / ml_sqrt (si * sj) : 0.0;
            corr [(size_t) i * (size_t) p + (size_t) j] = r;
            corr [(size_t) j * (size_t) p + (size_t) i] = r;
        }
    }
    return ML_SUCCESS;
}
/* ------------------------------------------------------ distributions ---- */

ML_API double ml_lognormal_pdf (double x, double mu, double sigma) {
    double z;
    if (x <= 0.0 || sigma <= 0.0) return 0.0;
    z = (ml_log (x) - mu) / sigma;
    return ml_exp (-0.5 * z * z) / (x * sigma * ml_sqrt (2.0 * SE_PI));
}
ML_API double ml_lognormal_cdf (double x, double mu, double sigma) {
    if (x <= 0.0 || sigma <= 0.0) return (x <= 0.0) ? 0.0 : ml_make_nan ();
    return ml_normal_cdf ((ml_log (x) - mu) / sigma, 0.0, 1.0);
}
ML_API double ml_lognormal_inv (double p, double mu, double sigma) {
    if (p <= 0.0 || p >= 1.0 || sigma <= 0.0) return ml_make_nan ();
    return ml_exp (mu + sigma * ml_normal_inv (p, 0.0, 1.0));
}
ML_API double ml_weibull_pdf (double x, double shape, double scale) {
    double t;
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0) return 0.0;
    t = x / scale;
    return (shape / scale) * ml_pow (t, shape - 1.0) * ml_exp (-ml_pow (t, shape));
}
ML_API double ml_weibull_cdf (double x, double shape, double scale) {
    if (x <= 0.0) return 0.0;
    if (shape <= 0.0 || scale <= 0.0) return ml_make_nan ();
    return 1.0 - ml_exp (-ml_pow (x / scale, shape));
}
ML_API double ml_cauchy_pdf (double x, double loc, double scale) {
    double z;
    if (scale <= 0.0) return 0.0;
    z = (x - loc) / scale;
    return 1.0 / (SE_PI * scale * (1.0 + z * z));
}
ML_API double ml_cauchy_cdf (double x, double loc, double scale) {
    if (scale <= 0.0) return ml_make_nan ();
    return 0.5 + ml_atan2 (x - loc, scale) / SE_PI;
}
ML_API double ml_uniform_cdf (double x, double a, double b) {
    if (b <= a) return ml_make_nan ();
    if (x <= a) return 0.0;
    if (x >= b) return 1.0;
    return (x - a) / (b - a);
}
ML_API double ml_beta_pdf (double x, double a, double b) {
    double lb;
    if (x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0) return 0.0;
    lb = ml_lgamma (a) + ml_lgamma (b) - ml_lgamma (a + b);
    if (x == 0.0) return (a < 1.0) ? ml_make_inf (1) : ((a == 1.0) ? ml_exp (-lb) : 0.0);
    if (x == 1.0) return (b < 1.0) ? ml_make_inf (1) : ((b == 1.0) ? ml_exp (-lb) : 0.0);
    return ml_exp ((a - 1.0) * ml_log (x) + (b - 1.0) * ml_log (1.0 - x) - lb);
}
ML_API double ml_gamma_pdf (double x, double shape, double rate) {
    if (x < 0.0 || shape <= 0.0 || rate <= 0.0) return 0.0;
    if (x == 0.0) return (shape < 1.0) ? ml_make_inf (1) : ((shape == 1.0) ? rate : 0.0);
    return ml_exp (shape * ml_log (rate) + (shape - 1.0) * ml_log (x) - rate * x - ml_lgamma (shape));
}
ML_API double ml_chi2_pdf (double x, double k) {
    return ml_gamma_pdf (x, k / 2.0, 0.5);
}
ML_API double ml_student_t_pdf (double x, double nu) {
    double lb;
    if (nu <= 0.0) return ml_make_nan ();
    lb = ml_lgamma ((nu + 1.0) / 2.0) - ml_lgamma (nu / 2.0) - 0.5 * ml_log (nu * SE_PI);
    return ml_exp (lb - ((nu + 1.0) / 2.0) * ml_log (1.0 + x * x / nu));
}
ML_API double ml_f_pdf (double x, double d1, double d2) {
    double lb;
    if (x < 0.0 || d1 <= 0.0 || d2 <= 0.0) return 0.0;
    lb = ml_lgamma ((d1 + d2) / 2.0) - ml_lgamma (d1 / 2.0) - ml_lgamma (d2 / 2.0) + 0.5 * d1 * ml_log (d1) + 0.5 * d2 * ml_log (d2);
    return ml_exp (lb + (d1 / 2.0 - 1.0) * ml_log (x) - ((d1 + d2) / 2.0) * ml_log (d2 + d1 * x));
}
ML_API double ml_poisson_cdf (int k, double lambda) {
    double sum = 0.0, lz;
    int i;
    if (k < 0) return 0.0;
    if (lambda <= 0.0) return ml_make_nan ();
    lz = ml_log (lambda);
    for (i = 0; i <= k; i++)
        sum += ml_exp (-lambda + (double) i * lz - ml_lgamma ((double) i + 1.0));
    return sum > 1.0 ? 1.0 : sum;
}
ML_API double ml_binomial_cdf (int n, int k, double p) {
    double sum = 0.0;
    int i;
    if (k < 0) return 0.0;
    if (k >= n) return 1.0;
    if (n < 0 || p < 0.0 || p > 1.0) return ml_make_nan ();
    if (p == 0.0) return 1.0;
    if (p == 1.0) return 0.0;
    /* Exact recurrence for small n: start at P(0), iterate P(i+1) = P(i) * ... */
    if (n <= 1000) {
        double term = ml_exp ((double) n * ml_log1p (-p));
        double s = term;
        for (i = 0; i < k; i++) {
            term *= ((double) (n - i) / (double) (i + 1)) * (p / (1.0 - p));
            s += term;
        }
        return s > 1.0 ? 1.0 : s;
    }
    for (i = 0; i <= k; i++) {
        double lp = se_log_choose (n, i);
        if (p > 0.0) lp += (double) i * ml_log (p);
        else if (i > 0) continue;
        if (p < 1.0) lp += (double) (n - i) * ml_log (1.0 - p);
        else if (i < n) continue;
        sum += ml_exp (lp);
    }
    return sum > 1.0 ? 1.0 : sum;
}
ML_API double ml_bernoulli_pmf (int k, double p) {
    if (p < 0.0 || p > 1.0) return ml_make_nan ();
    if (k == 1) return p;
    if (k == 0) return 1.0 - p;
    return 0.0;
}
ML_API double ml_negbinom_pmf (int k, double r, double p) {
    if (k < 0 || r <= 0.0 || p <= 0.0 || p > 1.0) return 0.0;
    return ml_exp (se_log_choose ((int) ((double) k + r - 1.0 + 0.5), k) + r * ml_log (p) + (double) k * ml_log (1.0 - p));
}
ML_API double ml_hypergeometric_pmf (int k, int N, int K, int n) {
    double lp;
    if (k < 0 || k > K || k > n || n - k > N - K) return 0.0;
    /* Exact integer path for small counts (avoids log-gamma cancellation). */
    if (N <= 50) {
        double num = 1.0, den = 1.0;
        for (int i = 0; i < k; i++) {
            num *= (double) (K - i);
            den *= (double) (i + 1);
        }
        double cKk = num / den;
        num = 1.0; den = 1.0;
        for (int i = 0; i < n - k; i++) {
            num *= (double) (N - K - i);
            den *= (double) (i + 1);
        }
        double cNKnk = num / den;
        num = 1.0; den = 1.0;
        for (int i = 0; i < n; i++) {
            num *= (double) (N - i);
            den *= (double) (i + 1);
        }
        double cNn = num / den;
        return (cKk * cNKnk) / cNn;
    }
    lp = se_log_choose (K, k) + se_log_choose (N - K, n - k) - se_log_choose (N, n);
    return ml_exp (lp);
}
ML_API double ml_multivariate_normal_pdf (const double *x, const double *mu,
                                          const double *Sigma, int n,
                                          double *work) {
    double det = 1.0, quad = 0.0, norm;
    int i, j, k;
    if (ML_UNLIKELY (!x || !mu || !Sigma || !work || n <= 0 || n > 16))
        return ml_make_nan ();
    /* Copy Sigma to work as A, invert in place via Gauss-Jordan. */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            work [(size_t) i * (size_t) n + (size_t) j] = Sigma [(size_t) i * (size_t) n + (size_t) j];
    {
        double inv [16 * 16];
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                inv [(size_t) i * (size_t) n + (size_t) j] = (i == j) ? 1.0 : 0.0;
        for (k = 0; k < n; k++) {
            int piv = k;
            double mx = ml_fabs (work [(size_t) k * (size_t) n + (size_t) k]);
            for (i = k + 1; i < n; i++) {
                double v = ml_fabs (work [(size_t) i * (size_t) n + (size_t) k]);
                if (v > mx) {
                    mx = v;
                    piv = i;
                }
            }
            if (mx == 0.0) return 0.0;
            if (piv != k) {
                for (j = 0; j < n; j++) {
                    double t = work [(size_t) k * (size_t) n + (size_t) j];
                    work [(size_t) k * (size_t) n + (size_t) j] = work [(size_t) piv * (size_t) n + (size_t) j];
                    work [(size_t) piv * (size_t) n + (size_t) j] = t;
                    t = inv [(size_t) k * (size_t) n + (size_t) j];
                    inv [(size_t) k * (size_t) n + (size_t) j] = inv [(size_t) piv * (size_t) n + (size_t) j];
                    inv [(size_t) piv * (size_t) n + (size_t) j] = t;
                }
                det = -det;
            }
            {
                double d = work [(size_t) k * (size_t) n + (size_t) k];
                det *= d;
                for (j = 0; j < n; j++) {
                    work [(size_t) k * (size_t) n + (size_t) j] /= d;
                    inv [(size_t) k * (size_t) n + (size_t) j] /= d;
                }
            }
            for (i = 0; i < n; i++) {
                double f;
                if (i == k) continue;
                f = work [(size_t) i * (size_t) n + (size_t) k];
                for (j = 0; j < n; j++) {
                    work [(size_t) i * (size_t) n + (size_t) j] -=
                        f * work [(size_t) k * (size_t) n + (size_t) j];
                    inv [(size_t) i * (size_t) n + (size_t) j] -=
                        f * inv [(size_t) k * (size_t) n + (size_t) j];
                }
            }
        }
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                double dx = x [i] - mu [i], dy = x [j] - mu [j];
                quad += dx * inv [(size_t) i * (size_t) n + (size_t) j] * dy;
            }
    }
    norm = 1.0 / ml_sqrt (ml_pow (2.0 * SE_PI, (double) n) * ml_fabs (det));
    return norm * ml_exp (-0.5 * quad);
}
/* -------------------------------------------------- hypothesis tests ----- */

ML_API ml_status_t ml_ttest_1samp (const double *d, int n, double mu0,
                                   double *t, double *p) {
    double m, s;
    if (ML_UNLIKELY (!d || n < 2 || !t || !p)) return ML_ERR_INVALID_ARG;
    m = ml_mean (d, n);
    s = ml_stddev (d, n);
    if (s == 0.0) return ML_ERR_SINGULAR;
    *t = (m - mu0) / (s / ml_sqrt ((double) n));
    *p = se_twosided_t (*t, n - 1);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ttest_2samp (const double *x, int nx, const double *y,
                                   int ny, double *t, double *p) {
    double mx, my, s2x, s2y, sp2, den;
    if (ML_UNLIKELY (!x || !y || nx < 2 || ny < 2 || !t || !p))
        return ML_ERR_INVALID_ARG;
    mx = ml_mean (x, nx);
    my = ml_mean (y, ny);
    s2x = ml_variance_s (x, nx);
    s2y = ml_variance_s (y, ny);
    sp2 = ((double) (nx - 1) * s2x + (double) (ny - 1) * s2y) / (double) (nx + ny - 2);
    den = ml_sqrt (sp2 * (1.0 / (double) nx + 1.0 / (double) ny));
    if (den == 0.0) return ML_ERR_SINGULAR;
    *t = (mx - my) / den;
    *p = se_twosided_t (*t, nx + ny - 2);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ttest_welch (const double *x, int nx, const double *y,
                                   int ny, double *t, double *p) {
    double mx, my, vx, vy, den, df;
    if (ML_UNLIKELY (!x || !y || nx < 2 || ny < 2 || !t || !p))
        return ML_ERR_INVALID_ARG;
    mx = ml_mean (x, nx);
    my = ml_mean (y, ny);
    vx = ml_variance_s (x, nx) / (double) nx;
    vy = ml_variance_s (y, ny) / (double) ny;
    den = ml_sqrt (vx + vy);
    if (den == 0.0) return ML_ERR_SINGULAR;
    *t = (mx - my) / den;
    df = (vx + vy) * (vx + vy) / (vx * vx / (double) (nx - 1) + vy * vy / (double) (ny - 1));
    *p = se_twosided_t (*t, (int) (df + 0.5));
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ttest_paired (const double *x, const double *y, int n,
                                    double *t, double *p) {
    double m = 0.0, m2 = 0.0, md, sd;
    int i;
    if (ML_UNLIKELY (!x || !y || n < 2 || !t || !p)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        double d = x [i] - y [i];
        double delta = d - m;
        m += delta / (double) (i + 1);
        m2 += delta * (d - m);
    }
    md = m;
    sd = ml_sqrt (m2 / (double) (n - 1));
    if (sd == 0.0) return ML_ERR_SINGULAR;
    *t = md / (sd / ml_sqrt ((double) n));
    *p = se_twosided_t (*t, n - 1);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ztest_mean (const double *d, int n, double mu0,
                                  double sigma, double *z, double *p) {
    double m;
    if (ML_UNLIKELY (!d || n < 1 || sigma <= 0.0 || !z || !p))
        return ML_ERR_INVALID_ARG;
    m = ml_mean (d, n);
    *z = (m - mu0) / (sigma / ml_sqrt ((double) n));
    *p = se_twosided_normal (*z);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_f_test_var (const double *x, int nx, const double *y,
                                  int ny, double *F, double *p) {
    double s2x, s2y, f;
    if (ML_UNLIKELY (!x || !y || nx < 2 || ny < 2 || !F || !p))
        return ML_ERR_INVALID_ARG;
    s2x = ml_variance_s (x, nx);
    s2y = ml_variance_s (y, ny);
    if (s2y == 0.0) return ML_ERR_SINGULAR;
    f = s2x / s2y;
    *F = f;
    *p = 2.0 * (1.0 - ml_f_cdf (f, nx - 1, ny - 1));
    if (*p > 1.0) *p = 1.0;
    if (f < 1.0) *p = 2.0 * (1.0 - ml_f_cdf (1.0 / f, ny - 1, nx - 1));
    return ML_SUCCESS;
}
ML_API ml_status_t ml_chi2_gof (const int *obs, const double *exp, int k,
                                double *chi2, double *p) {
    double c = 0.0;
    int i;
    if (ML_UNLIKELY (!obs || !exp || k < 2 || !chi2 || !p)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < k; i++) {
        double e = exp [i], o = (double) obs [i];
        if (e <= 0.0) return ML_ERR_INVALID_ARG;
        c += (o - e) * (o - e) / e;
    }
    *chi2 = c;
    *p = 1.0 - ml_chi2_cdf (c, k - 1);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_anova_oneway (const double *const *groups,
                                    const int *sizes, int g, double *F,
                                    double *p) {
    double grand = 0.0, ssb = 0.0, ssw = 0.0;
    int i, j, N = 0;
    if (ML_UNLIKELY (!groups || !sizes || g < 2 || !F || !p))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < g; i++) N += sizes [i];
    for (i = 0; i < g; i++)
        for (j = 0; j < sizes [i]; j++) grand += groups [i][j];
    grand /= (double) N;
    for (i = 0; i < g; i++) {
        double mi = ml_mean (groups [i], sizes [i]);
        ssb += (double) sizes [i] * (mi - grand) * (mi - grand);
        for (j = 0; j < sizes [i]; j++) {
            double d = groups [i][j] - mi;
            ssw += d * d;
        }
    }
    if (ssw == 0.0) return ML_ERR_SINGULAR;
    {
        double msb = ssb / (double) (g - 1);
        double msw = ssw / (double) (N - g);
        *F = msb / msw;
        *p = 1.0 - ml_f_cdf (*F, g - 1, N - g);
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_mann_whitney_u (const double *x, int nx, const double *y,
                                      int ny, double *U, double *p) {
    double r1 = 0.0, mu, sigma, u1;
    int i, j, N;
    if (ML_UNLIKELY (!x || !y || nx < 1 || ny < 1 || !U || !p))
        return ML_ERR_INVALID_ARG;
    N = nx + ny;
    for (i = 0; i < nx; i++) {
        double rank = 1.0;
        for (j = 0; j < nx; j++)
            if (j != i && x [j] < x [i]) rank += 1.0;
        for (j = 0; j < ny; j++)
            if (y [j] < x [i]) rank += 1.0;
        for (j = 0; j < nx; j++)
            if (j != i && x [j] == x [i] && j < i) rank += 0.5;
        for (j = 0; j < ny; j++)
            if (y [j] == x [i]) rank += 0.5;
        r1 += rank;
    }
    u1 = r1 - (double) nx * (double) (nx + 1) / 2.0;
    *U = u1;
    mu = (double) nx * (double) ny / 2.0;
    sigma = ml_sqrt ((double) nx * (double) ny * (double) (N + 1) / 12.0);
    if (sigma == 0.0) return ML_ERR_SINGULAR;
    *p = se_twosided_normal ((u1 - mu) / sigma);
    return ML_SUCCESS;
}
/* ------------------------------------------------- confidence intervals - */

ML_API ml_status_t ml_ci_mean (const double *d, int n, double level, double *lo,
                               double *hi) {
    double m, se, tc;
    if (ML_UNLIKELY (!d || n < 2 || !lo || !hi)) return ML_ERR_INVALID_ARG;
    if (level <= 0.0 || level >= 1.0) return ML_ERR_INVALID_ARG;
    m = ml_mean (d, n);
    se = ml_stddev (d, n) / ml_sqrt ((double) n);
    tc = ml_student_t_inv ((1.0 + level) / 2.0, n - 1);
    *lo = m - tc * se;
    *hi = m + tc * se;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ci_proportion (int k, int n, double level, double *lo,
                                     double *hi) {
    double z, ph, den, center, half;
    if (ML_UNLIKELY (n <= 0 || !lo || !hi)) return ML_ERR_INVALID_ARG;
    if (level <= 0.0 || level >= 1.0) return ML_ERR_INVALID_ARG;
    z = ml_normal_inv ((1.0 + level) / 2.0, 0.0, 1.0);
    ph = (double) k / (double) n;
    den = 1.0 + z * z / (double) n;
    center = (ph + z * z / (2.0 * (double) n)) / den;
    half = z * ml_sqrt (ph * (1.0 - ph) / (double) n + z * z / (4.0 * (double) n * (double) n)) / den;
    *lo = center - half;
    *hi = center + half;
    if (*lo < 0.0) *lo = 0.0;
    if (*hi > 1.0) *hi = 1.0;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ci_variance (const double *d, int n, double level,
                                   double *lo, double *hi) {
    double s2, chi_lo, chi_hi;
    if (ML_UNLIKELY (!d || n < 2 || !lo || !hi)) return ML_ERR_INVALID_ARG;
    if (level <= 0.0 || level >= 1.0) return ML_ERR_INVALID_ARG;
    s2 = ml_variance_s (d, n);
    chi_lo = ml_chi2_inv ((1.0 - level) / 2.0, n - 1);
    chi_hi = ml_chi2_inv ((1.0 + level) / 2.0, n - 1);
    *lo = (double) (n - 1) * s2 / chi_hi;
    *hi = (double) (n - 1) * s2 / chi_lo;
    return ML_SUCCESS;
}
/* ------------------------------------------------------- regression ------ */

ML_API ml_status_t ml_multiple_regression (const double *X, const double *y,
                                           int n, int p, double *beta,
                                           double *sigma2) {
    double A [64 * 64], b [64], XtX [64 * 64], Xty [64];
    int i, j, k;
    if (ML_UNLIKELY (!X || !y || !beta || n <= p || p <= 0 || p > 64))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < p; i++) {
        Xty [i] = 0.0;
        for (k = 0; k < n; k++) Xty [i] += X [(size_t) k * (size_t) p + (size_t) i] * y [k];
        for (j = 0; j < p; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++)
                s += X [(size_t) k * (size_t) p + (size_t) i] *
                     X [(size_t) k * (size_t) p + (size_t) j];
            XtX [(size_t) i * (size_t) p + (size_t) j] = s;
        }
    }
    for (i = 0; i < p * p; i++) A [i] = XtX [i];
    for (i = 0; i < p; i++) b [i] = Xty [i];
    if (se_solve (A, b, beta, p) != ML_SUCCESS) return ML_ERR_SINGULAR;
    if (sigma2) {
        double sse = 0.0;
        for (k = 0; k < n; k++) {
            double fit = 0.0;
            for (j = 0; j < p; j++)
                fit += X [(size_t) k * (size_t) p + (size_t) j] * beta [j];
            sse += (y [k] - fit) * (y [k] - fit);
        }
        *sigma2 = sse / (double) (n - p);
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ridge_regression (const double *X, const double *y, int n,
                                        int p, double lambda, double *beta) {
    double A [64 * 64], b [64];
    int i, j, k;
    if (ML_UNLIKELY (!X || !y || !beta || n <= 0 || p <= 0 || p > 64))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < p; i++) {
        b [i] = 0.0;
        for (k = 0; k < n; k++) b [i] += X [(size_t) k * (size_t) p + (size_t) i] * y [k];
        for (j = 0; j < p; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++)
                s += X [(size_t) k * (size_t) p + (size_t) i] *
                     X [(size_t) k * (size_t) p + (size_t) j];
            A [(size_t) i * (size_t) p + (size_t) j] = s + ((i == j) ? lambda : 0.0);
        }
    }
    return se_solve (A, b, beta, p);
}
ML_API ml_status_t ml_wls (const double *X, const double *y, const double *w,
                           int n, int p, double *beta) {
    double A [64 * 64], b [64];
    int i, j, k;
    if (ML_UNLIKELY (!X || !y || !w || !beta || n <= 0 || p <= 0 || p > 64))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < p; i++) {
        b [i] = 0.0;
        for (k = 0; k < n; k++)
            b [i] += w [k] * X [(size_t) k * (size_t) p + (size_t) i] * y [k];
        for (j = 0; j < p; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++)
                s += w [k] * X [(size_t) k * (size_t) p + (size_t) i] *
                     X [(size_t) k * (size_t) p + (size_t) j];
            A [(size_t) i * (size_t) p + (size_t) j] = s;
        }
    }
    return se_solve (A, b, beta, p);
}
ML_API ml_status_t ml_logistic_regression (const double *X, const double *y,
                                           int n, int p, double *beta,
                                           int max_iter, double tol) {
    int i, j, k, it;
    if (ML_UNLIKELY (!X || !y || !beta || n <= 0 || p <= 0 || p > 64))
        return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 50;
    if (tol <= 0.0) tol = 1e-8;
    for (i = 0; i < p; i++) beta [i] = 0.0;
    for (it = 0; it < max_iter; it++) {
        double A [64 * 64], g [64], delta [64], mx = 0.0;
        for (i = 0; i < p * p; i++) A [i] = 0.0;
        for (i = 0; i < p; i++) g [i] = 0.0;
        for (k = 0; k < n; k++) {
            double eta = 0.0, pr, wgt;
            for (j = 0; j < p; j++)
                eta += X [(size_t) k * (size_t) p + (size_t) j] * beta [j];
            pr = 1.0 / (1.0 + ml_exp (-eta));
            wgt = pr * (1.0 - pr);
            if (wgt < 1e-10) wgt = 1e-10;
            for (i = 0; i < p; i++) {
                g [i] += X [(size_t) k * (size_t) p + (size_t) i] * (y [k] - pr);
                for (j = 0; j < p; j++)
                    A [(size_t) i * (size_t) p + (size_t) j] +=
                        wgt * X [(size_t) k * (size_t) p + (size_t) i] *
                        X [(size_t) k * (size_t) p + (size_t) j];
            }
        }
        if (se_solve (A, g, delta, p) != ML_SUCCESS) return ML_ERR_SINGULAR;
        for (i = 0; i < p; i++) {
            beta [i] += delta [i];
            if (ml_fabs (delta [i]) > mx) mx = ml_fabs (delta [i]);
        }
        if (mx < tol) break;
    }
    return ML_SUCCESS;
}
/* ------------------------------------------------------------- MLE ------- */

ML_API ml_status_t ml_mle_normal (const double *d, int n, double *mu,
                                  double *sigma2) {
    double m, s = 0.0;
    int i;
    if (ML_UNLIKELY (!d || n < 1 || !mu)) return ML_ERR_INVALID_ARG;
    m = ml_mean (d, n);
    for (i = 0; i < n; i++) {
        double dd = d [i] - m;
        s += dd * dd;
    }
    *mu = m;
    if (sigma2) *sigma2 = s / (double) n;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_mle_exponential (const double *d, int n, double *lambda) {
    double m;
    if (ML_UNLIKELY (!d || n < 1 || !lambda)) return ML_ERR_INVALID_ARG;
    m = ml_mean (d, n);
    if (m <= 0.0) return ML_ERR_INVALID_ARG;
    *lambda = 1.0 / m;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_mle_poisson (const int *d, int n, double *lambda) {
    long sum = 0;
    int i;
    if (ML_UNLIKELY (!d || n < 1 || !lambda)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) sum += d [i];
    *lambda = (double) sum / (double) n;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_mle_bernoulli (const int *d, int n, double *p) {
    long sum = 0;
    int i;
    if (ML_UNLIKELY (!d || n < 1 || !p)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) sum += d [i];
    *p = (double) sum / (double) n;
    return ML_SUCCESS;
}
/* --------------------------------------------------------- bootstrap ----- */

static uint64_t se_rng (uint64_t *s) {
    uint64_t x = *s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
ML_API ml_status_t ml_bootstrap_ci_mean (const double *d, int n, int B,
                                         double level, uint64_t seed,
                                         double *lo, double *hi) {
    double means [2048];
    int i, b;
    if (ML_UNLIKELY (!d || n < 1 || !lo || !hi)) return ML_ERR_INVALID_ARG;
    if (level <= 0.0 || level >= 1.0) return ML_ERR_INVALID_ARG;
    if (B < 1) B = 1;
    if (B > 2048) B = 2048;
    for (b = 0; b < B; b++) {
        double s = 0.0;
        for (i = 0; i < n; i++) {
            uint64_t r = se_rng (&seed);
            int idx = (int) (r % (uint64_t) n);
            s += d [idx];
        }
        means [b] = s / (double) n;
    }
    {
        double a = (1.0 - level) / 2.0;
        *lo = ml_quantile (means, B, a, means);
        *hi = ml_quantile (means, B, 1.0 - a, means);
    }
    return ML_SUCCESS;
}
