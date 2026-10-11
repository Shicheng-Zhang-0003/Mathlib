#include "ml_compiler.h"
#include "ml_modeling.h"
#include "ml_core.h"
#include "ml_exp_log.h"
#include "ml_stats_inv.h"

/* Applied modeling: Markov, SIR, Leslie, M/M/1, conjugate Bayes, chaos.
 * Closed forms where they exist; tiny verified solvers elsewhere. */

#define MLM_NMAX 256

ML_API int ml_markov_is_stochastic(const double *P, int n, double tol) {
    int i;
    int j;
    if (!P || n <= 0 || n > MLM_NMAX) return 0;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-12;
    for (i = 0; i < n; i++) {
        long double rs = 0.0L;
        for (j = 0; j < n; j++) {
            double v = P[(size_t)i * (size_t)n + (size_t)j];
            if (ml_isnan(v) || !ml_isfinite(v) || v < 0.0) return 0;
            rs += (long double)v;
        }
        if (ml_fabs((double)(rs - 1.0L)) > tol * (1.0 + ml_fabs((double)rs))) return 0;
    }
    return 1;
}

ML_API ml_status_t ml_markov_step(const double *P, const double *pi,
                                  double *out, int n) {
    int i;
    int j;
    if (!P || !pi || !out || n <= 0 || n > MLM_NMAX) return ML_ERR_INVALID_ARG;
    if (!ml_markov_is_stochastic(P, n, 1e-9)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        if (ml_isnan(pi[i]) || !ml_isfinite(pi[i]) || pi[i] < 0.0) return ML_ERR_NAN_INPUT;
    }
    for (j = 0; j < n; j++) {
        long double s = 0.0L;
        for (i = 0; i < n; i++) s += (long double)pi[i] * (long double)P[(size_t)i * (size_t)n + (size_t)j];
        out[j] = (double)s;
        if (!ml_isfinite(out[j])) return ML_ERR_OVERFLOW;
    }
    return ML_SUCCESS;
}

ML_API ml_status_t ml_markov_stationary(const double *P, double *pi, int n,
                                        int max_iter, double tol) {
    /* Power iteration from uniform; relative 1-norm stop. Aperiodic
     * irreducible chains converge; periodic/degenerate ones exhaust and
     * report CONVERGENCE (never silent SUCCESS). */
    double v[256];
    double w[256];
    int i;
    int it;
    if (!P || !pi || n <= 0 || n > MLM_NMAX) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-12;
    if (max_iter <= 0) max_iter = 10000;
    if (max_iter > 1000000) max_iter = 1000000;
    if (!ml_markov_is_stochastic(P, n, 1e-9)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) v[i] = 1.0 / (double)n;
    for (it = 0; it < max_iter; it++) {
        long double diff = 0.0L;
        int j;
        for (j = 0; j < n; j++) {
            long double s = 0.0L;
            for (i = 0; i < n; i++) s += (long double)v[i] * (long double)P[(size_t)i * (size_t)n + (size_t)j];
            w[j] = (double)s;
            if (!ml_isfinite(w[j])) return ML_ERR_OVERFLOW;
        }
        for (j = 0; j < n; j++) diff += (w[j] > v[j]) ? (long double)(w[j] - v[j]) : (long double)(v[j] - w[j]);
        for (j = 0; j < n; j++) v[j] = w[j];
        if ((double)diff <= tol * (1.0 + 0.0)) {
            long double sum = 0.0L;
            for (j = 0; j < n; j++) sum += (long double)v[j];
            if (!(sum > 0.0L) || !ml_isfinite((double)sum)) return ML_ERR_CONVERGENCE;
            for (j = 0; j < n; j++) pi[j] = (double)((long double)v[j] / sum);
            return ML_SUCCESS;
        }
    }
    return ML_ERR_CONVERGENCE;
}

ML_API double ml_sir_r0(double beta, double gamma) {
    if (ml_isnan(beta) || ml_isnan(gamma)) return ml_make_nan();
    if (!ml_isfinite(beta) || !ml_isfinite(gamma)) return ml_make_nan();
    if (!(beta >= 0.0) || !(gamma > 0.0)) return ml_make_nan();
    return beta / gamma;
}

ML_API double ml_sir_final_size(double r0) {
    /* Attack ratio Z solves Z = 1 - exp(-R0 Z); Newton from Z = 0.5.
     * R0 <= 1: no epidemic, Z = 0. */
    double z;
    int i;
    if (ml_isnan(r0) || !ml_isfinite(r0)) return ml_make_nan();
    if (r0 <= 1.0) {
        if (r0 < 0.0) return ml_make_nan();
        return 0.0;
    }
    if (r0 > 1e12) return 1.0;
    z = 0.5;
    for (i = 0; i < 100; i++) {
        double e = ml_exp(-r0 * z);
        double f = 1.0 - e - z;
        double df = r0 * e - 1.0;
        double dz;
        if (!ml_isfinite(f)) return ml_make_nan();
        if (ml_fabs(f) < 1e-15) break;
        if (df == 0.0) return ml_make_nan();
        dz = f / df;
        z -= dz;
        if (z <= 0.0) z = 0.5 * (1.0 - e);
        if (z >= 1.0) z = 0.999;
        if (ml_fabs(dz) < 1e-15 * (1.0 + ml_fabs(z))) break;
    }
    if (!(z >= 0.0) || !(z <= 1.0) || !ml_isfinite(z)) return ml_make_nan();
    return z;
}

ML_API double ml_sir_peak(double r0) {
    /* I_max = 1 - (1 + ln R0)/R0 for R0 > 1, else 0. */
    double o;
    if (ml_isnan(r0) || !ml_isfinite(r0)) return ml_make_nan();
    if (r0 < 0.0) return ml_make_nan();
    if (r0 <= 1.0) return 0.0;
    o = 1.0 - (1.0 + ml_log(r0)) / r0;
    if (!ml_isfinite(o) || o < 0.0) return 0.0;
    return o;
}

ML_API ml_status_t ml_leslie_step(const double *L, const double *pop,
                                  double *out, int n) {
    int i;
    int j;
    if (!L || !pop || !out || n <= 0 || n > MLM_NMAX) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        if (ml_isnan(pop[i]) || !ml_isfinite(pop[i]) || pop[i] < 0.0) return ML_ERR_NAN_INPUT;
        for (j = 0; j < n; j++) {
            double v = L[(size_t)i * (size_t)n + (size_t)j];
            if (ml_isnan(v) || !ml_isfinite(v) || v < 0.0) return ML_ERR_NAN_INPUT;
            if (i > 0 && j != i - 1 && v != 0.0) return ML_ERR_INVALID_ARG;
            if (i > 0 && j == i - 1 && v > 1.0) return ML_ERR_INVALID_ARG;
        }
    }
    for (i = 0; i < n; i++) {
        long double s = 0.0L;
        for (j = 0; j < n; j++) s += (long double)L[(size_t)i * (size_t)n + (size_t)j] * (long double)pop[j];
        out[i] = (double)s;
        if (!ml_isfinite(out[i])) return ML_ERR_OVERFLOW;
    }
    return ML_SUCCESS;
}

ML_API ml_status_t ml_leslie_growth(const double *L, int n, double *lambda) {
    /* Dominant eigenvalue by power iteration (Leslie matrices are
     * primitive under mild conditions; exhaustion -> CONVERGENCE). */
    double v[256];
    double w[256];
    double lam = 0.0;
    int i;
    int it;
    if (!L || !lambda || n <= 0 || n > MLM_NMAX) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n * n; i++) {
        if (ml_isnan(L[i]) || !ml_isfinite(L[i]) || L[i] < 0.0) return ML_ERR_NAN_INPUT;
    }
    for (i = 0; i < n; i++) v[i] = 1.0;
    for (it = 0; it < 10000; it++) {
        long double nw = 0.0L;
        long double rq = 0.0L;
        int j;
        for (i = 0; i < n; i++) {
            long double s = 0.0L;
            for (j = 0; j < n; j++) s += (long double)L[(size_t)i * (size_t)n + (size_t)j] * (long double)v[j];
            w[i] = (double)s;
            if (!ml_isfinite(w[i])) return ML_ERR_OVERFLOW;
        }
        for (i = 0; i < n; i++) nw += (long double)w[i] * (long double)w[i];
        nw = __builtin_sqrtl(nw);
        if (!(nw > 0.0L)) return ML_ERR_SINGULAR;
        for (i = 0; i < n; i++) w[i] = (double)((long double)w[i] / nw);
        for (i = 0; i < n; i++) {
            long double s = 0.0L;
            for (j = 0; j < n; j++) s += (long double)L[(size_t)i * (size_t)n + (size_t)j] * (long double)w[j];
            rq += (long double)w[i] * s;
        }
        if (ml_fabs((double)rq - lam) <= 1e-10 * (1.0 + ml_fabs((double)rq))) {
            *lambda = (double)rq;
            return ML_SUCCESS;
        }
        lam = (double)rq;
        for (i = 0; i < n; i++) v[i] = w[i];
    }
    return ML_ERR_CONVERGENCE;
}

static int mlm_queue_ok(double lambda, double mu) {
    if (ml_isnan(lambda) || ml_isnan(mu)) return 0;
    if (!ml_isfinite(lambda) || !ml_isfinite(mu)) return 0;
    if (!(lambda >= 0.0) || !(mu > 0.0) || !(lambda < mu)) return 0;
    return 1;
}

ML_API double ml_mm1_rho(double lambda, double mu) {
    if (!mlm_queue_ok(lambda, mu)) return ml_make_nan();
    return lambda / mu;
}

ML_API double ml_mm1_L(double lambda, double mu) {
    /* L = rho/(1-rho). */
    double rho;
    if (!mlm_queue_ok(lambda, mu)) return ml_make_nan();
    rho = lambda / mu;
    return rho / (1.0 - rho);
}

ML_API double ml_mm1_W(double lambda, double mu) {
    /* W = 1/(mu-lambda). */
    if (!mlm_queue_ok(lambda, mu)) return ml_make_nan();
    return 1.0 / (mu - lambda);
}

ML_API double ml_mm1_wait(double lambda, double mu, double t) {
    /* P(wait > t) = exp(-(mu-lambda) t) for t >= 0. */
    double o;
    if (!mlm_queue_ok(lambda, mu)) return ml_make_nan();
    if (ml_isnan(t) || !ml_isfinite(t) || t < 0.0) return ml_make_nan();
    o = ml_exp(-(mu - lambda) * t);
    return ml_isfinite(o) ? o : ml_make_nan();
}

ML_API ml_status_t ml_beta_binomial_post(uint64_t a0, uint64_t b0,
                                         uint64_t succ, uint64_t trials,
                                         double *a1, double *b1) {
    long double na;
    long double nb;
    if (ML_UNLIKELY(!a1 || !b1)) return ML_ERR_INVALID_ARG;
    if (a0 < 1 || b0 < 1 || a0 > 1000000000000ULL || b0 > 1000000000000ULL)
        return ML_ERR_INVALID_ARG;
    if (succ > trials || trials > 1000000000000ULL) return ML_ERR_INVALID_ARG;
    na = (long double)a0 + (long double)succ;
    nb = (long double)b0 + (long double)(trials - succ);
    if (na > 1.7976931348623157e308L || nb > 1.7976931348623157e308L) return ML_ERR_OVERFLOW;
    *a1 = (double)na;
    *b1 = (double)nb;
    return ML_SUCCESS;
}

ML_API double ml_beta_post_mean(double a1, double b1) {
    double o;
    if (ml_isnan(a1) || ml_isnan(b1)) return ml_make_nan();
    if (!ml_isfinite(a1) || !ml_isfinite(b1) || !(a1 > 0.0) || !(b1 > 0.0))
        return ml_make_nan();
    o = a1 / (a1 + b1);
    return ml_isfinite(o) ? o : ml_make_nan();
}

ML_API ml_status_t ml_gamma_poisson_post(double a0, double b0, uint64_t k,
                                         double texp, double *a1, double *b1) {
    /* Gamma(a0,b0) prior + Poisson(k in exposure texp):
     * a1 = a0+k, b1 = b0+texp (rate form). */
    long double na;
    long double nb;
    if (ML_UNLIKELY(!a1 || !b1)) return ML_ERR_INVALID_ARG;
    if (ml_isnan(a0) || ml_isnan(b0) || ml_isnan(texp)) return ML_ERR_NAN_INPUT;
    if (!ml_isfinite(a0) || !ml_isfinite(b0) || !ml_isfinite(texp)) return ML_ERR_NAN_INPUT;
    if (!(a0 > 0.0) || !(b0 > 0.0) || !(texp >= 0.0)) return ML_ERR_INVALID_ARG;
    if (k > 1000000000000ULL) return ML_ERR_INVALID_ARG;
    na = (long double)a0 + (long double)k;
    nb = (long double)b0 + (long double)texp;
    if (!ml_isfinite((double)na) || !ml_isfinite((double)nb)) return ML_ERR_OVERFLOW;
    *a1 = (double)na;
    *b1 = (double)nb;
    return ML_SUCCESS;
}

ML_API ml_status_t ml_beta_credible(double a, double b, double alpha,
                                    double *lo, double *hi) {
    /* Equal-tailed interval via the staged beta_inv (reuses stat core). */
    double tail;
    if (ML_UNLIKELY(!lo || !hi)) return ML_ERR_INVALID_ARG;
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(alpha)) return ML_ERR_NAN_INPUT;
    if (!ml_isfinite(a) || !ml_isfinite(b)) return ML_ERR_NAN_INPUT;
    if (!(a > 0.0) || !(b > 0.0)) return ML_ERR_INVALID_ARG;
    if (!(alpha > 0.0) || !(alpha < 1.0) || !ml_isfinite(alpha)) return ML_ERR_INVALID_ARG;
    tail = alpha * 0.5;
    *lo = ml_beta_inv(tail, a, b);
    *hi = ml_beta_inv(1.0 - tail, a, b);
    if (!ml_isfinite(*lo) || !ml_isfinite(*hi)) return ML_ERR_CONVERGENCE;
    return ML_SUCCESS;
}

ML_API double ml_logistic_lyapunov(double r, double x0, int transient, int steps) {
    /* lambda = mean ln|r(1-2x)| over steps after transient burn-in. */
    long double acc = 0.0L;
    long double x;
    int i;
    if (ml_isnan(r) || ml_isnan(x0)) return ml_make_nan();
    if (!ml_isfinite(r) || !ml_isfinite(x0)) return ml_make_nan();
    if (r < 0.0 || r > 4.0 || x0 <= 0.0 || x0 >= 1.0) return ml_make_nan();
    if (transient < 0 || transient > 1000000 || steps < 1 || steps > 1000000)
        return ml_make_nan();
    x = (long double)x0;
    for (i = 0; i < transient; i++) {
        x = (long double)r * x * (1.0L - x);
        if (x <= 0.0L || x >= 1.0L) return ml_make_nan();
    }
    for (i = 0; i < steps; i++) {
        long double d = (long double)r * (1.0L - 2.0L * x);
        long double ad = (d < 0.0L) ? -d : d;
        if (!(ad > 0.0L) || !ml_isfinite((double)ad)) return ml_make_nan();
        acc += __builtin_logl(ad);
        x = (long double)r * x * (1.0L - x);
        if (x <= 0.0L || x >= 1.0L) {
            if (i < steps - 1) {
                /* Absorbed at a fixed point: exponent is -inf in theory;
                 * report a large negative finite value fail-loud. */
                return -745.0;
            }
            break;
        }
    }
    {
        double o = (double)(acc / (long double)steps);
        return ml_isfinite(o) ? o : ml_make_nan();
    }
}
