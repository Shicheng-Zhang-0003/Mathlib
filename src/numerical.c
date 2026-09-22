#include "ml_compiler.h"
#include "ml_numerical.h"
#include "ml_trig.h"
#include "ml_exp_log.h"

/* v11S CLOSURE IP-15: numerical methods robustness */

ML_API double ml_newton_raphson(ml_func_t f, ml_func_t df, double x0, double epsilon, int max_iter) {
    if (ML_UNLIKELY(f == NULL || df == NULL)) {
        return ml_make_nan();
    }

    /* Reject NaN/Inf/non-positive tolerance: NaN<=0 is false, so an
     * explicit finite-positive test is required. */
    if (ML_UNLIKELY(!(epsilon > 0.0) || !ml_isfinite(epsilon) || max_iter <= 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!ml_isfinite(x0))) {
        return ml_make_nan();
    }

    double x = x0;

    for (int i = 0; i < max_iter; i++) {
        double fx = f(x);

        if (ML_UNLIKELY(!ml_isfinite(fx))) {
            return ml_make_nan();
        }

        if (fx == 0.0) {
            return x;
        }

        double dfx = df(x);

        if (ML_UNLIKELY(!ml_isfinite(dfx))) {
            return ml_make_nan();
        }

        if (ML_UNLIKELY(dfx == 0.0)) {
            return ml_make_nan();
        }

        double x_next = x - fx / dfx;

        if (ML_UNLIKELY(!ml_isfinite(x_next))) {
            return ml_make_nan();
        }

        /* Scale-aware stop: absolute eps never triggers for |x|~1e12.
         * tol*(1+|x|) behaves absolute near 0, relative far away. */
        if (ml_fabs(x_next - x) <= epsilon * (1.0 + ml_fabs(x_next))) {
            return x_next;
        }

        x = x_next;
    }

    /* max_iter exhausted without convergence: NaN, not last iterate. */
    return ml_make_nan();
}

ML_API double ml_bisection(ml_func_t f, double a, double b, double epsilon, int max_iter) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!(epsilon > 0.0) || !ml_isfinite(epsilon) || max_iter <= 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isinf(a) || ml_isinf(b))) {
        return ml_make_nan();
    }

    if (a > b) {
        double t = a;
        a = b;
        b = t;
    }

    double fa = f(a);
    double fb = f(b);

    if (ML_UNLIKELY(ml_isnan(fa) || ml_isnan(fb))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!ml_isfinite(fa) || !ml_isfinite(fb))) {
        return ml_make_nan();
    }

    if (fa == 0.0) {
        return a;
    }

    if (fb == 0.0) {
        return b;
    }

    if (ml_signbit(fa) == ml_signbit(fb)) {
        return ml_make_nan();
    }

    double c = a;

    for (int i = 0; i < max_iter; i++) {
        /* Overflow-safe midpoint: a*0.5+b*0.5 never overflows, unlike
         * (a+b)*0.5 (same-sign overflow) or a+(b-a)/2 (wide-span overflow). */
        c = a * 0.5 + b * 0.5;

        /*
         * If the midpoint stops changing, further bisection is impossible
         * in double precision. Return the best current estimate.
         */
        if (c == a || c == b) {
            return c;
        }

        double fc = f(c);

        if (ML_UNLIKELY(ml_isnan(fc) || !ml_isfinite(fc))) {
            return ml_make_nan();
        }

        if (fc == 0.0 || ml_fabs(b - a) < epsilon) {
            return c;
        }

        if (ml_signbit(fa) != ml_signbit(fc)) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    /* Exhausted max_iter with bracket still wide: not converged. */
    return ml_make_nan();
}

ML_API double ml_brent(ml_func_t f, double a, double b, double tol, int max_iter) {
    /* Brent-Dekker: bisection safety + secant/IQI speed. No derivatives.
     * Requires f(a),f(b) opposite signs (or zero endpoint). */
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(!(tol > 0.0) || !ml_isfinite(tol) || max_iter <= 0)) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    {
        double fa = f(a), fb = f(b);
        if (ML_UNLIKELY(!ml_isfinite(fa) || !ml_isfinite(fb))) return ml_make_nan();
        if (fa == 0.0) return a;
        if (fb == 0.0) return b;
        if (ml_signbit(fa) == ml_signbit(fb)) return ml_make_nan();
        {
            double c = a, fc = fa, d = b - a, e = d;
            for (int i = 0; i < max_iter; i++) {
                if (ml_signbit(fb) == ml_signbit(fc)) {
                    c = a; fc = fa; d = e = b - a;
                }
                if (ml_fabs(fc) < ml_fabs(fb)) {
                    a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;
                }
                {
                    double m = 0.5 * (c - b);
                    double tol_act = 2.0 * 2.220446049250313e-16 * ml_fabs(b) + 0.5 * tol;
                    if (ml_fabs(m) <= tol_act || fb == 0.0) return b;
                    double p, q;
                    if (ml_fabs(e) >= tol_act && ml_fabs(fa) > ml_fabs(fb)) {
                        double s = fb / fa;
                        double r, t;
                        if (a == c) {
                            /* Secant. */
                            p = 2.0 * m * s; q = 1.0 - s;
                        } else {
                            /* Inverse quadratic interpolation. */
                            r = fa / fc; t = fb / fc;
                            p = s * (2.0 * m * r * (r - t) - (b - a) * (t - 1.0));
                            q = (r - 1.0) * (t - 1.0) * (s - 1.0);
                        }
                        if (p > 0.0) q = -q; else p = -p;
                        if (2.0 * p < 3.0 * m * q - ml_fabs(tol_act * q) &&
                            p < ml_fabs(0.5 * e * q)) {
                            e = d; d = p / q;
                        } else {
                            d = e = m;
                        }
                    } else {
                        d = e = m;
                    }
                    a = b; fa = fb;
                    if (ml_fabs(d) > tol_act) b += d;
                    else b += (m > 0.0) ? tol_act : -tol_act;
                    fb = f(b);
                    if (ML_UNLIKELY(!ml_isfinite(fb))) return ml_make_nan();
                }
            }
            /* Converged paths return b inside the loop; reaching here
             * means max_iter exhausted without tolerance. */
            return ml_make_nan();
        }
    }
}

ML_API double ml_kepler(double M, double e, double tol) {
    if (ML_UNLIKELY(ml_isnan(M) || ml_isnan(e) || ml_isnan(tol))) return ml_make_nan();
    if (ML_UNLIKELY(!(e >= 0.0) || !(e < 1.0))) return ml_make_nan();
    if (ML_UNLIKELY(!(tol > 0.0) || !ml_isfinite(tol))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(M))) return ml_make_nan();
    {
        double E = (e < 0.8) ? M : 3.14159265358979323846;
        for (int i = 0; i < 100; i++) {
            double f = E - e * ml_sin(E) - M;
            double fp = 1.0 - e * ml_cos(E);
            if (!ml_isfinite(f) || !ml_isfinite(fp) || fp == 0.0) return ml_make_nan();
            {
                double d = f / fp;
                E -= d;
                if (ml_fabs(d) < tol) return E;
            }
        }
        return ml_make_nan();
    }
}

ML_API double ml_derivative(ml_func_t f, double x, double h) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(h == 0.0 || ml_isnan(h) || ml_isinf(h))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(x) || ml_isinf(x))) {
        return ml_make_nan();
    }

    double xp = x + h;
    double xm = x - h;

    if (ML_UNLIKELY(!ml_isfinite(xp) || !ml_isfinite(xm))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(xp == xm)) {
        return ml_make_nan();
    }

    double fp = f(xp);
    double fm = f(xm);

    if (ML_UNLIKELY(!ml_isfinite(fp) || !ml_isfinite(fm))) {
        return ml_make_nan();
    }

    return (fp - fm) / (2.0 * h);
}

ML_API double ml_second_derivative(ml_func_t f, double x, double h) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(h == 0.0 || ml_isnan(h) || ml_isinf(h))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(x) || ml_isinf(x))) {
        return ml_make_nan();
    }

    double xp = x + h;
    double xm = x - h;

    if (ML_UNLIKELY(!ml_isfinite(xp) || !ml_isfinite(xm))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(xp == xm)) {
        return ml_make_nan();
    }

    double fp = f(xp);
    double f0 = f(x);
    double fm = f(xm);

    if (ML_UNLIKELY(!ml_isfinite(fp) || !ml_isfinite(f0) || !ml_isfinite(fm))) {
        return ml_make_nan();
    }

    double denom = h * h;

    if (ML_UNLIKELY(denom == 0.0 || !ml_isfinite(denom))) {
        return ml_make_nan();
    }

    return (fp - 2.0 * f0 + fm) / denom;
}

ML_API double ml_integral_simpson(ml_func_t f, double a, double b, int n) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    /*
     * Simpson's rule strictly requires an even number of subintervals.
     * Reject invalid grids explicitly instead of silently changing them.
     * Cap n to prevent CPU DoS (10M+ integrand evaluations).
     */
    if (ML_UNLIKELY(n < 2 || (n % 2) == 1 || n > 10000000)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isinf(a) || ml_isinf(b))) {
        return ml_make_nan();
    }

    if (a == b) {
        return 0.0;
    }

    double h = (b - a) / (double)n;

    if (ML_UNLIKELY(h == 0.0 || ml_isnan(h) || ml_isinf(h))) {
        return ml_make_nan();
    }

    double fa = f(a);
    double fb = f(b);

    if (ML_UNLIKELY(!ml_isfinite(fa) || !ml_isfinite(fb))) {
        return ml_make_nan();
    }

    /* Kahan-compensated Simpson accumulation: naive summation loses
     * O(n*eps) accuracy for large n; compensation keeps O(eps). */
    double result = fa + fb;
    double comp = 0.0;

    for (int i = 1; i < n; i++) {
        double x = a + (double)i * h;
        double fx = f(x);

        if (ML_UNLIKELY(!ml_isfinite(fx))) {
            return ml_make_nan();
        }

        {
            double w = ((i % 2) == 0) ? 2.0 * fx : 4.0 * fx;
            double y = w - comp;
            double t = result + y;
            comp = (t - result) - y;
            result = t;
        }
    }

    {
        double prod = result * h;
        if (ML_UNLIKELY(!ml_isfinite(prod))) {
            return ml_make_nan();
        }
        return prod / 3.0;
    }
}

static double ml_adaptive_simpson_rec(ml_func_t f, double a, double fa,
                                      double m, double fm, double b, double fb,
                                      double whole, double tol, int depth) {
    double h = (b - a) * 0.5;
    double lm = a * 0.5 + m * 0.5;
    double rm = m * 0.5 + b * 0.5;
    double flm = f(lm);
    double frm = f(rm);
    if (!ml_isfinite(flm) || !ml_isfinite(frm)) return ml_make_nan();
    {
        double left = (fa + 4.0 * flm + fm) * h / 6.0;
        double right = (fm + 4.0 * frm + fb) * h / 6.0;
        double delta = left + right - whole;
        if (depth <= 0 || ml_fabs(delta) <= 15.0 * tol) {
            return left + right + delta / 15.0;
        }
        {
            double l = ml_adaptive_simpson_rec(f, a, fa, lm, flm, m, fm,
                                               left, tol * 0.5, depth - 1);
            double r = ml_adaptive_simpson_rec(f, m, fm, rm, frm, b, fb,
                                               right, tol * 0.5, depth - 1);
            if (!ml_isfinite(l) || !ml_isfinite(r)) return ml_make_nan();
            return l + r;
        }
    }
}

ML_API double ml_integral_adaptive(ml_func_t f, double a, double b, double tol, int max_depth) {
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(!(tol > 0.0) || !ml_isfinite(tol))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    if (a == b) return 0.0;
    if (max_depth <= 0) max_depth = 12;
    if (max_depth > 20) max_depth = 20;
    {
        double m = a * 0.5 + b * 0.5;
        double fa = f(a), fm = f(m), fb = f(b);
        if (ML_UNLIKELY(!ml_isfinite(fa) || !ml_isfinite(fm) || !ml_isfinite(fb))) {
            return ml_make_nan();
        }
        {
            double whole = (fa + 4.0 * fm + fb) * (b - a) / 6.0;
            if (!ml_isfinite(whole)) return ml_make_nan();
            return ml_adaptive_simpson_rec(f, a, fa, m, fm, b, fb, whole, tol, max_depth);
        }
    }
}

ML_API double ml_integral_tanhsinh(ml_func_t f, double a, double b, double tol) {
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(!(tol > 0.0) || !ml_isfinite(tol))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    if (a == b) return 0.0;
    if (ML_UNLIKELY(!(a < b))) return ml_make_nan();
    {
        double half = (b - a) * 0.5, mid = a * 0.5 + b * 0.5;
        double prev = 0.0;
        double h = 1.0;
        for (int lev = 0; lev < 10; lev++) {
            double s = 0.0, comp = 0.0;
            int K = (int)(4.0 / h) + 20;
            if (K > 2000) K = 2000;
            for (int k = -K; k <= K; k++) {
                double t = (double)k * h;
                double sh = (ml_exp(t) - ml_exp(-t)) * 0.5;
                double ch = (ml_exp(t) + ml_exp(-t)) * 0.5;
                double u = sh * 1.57079632679489661923;
                double eu = ml_exp(u), enu = ml_exp(-u);
                double tanh_u, w;
                if (!ml_isfinite(eu)) continue;
                tanh_u = (eu - enu) / (eu + enu);
                {
                    double ch_u = (eu + enu) * 0.5;
                    w = half * 1.57079632679489661923 * ch / (ch_u * ch_u);
                }
                {
                    double x = mid + half * tanh_u;
                    if (!(x > a && x < b)) continue;
                    {
                        double fx = f(x);
                        if (!ml_isfinite(fx)) {
                            if (w * ml_fabs(fx) < tol * 1e-3) continue;
                            return ml_make_nan();
                        }
                        {
                            double term = fx * w * h;
                            double yy = term - comp;
                            double tt = s + yy;
                            comp = (tt - s) - yy;
                            s = tt;
                        }
                    }
                }
            }
            if (lev > 2 && ml_fabs(s - prev) <= tol) return s;
            prev = s;
            h *= 0.5;
        }
        return prev;
    }
}

ML_API double ml_integral_trapezoid(ml_func_t f, double a, double b, int n) {
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(n < 1 || n > 10000000)) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b))) return ml_make_nan();
    if (a == b) return 0.0;
    {
        double h = (b - a) / (double)n;
        if (ML_UNLIKELY(!ml_isfinite(h) || h == 0.0)) return ml_make_nan();
        double fa = f(a), fb = f(b);
        if (ML_UNLIKELY(!ml_isfinite(fa) || !ml_isfinite(fb))) return ml_make_nan();
        {
            double s = (fa + fb) * 0.5, comp = 0.0;
            for (int i = 1; i < n; i++) {
                double fx = f(a + (double)i * h);
                if (ML_UNLIKELY(!ml_isfinite(fx))) return ml_make_nan();
                {
                    double yy = fx - comp;
                    double tt = s + yy;
                    comp = (tt - s) - yy;
                    s = tt;
                }
            }
            {
                double r = s * h;
                return ml_isfinite(r) ? r : ml_make_nan();
            }
        }
    }
}

ML_API double ml_arith_nth(double a1, double d, long long n) {
    if (ml_isnan(a1) || ml_isnan(d) || n <= 0) return ml_make_nan();
    if (!ml_isfinite(a1) || !ml_isfinite(d)) return ml_make_nan();
    return a1 + (double)(n - 1) * d;
}

ML_API double ml_arith_sum(double a1, double d, long long n) {
    if (ml_isnan(a1) || ml_isnan(d) || n <= 0) return ml_make_nan();
    if (!ml_isfinite(a1) || !ml_isfinite(d)) return ml_make_nan();
    {
        double nn = (double)n;
        double an = a1 + (double)(n - 1) * d;
        return nn * (a1 + an) * 0.5;
    }
}

ML_API double ml_geom_nth(double a1, double r, long long n) {
    if (ml_isnan(a1) || ml_isnan(r) || n <= 0) return ml_make_nan();
    if (!ml_isfinite(a1) || !ml_isfinite(r)) return ml_make_nan();
    return a1 * ml_pow(r, (double)(n - 1));
}

ML_API double ml_geom_sum(double a1, double r, long long n) {
    if (ml_isnan(a1) || ml_isnan(r) || n <= 0) return ml_make_nan();
    if (!ml_isfinite(a1) || !ml_isfinite(r)) return ml_make_nan();
    if (r == 1.0) return a1 * (double)n;
    if (ml_fabs(r - 1.0) < 1e-8) {
        /* (1-r^n)/(1-r) cancels for r~1: use expm1/log1p path.
         * r^n = exp(n*log(r)); 1-r^n = -expm1(n*log1p(r-1)). */
        double lr = ml_log1p(r - 1.0);
        if (!ml_isfinite(lr)) return a1 * (double)n;
        {
            double e = ml_expm1((double)n * lr);
            if (!ml_isfinite(e)) return a1 * (1.0 - ml_pow(r, (double)n)) / (1.0 - r);
            return a1 * (-e) / (1.0 - r);
        }
    }
    return a1 * (1.0 - ml_pow(r, (double)n)) / (1.0 - r);
}

ML_API double ml_sum_k(long long n) {
    if (n <= 0) return ml_make_nan();
    {
        double nn = (double)n;
        if (nn > 4294967295.0) return ml_make_nan();
        return nn * (nn + 1.0) * 0.5;
    }
}

ML_API double ml_sum_k2(long long n) {
    if (n <= 0) return ml_make_nan();
    {
        double nn = (double)n;
        if (nn > 2097151.0) return ml_make_nan();
        return nn * (nn + 1.0) * (2.0 * nn + 1.0) / 6.0;
    }
}

ML_API double ml_sum_k3(long long n) {
    if (n <= 0) return ml_make_nan();
    {
        double s = ml_sum_k(n);
        if (!ml_isfinite(s)) return ml_make_nan();
        return s * s;
    }
}

ML_API int ml_nim_win(const uint64_t *piles, int n) {
    if (!piles || n <= 0) return 0;
    {
        uint64_t x = 0;
        for (int i = 0; i < n; i++) x ^= piles[i];
        return (x != 0) ? 1 : 0;
    }
}

ML_API int ml_majorizes(const double *a, const double *b, int n, double tol) {
    /* Karamata: sorted-desc prefix sums A>=B, totals equal. O(n^2) select. */
    if (!a || !b || n <= 0) return 0;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-12;
    {
        double sa = 0.0, sb = 0.0;
        for (int i = 0; i < n; i++) {
            if (!ml_isfinite(a[i]) || !ml_isfinite(b[i])) return 0;
            sa += a[i]; sb += b[i];
        }
        if (ml_fabs(sa - sb) > tol * (1.0 + ml_fabs(sa))) return 0;
        for (int k = 1; k < n; k++) {
            /* k-th largest via partial select on temp copies. */
            double ta[256], tb[256];
            int nn = (n > 256) ? 256 : n;
            if (n > 256) return 0;
            for (int i = 0; i < n; i++) { ta[i] = a[i]; tb[i] = b[i]; }
            for (int i = 0; i < k; i++) {
                int ma = i, mb = i;
                for (int j = i + 1; j < n; j++) {
                    if (ta[j] > ta[ma]) ma = j;
                    if (tb[j] > tb[mb]) mb = j;
                }
                {
                    double t = ta[i]; ta[i] = ta[ma]; ta[ma] = t;
                    t = tb[i]; tb[i] = tb[mb]; tb[mb] = t;
                }
            }
            {
                double pa = 0.0, pb = 0.0;
                for (int i = 0; i < k; i++) { pa += ta[i]; pb += tb[i]; }
                if (pa + tol < pb) return 0;
            }
            (void)nn;
        }
        return 1;
    }
}
