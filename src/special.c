#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_integral.h"
#include "ml_math_extra.h"
#include "ml_special.h"
#include <limits.h>

#define SP_PI ML_PI
#define SP_EULER 0.577215664901532860606512090082402431
#define SP_MAX_ORDER 4096
static double sp_floor (double x) { return ml_floor (x); }
static int sp_is_int (double x) { return ml_isfinite (x) && x == sp_floor (x); }
static int sp_integer_parity (double x, int *parity) {
    if (!sp_is_int (x) || parity == NULL) return 0;
    /* Every representable binary64 integer at this magnitude is even. */
    if (ml_fabs (x) >= 0x1p53) *parity = 0;
    else *parity = ((long long) x % 2LL) != 0;
    return 1;
}
static int sp_near_integer (double x) {
    return !sp_is_int (x) && ml_fabs (x - ml_round (x)) < 1e-7;
}
static int sp_asymptotic_ok (double nu, double x) {
    double an = ml_fabs (nu);
    /* The fixed-order Hankel expansions require x to dominate nu^2. */
    return x >= 18.0 && x >= 4.0 * an * an;
}
/* ------------------------------------------------------------------ Bessel */

static double j_series (double nu, double x) {
    int sign = 1, k, converged = 0;
    double log_gamma = ml_lgamma_r (nu + 1.0, &sign);
    long double half = (long double) x * 0.5L;
    long double half2 = half * half;
    long double term = (long double) sign * __builtin_expl ((long double) nu * __builtin_logl (half) - (long double) log_gamma);
    long double sum = term, comp = 0.0L;
    for (k = 1; k < 400; k++) {
        long double y, next;
        term *= -half2 / ((long double) k * ((long double) k + (long double) nu));
        y = term - comp;
        next = sum + y;
        comp = (next - sum) - y;
        sum = next;
        if (__builtin_fabsl (term) < 1e-20L * __builtin_fabsl (sum)) {
            converged = 1;
            break;
        }
    }
    return converged ? (double) sum : ml_make_nan ();
}
static double i_series (double nu, double x) {
    int sign = 1, k, converged = 0;
    double log_gamma = ml_lgamma_r (nu + 1.0, &sign);
    long double half = (long double) x * 0.5L;
    long double half2 = half * half;
    long double term = (long double) sign * __builtin_expl ((long double) nu * __builtin_logl (half) - (long double) log_gamma);
    long double sum = term, comp = 0.0L;
    for (k = 1; k < 2000; k++) {
        long double y, next;
        term *= half2 / ((long double) k * ((long double) k + (long double) nu));
        y = term - comp;
        next = sum + y;
        comp = (next - sum) - y;
        sum = next;
        if (__builtin_fabsl (term) < 1e-20L * (__builtin_fabsl (sum) + 1e-300L)) {
            converged = 1;
            break;
        }
    }
    return converged ? (double) sum : ml_make_nan ();
}
static double j_asym (double nu, double x) {
    double mu = 4.0 * nu * nu;
    double a [11], px = 1.0 / x, xp = 1.0, P = 0.0, Q = 0.0, w, amp;
    int k;
    a [0] = 1.0;
    for (k = 1; k < 11; k++) {
        double t = (double) (2 * k - 1);
        a [k] = a [k - 1] * (mu - t * t) / ((double) k * 8.0);
    }
    for (k = 0; k < 11; k++) {
        double val = a [k] * xp;
        double sgn = (((k / 2) % 2) == 0) ? 1.0 : -1.0;
        if ((k % 2) == 0) P += sgn * val;
        else Q += sgn * val;
        xp *= px;
    }
    w = x - nu * SP_PI * 0.5 - SP_PI * 0.25;
    amp = ml_sqrt (2.0 / (SP_PI * x));
    return amp * (ml_cos (w) * P - ml_sin (w) * Q);
}
static double y_asym (double nu, double x) {
    double mu = 4.0 * nu * nu;
    double a [11], px = 1.0 / x, xp = 1.0, P = 0.0, Q = 0.0, w, amp;
    int k;
    a [0] = 1.0;
    for (k = 1; k < 11; k++) {
        double t = (double) (2 * k - 1);
        a [k] = a [k - 1] * (mu - t * t) / ((double) k * 8.0);
    }
    for (k = 0; k < 11; k++) {
        double val = a [k] * xp;
        double sgn = (((k / 2) % 2) == 0) ? 1.0 : -1.0;
        if ((k % 2) == 0) P += sgn * val;
        else Q += sgn * val;
        xp *= px;
    }
    w = x - nu * SP_PI * 0.5 - SP_PI * 0.25;
    amp = ml_sqrt (2.0 / (SP_PI * x));
    return amp * (ml_sin (w) * P + ml_cos (w) * Q);
}
static double i_asym (double nu, double x) {
    double mu = 4.0 * nu * nu, sum = 1.0, term = 1.0;
    int k;
    for (k = 1; k < 20; k++) {
        double t = (double) (2 * k - 1);
        term *= -(mu - t * t) / ((double) k * 8.0 * x);
        sum += term;
        if (ml_fabs (term) < 1e-18 * ml_fabs (sum)) break;
    }
    return ml_exp (x) / ml_sqrt (2.0 * SP_PI * x) * sum;
}
static double k_asym (double nu, double x) {
    double mu = 4.0 * nu * nu, sum = 1.0, term = 1.0;
    int k;
    for (k = 1; k < 20; k++) {
        double t = (double) (2 * k - 1);
        term *= (mu - t * t) / ((double) k * 8.0 * x);
        sum += term;
        if (ml_fabs (term) < 1e-18 * ml_fabs (sum)) break;
    }
    return ml_sqrt (SP_PI / (2.0 * x)) * ml_exp (-x) * sum;
}
static double j_integer (int n, double x) {
    int k;
    if (n < 0) return (n % 2) ? -j_integer (-n, x) : j_integer (-n, x);
    if (n == 0) return ml_bessel_j0 (x);
    if (n == 1) return ml_bessel_j1 (x);
    /* Direct series for small orders and small x: converges rapidly,
     * avoids Miller recurrence accumulation error (~3 ULP at J2(1)). */
    if (n <= 5 && x <= 10.0) {
        long double half = (long double) x * 0.5L;
        long double half2 = half * half;
        long double term = __builtin_powl (half, (long double) n);
        /* term = (x/2)^n / n! */
        long double fact = 1.0L;
        for (int f = 2; f <= n; f++) fact *= (long double) f;
        term /= fact;
        long double sum = term;
        int m;
        for (m = 1; m < 200; m++) {
            term *= -half2 / ((long double) m * (long double) (m + n));
            sum += term;
            if (__builtin_fabsl (term) < 1e-20L * __builtin_fabsl (sum)) break;
        }
        double r = (double) sum;
        return ml_isfinite (r) ? r : ml_make_nan ();
    }
    if (x > (double) n + 10.0) {
        double jm = ml_bessel_j0 (x), j = ml_bessel_j1 (x);
        for (k = 1; k < n; k++) {
            double jp = (2.0 * (double) k / x) * j - jm;
            jm = j;
            j = jp;
        }
        return j;
    }
    /* Miller downward recurrence, normalized by J0 + 2*sum(J_2k) = 1. */
    {
        int top = n + 32 + (int) ml_sqrt (40.0 * x);
        long double next = 0.0L, curr = 1.0L, sum_even = 0.0L, target = 0.0L;
        for (k = top; k >= 1; k--) {
            long double prev = (2.0L * (long double) k / (long double) x) * curr - next;
            int index = k - 1;
            long double mag = __builtin_fabsl (prev);
            if (__builtin_fabsl (curr) > mag) mag = __builtin_fabsl (curr);
            if (__builtin_fabsl (next) > mag) mag = __builtin_fabsl (next);
            if (index == n) target = prev;
            if (index >= 2 && (index % 2) == 0) sum_even += prev;
            if (mag > 1e100L) {
                prev *= 1e-100L;
                curr *= 1e-100L;
                next *= 1e-100L;
                sum_even *= 1e-100L;
                target *= 1e-100L;
            }
            next = curr;
            curr = prev;
        }
        {
            long double norm = curr + 2.0L * sum_even;
            if (norm == 0.0L || target != target) return ml_make_nan ();
            return (double) (target / norm);
        }
    }
}
static double y_integer (int n, double x) {
    int sgn = 1, k;
    double ym, y0;
    if (n < 0) {
        if (n % 2) sgn = -1;
        n = -n;
    }
    if (n == 0) return sgn * ml_bessel_y0 (x);
    if (n == 1) return sgn * ml_bessel_y1 (x);
    ym = ml_bessel_y0 (x);
    y0 = ml_bessel_y1 (x);
    for (k = 1; k < n; k++) {
        double y1 = (2.0 * (double) k / x) * y0 - ym;
        ym = y0;
        y0 = y1;
    }
    return sgn * y0;
}
static double k_integer (int n, double x) {
    int k;
    double km, k0;
    if (n < 0) n = -n;
    if (n == 0) return ml_bessel_k0 (x);
    if (n == 1) return ml_bessel_k1 (x);
    km = ml_bessel_k0 (x);
    k0 = ml_bessel_k1 (x);
    for (k = 1; k < n; k++) {
        double k1 = km + (2.0 * (double) k / x) * k0;
        km = k0;
        k0 = k1;
    }
    return k0;
}
static long double kv_integrand (long double t, long double nu, long double x) {
    long double c = __builtin_coshl (t);
    long double v = -x * c + nu * t - 0.693147180559945309417232121458176568L;
    v += __builtin_log1pl (__builtin_expl (-2.0L * nu * t));
    if (v < -11300.0L) return 0.0L;
    return __builtin_expl (v);
}
static long double kv_integral_rec (long double a, long double b,
                                    long double fa, long double fm, long double fb,
                                    long double whole, long double nu, long double x,
                                    long double tol, int depth) {
    long double m = (a + b) * 0.5L;
    long double lm = (a + m) * 0.5L, rm = (m + b) * 0.5L;
    long double flm = kv_integrand (lm, nu, x), frm = kv_integrand (rm, nu, x);
    long double left = (m - a) / 6.0L * (fa + 4.0L * flm + fm);
    long double right = (b - m) / 6.0L * (fm + 4.0L * frm + fb);
    long double both = left + right;
    if (depth <= 0 || __builtin_fabsl (both - whole) <= 15.0L * tol)
        return both + (both - whole) / 15.0L;
    return kv_integral_rec (a, m, fa, flm, fm, left, nu, x, tol * 0.5L, depth - 1) + kv_integral_rec (m, b, fm, frm, fb, right, nu, x, tol * 0.5L, depth - 1);
}
static double k_noninteger_quad (double nu, double x) {
    long double an = (long double) ml_fabs (nu), xx = (long double) x;
    long double log_est = an > 1.0L
                              ? an * (__builtin_logl (2.0L * an / xx) - 1.0L)
                              : -xx;
    long double top, a, fa, fm, fb, whole;
    int i;
    if (log_est > 710.0L) return ml_make_inf (0);
    if (xx > 11000.0L) return 0.0;
    top = __builtin_logl (2.0L * (an + 80.0L) / xx);
    for (i = 0; i < 12; i++) {
        long double next = __builtin_logl (2.0L * (an * top + 80.0L) / xx);
        if (__builtin_fabsl (next - top) < 1e-18L) break;
        top = next;
    }
    if (!(top > 0.0L) || top > 12000.0L) return ml_make_nan ();
    a = 0.0L;
    fa = kv_integrand (a, an, xx);
    fm = kv_integrand (top * 0.5L, an, xx);
    fb = kv_integrand (top, an, xx);
    whole = top / 6.0L * (fa + 4.0L * fm + fb);
    return (double) kv_integral_rec (a, top, fa, fm, fb, whole, an, xx,
                                     1e-18L * __builtin_fabsl (whole) + 1e-30L,
                                     26);
}
ML_API double ml_bessel_jv (double nu, double x) {
    int parity;
    if (ml_isnan (nu) || ml_isnan (x)) return ml_make_nan ();
    if (!ml_isfinite (nu) || ml_fabs (nu) > (double) SP_MAX_ORDER) return ml_make_nan ();
    if (nu < 0.0 && sp_integer_parity (nu, &parity)) {
        double v = ml_bessel_jv (-nu, x);
        return parity ? -v : v;
    }
    if (x < 0.0) {
        if (sp_integer_parity (nu, &parity)) {
            double v = ml_bessel_jv (nu, -x);
            return parity ? -v : v;
        }
        return ml_make_nan ();
    }
    if (ml_isinf (x)) return 0.0;
    if (x == 0.0) return (nu == 0.0) ? 1.0 : ((nu > 0.0 || sp_is_int (nu)) ? 0.0 : ml_make_nan ());
    if (sp_is_int (nu)) return j_integer ((int) nu, x);
    if (x < 18.0) return j_series (nu, x);
    if (!sp_asymptotic_ok (nu, x)) return ml_make_nan ();
    return j_asym (nu, x);
}
ML_API double ml_bessel_iv (double nu, double x) {
    int parity;
    if (ml_isnan (nu) || ml_isnan (x)) return ml_make_nan ();
    if (!ml_isfinite (nu) || ml_fabs (nu) > (double) SP_MAX_ORDER) return ml_make_nan ();
    if (nu < 0.0 && sp_is_int (nu)) return ml_bessel_iv (-nu, x);
    if (x < 0.0) {
        if (sp_integer_parity (nu, &parity)) {
            double v = ml_bessel_iv (nu, -x);
            return parity ? -v : v;
        }
        return ml_make_nan ();
    }
    if (ml_isinf (x)) return ml_make_inf (0);
    if (x == 0.0) return (nu == 0.0) ? 1.0 : ((nu > 0.0) ? 0.0 : ml_make_nan ());
    if (x < 710.0) return i_series (nu, x);
    if (!sp_asymptotic_ok (nu, x)) return ml_make_nan ();
    return i_asym (nu, x);
}
ML_API double ml_bessel_yv (double nu, double x) {
    if (ml_isnan (nu) || ml_isnan (x)) return ml_make_nan ();
    if (!ml_isfinite (nu) || ml_fabs (nu) > (double) SP_MAX_ORDER) return ml_make_nan ();
    if (ml_isinf (x) && x > 0.0) return 0.0;
    if (x <= 0.0) return ml_make_nan ();
    if (sp_is_int (nu)) return y_integer ((int) nu, x);
    if (sp_near_integer (nu)) return ml_make_nan ();
    if (sp_asymptotic_ok (nu, x)) return y_asym (nu, x);
    if (x >= 18.0) return ml_make_nan ();
    return (ml_cos (nu * SP_PI) * ml_bessel_jv (nu, x) - ml_bessel_jv (-nu, x)) / ml_sin (nu * SP_PI);
}
ML_API double ml_bessel_kv (double nu, double x) {
    if (ml_isnan (nu) || ml_isnan (x)) return ml_make_nan ();
    if (!ml_isfinite (nu) || ml_fabs (nu) > (double) SP_MAX_ORDER) return ml_make_nan ();
    if (ml_isinf (x) && x > 0.0) return 0.0;
    if (x <= 0.0) return ml_make_nan ();
    if (sp_is_int (nu)) return k_integer ((int) nu, x);
    if (sp_near_integer (nu)) return ml_make_nan ();
    if (sp_asymptotic_ok (nu, x)) return k_asym (nu, x);
    return k_noninteger_quad (nu, x);
}
ML_API double ml_bessel_jn (int n, double x) { return ml_bessel_jv ((double) n, x); }
ML_API double ml_bessel_yn (int n, double x) { return ml_bessel_yv ((double) n, x); }
ML_API double ml_sph_bessel_j (int l, double x) {
    if (l < 0 || l > SP_MAX_ORDER) return ml_make_nan ();
    if (x < 0.0) {
        double v = ml_sph_bessel_j (l, -x);
        return (l % 2) ? -v : v;
    }
    double a = (double) l + 0.5;
    if (x <= 0.0) return (x == 0.0 && l == 0) ? 1.0 : ((x == 0.0 && l > 0) ? 0.0 : ml_make_nan ());
    return ml_sqrt (SP_PI / (2.0 * x)) * ml_bessel_jv (a, x);
}
ML_API double ml_sph_bessel_y (int l, double x) {
    if (l < 0 || l > SP_MAX_ORDER) return ml_make_nan ();
    if (x < 0.0) {
        double v = ml_sph_bessel_y (l, -x);
        return (l % 2) ? v : -v;
    }
    double a = (double) l + 0.5;
    if (x == 0.0) return ml_make_inf (1);
    return ml_sqrt (SP_PI / (2.0 * x)) * ml_bessel_yv (a, x);
}
ML_API double ml_sph_bessel_i (int l, double x) {
    if (l < 0 || l > SP_MAX_ORDER) return ml_make_nan ();
    if (x < 0.0) {
        double v = ml_sph_bessel_i (l, -x);
        return (l % 2) ? -v : v;
    }
    if (ml_isinf (x)) return ml_make_inf (0);
    double a = (double) l + 0.5;
    if (x <= 0.0) return (x == 0.0) ? (l == 0 ? 1.0 : 0.0) : ml_make_nan ();
    return ml_sqrt (SP_PI / (2.0 * x)) * ml_bessel_iv (a, x);
}
ML_API double ml_sph_bessel_k (int l, double x) {
    if (l < 0 || l > SP_MAX_ORDER) return ml_make_nan ();
    double a = (double) l + 0.5;
    if (x == 0.0) return ml_make_inf (0);
    if (x < 0.0) return ml_make_nan ();
    return ml_sqrt (SP_PI / (2.0 * x)) * ml_bessel_kv (a, x);
}
/* -------------------------------------------------------------------- Airy */

#define SP_AIRY_BI0 0.6149266274460007
#define SP_AIRY_BIP0 0.4482883573538264
#define SP_AIRY_AIP0 -0.2588194037928068
ML_API double ml_airy_bi (double x) {
    double ax, zeta;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return x > 0.0 ? ml_make_inf (0) : ml_make_nan ();
    if (ml_fabs (x) < 1e-100) return SP_AIRY_BI0 + SP_AIRY_BIP0 * x;
    ax = ml_fabs (x);
    zeta = (2.0 / 3.0) * ax * ml_sqrt (ax);
    if (!ml_isfinite (zeta)) return x > 0.0 ? ml_make_inf (0) : ml_make_nan ();
    if (x > 0.0)
        return ml_sqrt (x / 3.0) * (ml_bessel_iv (-1.0 / 3.0, zeta) + ml_bessel_iv (1.0 / 3.0, zeta));
    return ml_sqrt (ax / 3.0) * (ml_bessel_jv (-1.0 / 3.0, zeta) - ml_bessel_jv (1.0 / 3.0, zeta));
}
ML_API double ml_airy_ai_prime (double x) {
    double ax, zeta;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return x > 0.0 ? 0.0 : ml_make_nan ();
    if (ml_fabs (x) < 1e-100) return SP_AIRY_AIP0;
    if (x > 0.0) {
        zeta = (2.0 / 3.0) * x * ml_sqrt (x);
        return -x / (SP_PI * ml_sqrt (3.0)) * ml_bessel_kv (2.0 / 3.0, zeta);
    }
    ax = -x;
    zeta = (2.0 / 3.0) * ax * ml_sqrt (ax);
    if (!ml_isfinite (zeta)) return ml_make_nan ();
    return -ax / 3.0 * (ml_bessel_jv (-2.0 / 3.0, zeta) - ml_bessel_jv (2.0 / 3.0, zeta));
}
ML_API double ml_airy_bi_prime (double x) {
    double ax, zeta;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return x > 0.0 ? ml_make_inf (0) : ml_make_nan ();
    if (ml_fabs (x) < 1e-100) return SP_AIRY_BIP0;
    ax = ml_fabs (x);
    zeta = (2.0 / 3.0) * ax * ml_sqrt (ax);
    if (!ml_isfinite (zeta)) return x > 0.0 ? ml_make_inf (0) : ml_make_nan ();
    if (x > 0.0)
        return x / ml_sqrt (3.0) * (ml_bessel_iv (-2.0 / 3.0, zeta) + ml_bessel_iv (2.0 / 3.0, zeta));
    return ax / ml_sqrt (3.0) * (ml_bessel_jv (-2.0 / 3.0, zeta) + ml_bessel_jv (2.0 / 3.0, zeta));
}
/* ----------------------------------------------- Exponential / trig integrals */

static double e1_series (double x) {
    long double sum = 0.0L, tk = 1.0L;
    long double xl = (long double) x;
    int k;
    for (k = 1; k < 400; k++) {
        tk *= xl / (long double) k;
        long double term = ((k % 2) ? 1.0L : -1.0L) * tk / (long double) k;
        sum += term;
        if (__builtin_fabsl (tk / (long double) k) < 1e-18L * __builtin_fabsl (sum)) break;
    }
    long double r = -((long double) SP_EULER) - __builtin_logl (xl) + sum;
    double o = (double) r;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
static double e1_cf (double x) {
    double b = x + 1.0;
    double c = 1e300, d = 1.0 / b, h = d;
    int n;
    for (n = 1; n < 300; n++) {
        double a = -(double) n * (double) n;
        double del;
        b += 2.0;
        d = 1.0 / (a * d + b);
        c = b + a / c;
        del = c * d;
        h *= del;
        if (ml_fabs (del - 1.0) < 1e-16) break;
    }
    return ml_exp (-x) * h;
}
ML_API double ml_expint_e1 (double x) {
    if (ml_isnan (x)) return ml_make_nan ();
    if (x < 0.0) return ml_make_nan ();
    if (x == 0.0) return ml_make_inf (0);
    if (ml_isinf (x)) return 0.0;
    if (x <= 1.0) return e1_series (x);
    return e1_cf (x);
}
ML_API double ml_si (double x) {
    double ax = ml_fabs (x), s;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return ml_copysign (SP_PI * 0.5, x);
    if (x == 0.0) return 0.0;
    if (ax < 21.0) {
        double term = x, sum = x;
        int k;
        for (k = 1; k < 300; k++) {
            term *= -x * x / (((double) (2 * k)) * ((double) (2 * k + 1)));
            sum += term / (double) (2 * k + 1);
            if (ml_fabs (term / (double) (2 * k + 1)) < 1e-18 * ml_fabs (sum)) break;
        }
        return sum;
    }
    {
        double f = 1.0 / ax, g = 1.0 / (ax * ax), fs = f, gs = g;
        int k;
        for (k = 1; k < 50; k++) {
            double nf = f * (-(double) (2 * k) * (double) (2 * k - 1) / (ax * ax));
            double ng = g * (-(double) (2 * k + 1) * (double) (2 * k) / (ax * ax));
            if (ml_fabs (nf) >= ml_fabs (f) || ml_fabs (ng) >= ml_fabs (g)) break;
            f = nf;
            g = ng;
            fs += f;
            gs += g;
        }
        s = SP_PI * 0.5 - fs * ml_cos (ax) - gs * ml_sin (ax);
    }
    return (x < 0.0) ? -s : s;
}
ML_API double ml_ci (double x) {
    double ax = ml_fabs (x);
    if (ml_isnan (x)) return ml_make_nan ();
    if (x < 0.0) return ml_make_nan ();
    if (x == 0.0) return ml_make_inf (1);
    if (ml_isinf (x)) return 0.0;
    if (ax < 2.0) {
        double sum = 0.0, term = 1.0;
        int k;
        for (k = 1; k < 300; k++) {
            term *= -x * x / (((double) (2 * k - 1)) * ((double) (2 * k)));
            sum += term / (double) (2 * k);
            if (ml_fabs (term / (double) (2 * k)) < 1e-18 * ml_fabs (sum + 1.0)) break;
        }
        return SP_EULER + ml_log (ax) + sum;
    }
    if (ax < 21.0) {
        double sum = 0.0, term = 1.0;
        int k;
        for (k = 1; k < 500; k++) {
            term *= -ax * ax / (((double) (2 * k - 1)) * ((double) (2 * k)));
            sum += term / (double) (2 * k);
            if (ml_fabs (term / (double) (2 * k)) < 1e-18 * ml_fabs (sum + 1.0)) break;
        }
        return SP_EULER + ml_log (ax) + sum;
    }
    {
        double f = 1.0 / ax, g = 1.0 / (ax * ax), fs = f, gs = g;
        int k;
        for (k = 1; k < 50; k++) {
            double nf = f * (-(double) (2 * k) * (double) (2 * k - 1) / (ax * ax));
            double ng = g * (-(double) (2 * k + 1) * (double) (2 * k) / (ax * ax));
            if (ml_fabs (nf) >= ml_fabs (f) || ml_fabs (ng) >= ml_fabs (g)) break;
            f = nf;
            g = ng;
            fs += f;
            gs += g;
        }
        return fs * ml_sin (ax) - gs * ml_cos (ax);
    }
}
ML_API double ml_shi (double x) {
    double term = x, sum = x;
    int k;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return x;
    if (x == 0.0) return 0.0;
    for (k = 1; k < 400; k++) {
        term *= x * x / (((double) (2 * k)) * ((double) (2 * k + 1)));
        sum += term / (double) (2 * k + 1);
        if (ml_fabs (term / (double) (2 * k + 1)) < 1e-18 * ml_fabs (sum)) break;
    }
    return sum;
}
ML_API double ml_chi (double x) {
    double ax = ml_fabs (x), sum = 0.0, term = 1.0;
    int k;
    if (ml_isnan (x)) return ml_make_nan ();
    if (x < 0.0) return ml_make_nan ();
    if (x == 0.0) return ml_make_inf (1);
    if (ml_isinf (x)) return ml_make_inf (0);
    for (k = 1; k < 500; k++) {
        term *= ax * ax / (((double) (2 * k - 1)) * ((double) (2 * k)));
        sum += term / (double) (2 * k);
        if (ml_fabs (term / (double) (2 * k)) < 1e-18 * ml_fabs (sum + 1.0)) break;
    }
    return SP_EULER + ml_log (ax) + sum;
}
/* --------------------------------------------------- Fresnel / Dawson / erfcx */

ML_API double ml_fresnel_s (double x) {
    double ax, sum, term;
    int k;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x) || ml_fabs (x) >= 0x1p55) return ml_copysign (0.5, x);
    ax = ml_fabs (x);
    if (ax < 4.0) {
        double a2 = (SP_PI * 0.5) * (SP_PI * 0.5);
        sum = term = (SP_PI * 0.5) * ax * ax * ax / 3.0;
        for (k = 0; k < 120; k++) {
            double dk = (double) k;
            term *= -a2 * ax * ax * ax * ax * (4.0 * dk + 3.0) /
                    ((2.0 * dk + 3.0) * (2.0 * dk + 2.0) * (4.0 * dk + 7.0));
            sum += term;
            if (ml_fabs (term) < 1e-18 * ml_fabs (sum)) break;
        }
        return (x < 0.0) ? -sum : sum;
    }
    {
        double f = 1.0 / (SP_PI * ax), g = 1.0 / (SP_PI * SP_PI * ax * ax * ax);
        double tf = f, tg = g, xx = ax * ax, phase;
        for (k = 1; k < 20; k++) {
            double nf = tf * (-(double) (4 * k - 1) * (double) (4 * k - 3)) /
                        (SP_PI * SP_PI * xx * xx);
            double ng = tg * (-(double) (4 * k + 1) * (double) (4 * k - 1)) /
                        (SP_PI * SP_PI * xx * xx);
            if (ml_fabs (nf) >= ml_fabs (tf) || ml_fabs (ng) >= ml_fabs (tg)) break;
            tf = nf;
            tg = ng;
            f += tf;
            g += tg;
        }
        phase = (SP_PI * 0.5) * xx;
        sum = 0.5 - f * ml_cos (phase) - g * ml_sin (phase);
    }
    return (x < 0.0) ? -sum : sum;
}
ML_API double ml_fresnel_c (double x) {
    double ax, sum, term;
    int k;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x) || ml_fabs (x) >= 0x1p55) return ml_copysign (0.5, x);
    ax = ml_fabs (x);
    if (ax < 4.0) {
        double a2 = (SP_PI * 0.5) * (SP_PI * 0.5);
        sum = term = ax;
        for (k = 0; k < 120; k++) {
            double dk = (double) k;
            term *= -a2 * ax * ax * ax * ax * (4.0 * dk + 1.0) /
                    ((2.0 * dk + 2.0) * (2.0 * dk + 1.0) * (4.0 * dk + 5.0));
            sum += term;
            if (ml_fabs (term) < 1e-18 * ml_fabs (sum)) break;
        }
        return (x < 0.0) ? -sum : sum;
    }
    {
        double f = 1.0 / (SP_PI * ax), g = 1.0 / (SP_PI * SP_PI * ax * ax * ax);
        double tf = f, tg = g, xx = ax * ax, phase;
        for (k = 1; k < 20; k++) {
            double nf = tf * (-(double) (4 * k - 1) * (double) (4 * k - 3)) /
                        (SP_PI * SP_PI * xx * xx);
            double ng = tg * (-(double) (4 * k + 1) * (double) (4 * k - 1)) /
                        (SP_PI * SP_PI * xx * xx);
            if (ml_fabs (nf) >= ml_fabs (tf) || ml_fabs (ng) >= ml_fabs (tg)) break;
            tf = nf;
            tg = ng;
            f += tf;
            g += tg;
        }
        phase = (SP_PI * 0.5) * xx;
        sum = 0.5 + f * ml_sin (phase) - g * ml_cos (phase);
    }
    return (x < 0.0) ? -sum : sum;
}
ML_API double ml_dawson (double x) {
    double ax = ml_fabs (x), sum;
    int k;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ax < 4.0) {
        long double xx = (long double) x;
        long double lt = xx, ls = xx, comp = 0.0L;
        for (k = 1; k < 300; k++) {
            long double y, next;
            lt *= (-2.0L * xx * xx) / (long double) (2 * k + 1);
            y = lt - comp;
            next = ls + y;
            comp = (next - ls) - y;
            ls = next;
            if (__builtin_fabsl (lt) < 1e-20L * __builtin_fabsl (ls)) break;
        }
        return (double) ls;
    }
    {
        double y = 1.0 / (2.0 * ax), yp = y;
        sum = y;
        for (k = 1; k < 50; k++) {
            double next = yp * (double) (2 * k - 1) / (2.0 * ax * ax);
            if (ml_fabs (next) >= ml_fabs (yp)) break;
            yp = next;
            sum += yp;
            if (ml_fabs (yp) < 1e-18 * ml_fabs (sum)) break;
        }
        return (x < 0.0) ? -sum : sum;
    }
}
ML_API double ml_erfcx (double x) {
    if (ml_isnan (x)) return ml_make_nan ();
    if (x < 0.0) {
        double p = ml_erfcx (-x);
        return 2.0 * ml_exp (x * x) - p;
    }
    if (x < 1.0) {
        return ml_exp (x * x) * ml_erfc (x);
    }
    {
        double b = x;
        double c = 1e300, d = 1.0 / b, h = d;
        int n;
        for (n = 1; n < 300; n++) {
            double a = (double) n * 0.5;
            double del;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            del = c * d;
            h *= del;
            if (ml_fabs (del - 1.0) < 1e-16) break;
        }
        return h * 0.5641895835477563; /* CF gives sqrt(pi)*erfcx */
    }
}
/* ---------------------------------------------------------------- Polygamma */

ML_API double ml_polygamma (int n, double x) {
    static const double B [] = {0.0, 1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0, 5.0 / 66.0, -691.0 / 2730.0};
    if (n < 0 || n > 100 || ml_isnan (x)) return ml_make_nan ();
    if (n == 0) return ml_digamma (x);
    if (x <= 0.0) return ml_make_nan ();
    {
        double sum = 0.0, xx = x;
        int j;
        /* recurrence shift x up to >= 12 */
        while (xx < 12.0) {
            /* add (-1)^(n+1) n! / xx^(n+1) */
            {
                double pw = ml_pow (xx, (double) n + 1.0);
                double fact = ml_gamma_new ((double) n + 1.0);
                double sgn = (n % 2) ? 1.0 : -1.0; /* (-1)^(n+1) */
                sum += sgn * fact / pw;
            }
            xx += 1.0;
        }
        if (n >= 8) {
            double tail = 0.0, term;
            double fact_n = ml_gamma_new ((double) n + 1.0);
            double power = (double) n + 1.0;
            for (j = 0; j < 100000; j++) {
                term = 1.0 / ml_pow (xx + (double) j, power);
                tail += term;
                if (term < 1e-18 * tail) break;
            }
            if (j == 100000) return ml_make_nan ();
            {
                double sgn = (n % 2) ? 1.0 : -1.0;
                return sum + sgn * fact_n * tail;
            }
        }
        /* asymptotic */
        {
            double aw = 0.0;
            double fact_nm1 = ml_gamma_new ((double) n);     /* (n-1)! */
            double fact_n = ml_gamma_new ((double) n + 1.0); /* n! */
            double base = ml_pow (xx, (double) n);
            aw = fact_nm1 / base + fact_n / (2.0 * base * xx);
            for (j = 1; j < 7; j++) {
                int k = 2 * j;
                double num = ml_gamma_new ((double) k + (double) n);
                double den = ml_gamma_new ((double) k + 1.0) * ml_pow (xx, (double) (k + n));
                aw += B [j] * num / den;
            }
            {
                double sgn = (n % 2) ? 1.0 : -1.0; /* (-1)^(n+1) */
                sum += sgn * aw;
            }
        }
        return sum;
    }
}
ML_API double ml_trigamma (double x) { return ml_polygamma (1, x); }
/* --------------------------------------------------- Orthogonal polynomials */

ML_API double ml_assoc_legendre (int l, int m, double x) {
    int mm;
    double pmm, plm1, pl;
    int i;
    if (m == INT_MIN || l < 0 || l > SP_MAX_ORDER) return ml_make_nan ();
    mm = m < 0 ? -m : m;
    if (mm > l || ml_fabs (x) > 1.0) return ml_make_nan ();
    /* P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2} */
    pmm = 1.0;
    if (mm > 0) {
        double somx2 = ml_sqrt ((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (i = 1; i <= mm; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == mm) pl = pmm;
    else {
        plm1 = x * (2.0 * (double) mm + 1.0) * pmm;
        if (l == mm + 1) pl = plm1;
        else {
            pl = 0.0;
            for (i = mm + 2; i <= l; i++) {
                pl = (x * (2.0 * (double) i - 1.0) * plm1 - ((double) i + (double) mm - 1.0) * pmm) / ((double) i - (double) mm);
                pmm = plm1;
                plm1 = pl;
            }
        }
    }
    if (m < 0) {
        double ratio = ml_gamma_new ((double) (l - mm) + 1.0) / ml_gamma_new ((double) (l + mm) + 1.0);
        double sgn = (mm % 2) ? -1.0 : 1.0;
        pl *= sgn * ratio;
    }
    return pl;
}
ML_API double ml_spherical_harmonic (int l, int m, double theta, double phi) {
    int mm;
    double norm, pl, ct;
    if (m == INT_MIN || l < 0 || l > SP_MAX_ORDER) return ml_make_nan ();
    mm = m < 0 ? -m : m;
    if (mm > l || !ml_isfinite (theta) || !ml_isfinite (phi)) return ml_make_nan ();
    ct = ml_cos (theta);
    pl = ml_assoc_legendre (l, mm, ct);
    norm = ml_sqrt ((2.0 * (double) l + 1.0) / (4.0 * SP_PI) * ml_exp (ml_lgamma ((double) (l - mm) + 1.0) - ml_lgamma ((double) (l + mm) + 1.0)));
    /* Convention: shake off Condon-Shortley (-1)^m for real harmonics. */
    if (mm % 2) pl = -pl;
    if (m > 0) return ml_sqrt (2.0) * norm * pl * ml_cos ((double) mm * phi);
    if (m < 0) return ml_sqrt (2.0) * norm * pl * ml_sin ((double) mm * phi);
    return norm * pl;
}
ML_API double ml_gegenbauer (int n, double alpha, double x) {
    double c0, c1, c2;
    int k;
    if (n < 0 || n > SP_MAX_ORDER || !ml_isfinite (alpha) || !ml_isfinite (x)) return ml_make_nan ();
    if (n == 0) return 1.0;
    if (n == 1) return 2.0 * alpha * x;
    c0 = 1.0;
    c1 = 2.0 * alpha * x;
    for (k = 1; k < n; k++) {
        double dk = (double) k;
        c2 = (2.0 * (dk + alpha) * x * c1 - (dk + 2.0 * alpha - 1.0) * c0) / (dk + 1.0);
        c0 = c1;
        c1 = c2;
    }
    return c1;
}
ML_API double ml_jacobi_p (int n, double a, double b, double x) {
    double p0, p1, p2;
    int k;
    if (n < 0 || n > SP_MAX_ORDER || !ml_isfinite (a) || !ml_isfinite (b) ||
        !ml_isfinite (x) || a <= -1.0 || b <= -1.0) return ml_make_nan ();
    if (n == 0) return 1.0;
    p0 = 1.0;
    p1 = 0.5 * (a - b + (a + b + 2.0) * x);
    if (n == 1) return p1;
    for (k = 1; k < n; k++) {
        double dk = (double) k;
        double n1 = dk + 1.0;
        double A = 2.0 * n1 * (dk + a + b + 1.0) * (2.0 * dk + a + b);
        double B = (2.0 * dk + a + b + 1.0) * ((2.0 * dk + a + b) * (2.0 * dk + a + b + 2.0) * x + a * a - b * b);
        double C = 2.0 * (dk + a) * (dk + b) * (2.0 * dk + a + b + 2.0);
        p2 = (B * p1 - C * p0) / A;
        p0 = p1;
        p1 = p2;
    }
    return p1;
}
ML_API double ml_zernike_r (int n, int m, double rho) {
    int k;
    double sum = 0.0;
    if (n < 0 || n > SP_MAX_ORDER || m < 0 || m > n ||
        ((n - m) % 2) != 0 || !ml_isfinite (rho)) return ml_make_nan ();
    for (k = 0; k <= (n - m) / 2; k++) {
        double num = ml_gamma_new ((double) (n - k) + 1.0);
        double den = ml_gamma_new ((double) k + 1.0) * ml_gamma_new ((double) ((n + m) / 2 - k) + 1.0) * ml_gamma_new ((double) ((n - m) / 2 - k) + 1.0);
        double sgn = (k % 2) ? -1.0 : 1.0;
        sum += sgn * num / den * ml_pow (rho, (double) (n - 2 * k));
    }
    return sum;
}
/* ---------------------------------------------------- Polylog / Lerch */

ML_API double ml_polylog (int s, double x) {
    double sum = 0.0, xp = 1.0;
    int k, converged = 0;
    if (ml_isnan (x)) return ml_make_nan ();
    if (s < 1) return ml_make_nan ();
    if (x == 1.0) return (s == 1) ? ml_make_inf (0) : ml_zeta ((double) s);
    if (ml_fabs (x) > 1.0) return ml_make_nan ();
    if (x == 0.0) return x;
    if (x == -1.0) {
        if (s == 1) return -ml_log (2.0);
        return (ml_exp2 (1.0 - (double) s) - 1.0) * ml_zeta ((double) s);
    }
    for (k = 1; k < 100000; k++) {
        double term;
        xp *= x;
        term = xp / ml_pow ((double) k, (double) s);
        sum += term;
        if (ml_fabs (term) < 1e-18 * (ml_fabs (sum) + 1e-300)) {
            converged = 1;
            break;
        }
    }
    return converged ? sum : ml_make_nan ();
}
ML_API double ml_lerch_phi (double z, double s, double a) {
    double sum = 0.0, zp = 1.0;
    int k, converged = 0;
    if (ml_isnan (z) || ml_isnan (s) || ml_isnan (a)) return ml_make_nan ();
    if (ml_fabs (z) >= 1.0 || a <= 0.0) return ml_make_nan ();
    if (z == 0.0) return ml_pow (a, -s);
    for (k = 0; k < 200000; k++) {
        double term;
        if (k > 0) zp *= z;
        term = zp / ml_pow ((double) k + a, s);
        sum += term;
        if (ml_fabs (term) < 1e-18 * (ml_fabs (sum) + 1e-300)) {
            converged = 1;
            break;
        }
    }
    return converged ? sum : ml_make_nan ();
}
/* ------------------------------------------------------ Jacobi elliptic */

ML_API double ml_ellip_pi (double n, double k) {
    double m = k * k;
    if (!ml_isfinite (n) || !ml_isfinite (k) || k < 0.0 || k >= 1.0 || n >= 1.0)
        return ml_make_nan ();
    return ml_carlson_rf (0.0, 1.0 - m, 1.0) + (n / 3.0) * ml_carlson_rj (0.0, 1.0 - m, 1.0, 1.0 - n);
}
static int jacobi_solve (double u, double m, double *sn, double *cn, double *dn) {
    double k, K, per, uu, phi;
    int sn_sign = 1, cn_sign = 1, i;
    if (!ml_isfinite (u) || ml_isnan (m) || m < 0.0 || m > 1.0) return -1;
    if (m == 0.0) {
        *sn = ml_sin (u);
        *cn = ml_cos (u);
        *dn = 1.0;
        return 0;
    }
    if (m == 1.0) {
        double s = ml_tanh (u);
        *sn = s;
        *cn = 1.0 / ml_cosh (u);
        *dn = 1.0 / ml_cosh (u);
        return 0;
    }
    k = ml_sqrt (m);
    K = ml_ellip_k (k);
    per = 4.0 * K;
    uu = ml_fmod (u, per);
    if (uu < 0.0) uu += per;
    if (uu > 2.0 * K) {
        uu -= 2.0 * K;
        sn_sign = -1;
        cn_sign = -1;
    }
    if (uu > K) {
        uu = 2.0 * K - uu;
        cn_sign = -cn_sign;
    }
    if (uu >= K) phi = 0.5 * SP_PI;
    else {
        double lo = 0.0, hi = 0.5 * SP_PI;
        phi = 0.25 * SP_PI;
        for (i = 0; i < 80; i++) {
            double s = ml_sin (phi), c = ml_cos (phi);
            double F = s * ml_carlson_rf (c * c, 1.0 - m * s * s, 1.0);
            if (F < uu) lo = phi;
            else hi = phi;
            phi = 0.5 * (lo + hi);
            if (hi - lo < 2e-16) break;
        }
    }
    *sn = (double) sn_sign * ml_sin (phi);
    *cn = (double) cn_sign * ml_cos (phi);
    *dn = ml_sqrt (1.0 - m * ml_sin (phi) * ml_sin (phi));
    return 0;
}
ML_API double ml_jacobi_sn (double u, double m) {
    double s, c, d;
    if (jacobi_solve (u, m, &s, &c, &d) != 0) return ml_make_nan ();
    return s;
}
ML_API double ml_jacobi_cn (double u, double m) {
    double s, c, d;
    if (jacobi_solve (u, m, &s, &c, &d) != 0) return ml_make_nan ();
    return c;
}
ML_API double ml_jacobi_dn (double u, double m) {
    double s, c, d;
    if (jacobi_solve (u, m, &s, &c, &d) != 0) return ml_make_nan ();
    return d;
}
/* ---------------------------------------------------------------- Struve */

static long double struve_h0_integrand (long double theta, long double x) {
    return __builtin_sinl (x * __builtin_cosl (theta));
}
static long double struve_h0_rec (long double a, long double b,
                                  long double fa, long double fm, long double fb,
                                  long double whole, long double x,
                                  long double tol, int depth) {
    long double m = (a + b) * 0.5L;
    long double lm = (a + m) * 0.5L, rm = (m + b) * 0.5L;
    long double flm = struve_h0_integrand (lm, x), frm = struve_h0_integrand (rm, x);
    long double left = (m - a) / 6.0L * (fa + 4.0L * flm + fm);
    long double right = (b - m) / 6.0L * (fm + 4.0L * frm + fb);
    long double both = left + right;
    if (depth <= 0 || __builtin_fabsl (both - whole) <= 15.0L * tol)
        return both + (both - whole) / 15.0L;
    return struve_h0_rec (a, m, fa, flm, fm, left, x, tol * 0.5L, depth - 1) + struve_h0_rec (m, b, fm, frm, fb, right, x, tol * 0.5L, depth - 1);
}
ML_API double ml_struve_h0 (double x) {
    const long double half_pi = 1.57079632679489661923132169163975144L;
    long double xx, fa, fm, fb, whole, integral;
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return ml_copysign (0.0, x);
    if (ml_fabs (x) > 10000.0) return ml_make_nan ();
    xx = (long double) x;
    fa = struve_h0_integrand (0.0L, xx);
    fm = struve_h0_integrand (half_pi * 0.5L, xx);
    fb = 0.0L;
    whole = half_pi / 6.0L * (fa + 4.0L * fm + fb);
    integral = struve_h0_rec (0.0L, half_pi, fa, fm, fb, whole, xx,
                              1e-18L, 26);
    return (double) (2.0L * integral / 3.14159265358979323846264338327950288L);
}
