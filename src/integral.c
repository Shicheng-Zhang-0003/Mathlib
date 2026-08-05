#include "ml_compiler.h"
#include "ml_integral.h"
#include "ml_trig.h"

/*
* v11S contract:
* ml_factorial_float(x) means x! = Gamma(x + 1).
*/
ML_API double ml_factorial_float(double x) {
    if (ml_isnan(x)) return x;
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x == 0.0) return 1.0;
    return ml_gamma_new(x + 1.0);
}

/*
* Educational / experimental numerical integrator.
* Retained for compatibility.
*/
ML_API double ml_integral_traditional(double a, double b, double exponent, double additive, double d) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(exponent) ||
        ml_isnan(additive) || ml_isnan(d)) {
        return ml_make_nan();
    }
    if (ml_isinf(a) || ml_isinf(b) || ml_isinf(exponent) ||
        ml_isinf(additive) || ml_isinf(d)) {
        return ml_make_nan();
    }
    if (d == 0.0) {
        return ml_make_nan();
    }
    if ((d > 0.0 && a >= b) || (d < 0.0 && a <= b)) {
        return 0.0;
    }
    double result = 0.0;
    double x = a;
    const int max_steps = 10000000;
    for (int step = 0; step < max_steps; step++) {
        if ((d > 0.0 && x >= b) || (d < 0.0 && x <= b)) {
            return result;
        }
        double term = ml_pow(x, exponent) + additive;
        if (ml_isnan(term)) {
            return ml_make_nan();
        }
        result += term * d;
        double next_x = x + d;
        if (next_x == x) {
            return ml_make_nan();
        }
        x = next_x;
    }
    return ml_make_nan();
}

/* MATHLIB_V12A1_GAMMA_STIRLING_V4 */
/*
* Gamma / log-gamma via double-double Stirling expansion.
*
* Changes vs V3:
*   - Recurrence depth adapts: x < 1 uses 12 steps, 1 <= x < 8 uses 8.
*   - Reflection uses compensated sin(pi*x) via ml_sinpi().
*   - Stirling correction carries 12 terms (B_2 .. B_24).
*   - Half-integer reflection uses |sin(pi*x)| = 1 exactly.
*   - Exact integer shortcuts for n <= 23 retained.
*/

#ifndef ML_LN2_HI
#define ML_LN2_HI 6.93147180369123816490e-01
#endif
#ifndef ML_LN2_LO
#define ML_LN2_LO 1.90821492927058500170e-10
#endif
#define ML_HALF_LOG_2PI 0.91893853320467274178
#define ML_GAMMA_OVERFLOW 171.6243769563027
#define ML_GAMMA_EXP_OVERFLOW 709.782712893384
#define ML_GAMMA_EXP_UNDERFLOW (-745.133219101941)

/* Double-double pi */
#define ML_PI_HI_D 0x1.921fb54442d18p+1
#define ML_PI_LO_D 0x1.1a62633145c07p-53

/* ---- DD primitives ---- */
typedef struct { double hi; double lo; } ml_dd_t;

static inline ml_dd_t ml_dd_from_d(double a) {
    ml_dd_t r; r.hi = a; r.lo = 0.0; return r;
}
static inline ml_dd_t ml_dd_two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    ml_dd_t r;
    r.hi = s;
    r.lo = (a - (s - v)) + (b - v);
    return r;
}
static inline ml_dd_t ml_dd_renorm(double hi, double lo) {
    return ml_dd_two_sum(hi, lo);
}
static inline ml_dd_t ml_dd_add(ml_dd_t a, ml_dd_t b) {
    ml_dd_t s = ml_dd_two_sum(a.hi, b.hi);
    return ml_dd_renorm(s.hi, s.lo + a.lo + b.lo);
}
static inline ml_dd_t ml_dd_add_d(ml_dd_t a, double b) {
    return ml_dd_add(a, ml_dd_from_d(b));
}
static inline ml_dd_t ml_dd_sub(ml_dd_t a, ml_dd_t b) {
    ml_dd_t nb; nb.hi = -b.hi; nb.lo = -b.lo;
    return ml_dd_add(a, nb);
}
static inline ml_dd_t ml_dd_mul_d(ml_dd_t a, double b) {
    double p = a.hi * b;
    double e = ML_FMA(a.hi, b, -p) + a.lo * b;
    return ml_dd_renorm(p, e);
}
static inline ml_dd_t ml_dd_mul(ml_dd_t a, ml_dd_t b) {
    double p = a.hi * b.hi;
    double e = ML_FMA(a.hi, b.hi, -p) + (a.hi * b.lo + a.lo * b.hi);
    return ml_dd_renorm(p, e);
}

/* ---- DD log ---- */
static ml_dd_t ml_log_dd(double x) {
    if (ml_isnan(x) || x <= 0.0) return ml_dd_from_d(ml_make_nan());
    if (ml_isinf(x))             return ml_dd_from_d(x);
    if (x == 1.0)                return ml_dd_from_d(0.0);

    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + (double)adjust);
    e -= adjust;

    double num = m - 1.0;
    ml_dd_t den = ml_dd_two_sum(m, 1.0);
    double q  = num / den.hi;
    double r  = ML_FMA(-q, den.hi, num);
    r         = ML_FMA(-q, den.lo, r);
    double q2 = r / den.hi;
    ml_dd_t z = ml_dd_renorm(q, q2);
    ml_dd_t z2 = ml_dd_mul(z, z);

    static const double lc[11] = {
        2.0,
        0.6666666666666666,
        0.4,
        0.2857142857142857,
        0.2222222222222222,
        0.18181818181818182,
        0.15384615384615385,
        0.13333333333333333,
        0.11764705882352941,
        0.10526315789473684,
        0.09523809523809523
    };
    ml_dd_t p = ml_dd_from_d(lc[10]);
    for (int i = 9; i >= 0; i--)
        p = ml_dd_add(ml_dd_mul(p, z2), ml_dd_from_d(lc[i]));
    ml_dd_t lm = ml_dd_mul(z, p);

    double ed   = (double)e;
    double ehi  = ed * ML_LN2_HI;
    double elo  = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
    ml_dd_t eln2 = ml_dd_renorm(ehi, elo);
    return ml_dd_add(lm, eln2);
}

/* DD log(pi) */
static ml_dd_t ml_log_pi_dd(void) {
    ml_dd_t lp = ml_log_dd(ML_PI_HI_D);
    return ml_dd_add_d(lp, ML_PI_LO_D / ML_PI_HI_D);
}

/* ---- Stirling expansion, 12 terms (B_2 .. B_24) ---- */
static ml_dd_t ml_stirling_lgamma_dd(double x) {
    if (!ml_isfinite(x))  return ml_dd_from_d(x);
    if (x <= 0.0)         return ml_dd_from_d(ml_make_nan());

    static const double sc[12] = {
         1.0 / 12.0,
        -1.0 / 360.0,
         1.0 / 1260.0,
        -1.0 / 1680.0,
         1.0 / 1188.0,
        -691.0 / 360360.0,
         1.0 / 156.0,
        -3617.0 / 122400.0,
         43867.0 / 244188.0,
        -174611.0 / 125400.0,
         854513.0 / 138600.0,
        -236364091.0 / 6547290.0
    };
    double invx  = 1.0 / x;
    double invx2 = invx * invx;
    double corr  = invx * sc[0];
    double p     = invx * invx2;
    for (int i = 1; i < 12; i++) {
        corr += p * sc[i];
        p    *= invx2;
    }

    ml_dd_t lx   = ml_log_dd(x);
    ml_dd_t w    = ml_dd_from_d(x - 0.5);
    ml_dd_t prod = ml_dd_mul(w, lx);
    ml_dd_t L    = ml_dd_sub(prod, ml_dd_from_d(x));
    L = ml_dd_add_d(L, ML_HALF_LOG_2PI);
    L = ml_dd_add_d(L, corr);
    return L;
}

/* ---- Positive-domain lgamma with adaptive recurrence ---- */
static ml_dd_t ml_lgamma_positive_dd(double x) {
    if (x >= 8.0)
        return ml_stirling_lgamma_dd(x);

    /*
    * Adaptive depth:
    *   x < 1   -> 12 steps  (avoids cancellation for tiny x)
    *   x >= 1  ->  8 steps  (sufficient for x in [1,8))
    */
    int depth = (x < 1.0) ? 12 : 8;
    ml_dd_t L = ml_stirling_lgamma_dd(x + (double)depth);
    for (int k = 0; k < depth; k++) {
        ml_dd_t lv = ml_log_dd(x + (double)k);
        L = ml_dd_sub(L, lv);
    }
    return L;
}

/* ---- exp(hi+lo) ---- */
static double ml_exp_dd(ml_dd_t L) {
    if (L.hi > ML_GAMMA_EXP_OVERFLOW)  return ml_make_inf(0);
    if (L.hi < ML_GAMMA_EXP_UNDERFLOW) return 0.0;
    double g = ml_exp(L.hi);
    if (!ml_isfinite(g) || g == 0.0) return g;
    return ML_FMA(g, L.lo, g);
}

/* ---- Half-integer helpers ---- */
static int ml_is_half_integer(double x) {
    if (!ml_isfinite(x))     return 0;
    if (x == ml_round(x))    return 0;
    double t = 2.0 * x;
    if (!ml_isfinite(t))     return 0;
    return t == ml_round(t);
}
static double ml_half_sin_sign(double x) {
    double nd = ml_round(x - 0.5);
    if (!ml_isfinite(nd) || ml_fabs(nd) >= 9007199254740992.0)
        return 1.0;
    long long n = (long long)nd;
    if ((n % 2LL) != 0LL) return -1.0;
    return 1.0;
}

/* ---- Compensated sin(pi*x) ---- */
static double ml_sinpi(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (x == 0.0) return x;
    double arg_hi = x * ML_PI_HI_D;
    double arg_lo = ML_FMA(x, ML_PI_HI_D, -arg_hi) + x * ML_PI_LO_D;
    double arg    = arg_hi + arg_lo;
    if (!ml_isfinite(arg))
        return ml_sin(ML_PI * x);
    return ml_sin(arg);
}

/* ---- Exact small factorials ---- */
static double ml_factorial_exact_small(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++)
        f *= (double)k;
    return f;
}

/* ---- Public APIs ---- */
ML_API double ml_lgamma(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x <= 0.0 && x == ml_round(x)) return ml_make_inf(0);
    if (x == 1.0 || x == 2.0) return 0.0;

    if (x > 0.0) {
        if (x == ml_round(x) && x <= 23.0)
            return ml_log(ml_factorial_exact_small((int)x));
        ml_dd_t L = ml_lgamma_positive_dd(x);
        return L.hi + L.lo;
    }

    /* Reflection: log|Gamma(x)| = log(pi) - log|sin(pi x)| - lgamma(1-x) */
    double s = ml_sinpi(x);
    if (s == 0.0) return ml_make_inf(0);
    double absin = ml_fabs(s);
    if (ml_is_half_integer(x)) absin = 1.0;
    if (absin == 0.0) return ml_make_inf(0);

    ml_dd_t logterm = ml_log_pi_dd();
    if (absin != 1.0)
        logterm = ml_dd_sub(logterm, ml_log_dd(absin));
    ml_dd_t Lpos = ml_lgamma_positive_dd(1.0 - x);
    ml_dd_t r    = ml_dd_sub(logterm, Lpos);
    return r.hi + r.lo;
}

ML_API double ml_gamma_new(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : ml_make_nan();
    if (x <= 0.0 && x == ml_round(x)) return ml_make_nan();

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) return ml_make_inf(0);
        if (x == ml_round(x) && x <= 23.0)
            return ml_factorial_exact_small((int)x);
        ml_dd_t L = ml_lgamma_positive_dd(x);
        return ml_exp_dd(L);
    }

    /* Reflection: Gamma(x) = pi / (sin(pi x) * Gamma(1-x)) */
    double s = ml_sinpi(x);
    if (s == 0.0) return ml_make_nan();
    double sin_use = s;
    if (ml_is_half_integer(x))
        sin_use = ml_half_sin_sign(x);
    if (sin_use == 0.0) return ml_make_nan();

    ml_dd_t Lpos = ml_lgamma_positive_dd(1.0 - x);
    double  G    = ml_exp_dd(Lpos);
    if (ml_isinf(G))  return ml_copysign(0.0, sin_use);
    if (G == 0.0)     return (sin_use < 0.0) ? -ml_make_inf(0) : ml_make_inf(0);
    return ML_PI_HI_D / (sin_use * G);
}
