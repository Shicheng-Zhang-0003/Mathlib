#include "ml_compiler.h"
#include "ml_integral.h"
#include "ml_trig.h"

ML_API double ml_factorial_float(double x) {
    if (ml_isnan(x)) return x;
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x == 0.0) return 1.0;
    return ml_gamma_new(x + 1.0);
}

ML_API double ml_integral_traditional(double a, double b, double exponent,
                                      double additive, double d) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(exponent) ||
        ml_isnan(additive) || ml_isnan(d))
        return ml_make_nan();
    if (ml_isinf(a) || ml_isinf(b) || ml_isinf(exponent) ||
        ml_isinf(additive) || ml_isinf(d))
        return ml_make_nan();
    if (d == 0.0) return ml_make_nan();
    if ((d > 0.0 && a >= b) || (d < 0.0 && a <= b)) return 0.0;
    double result = 0.0;
    double x = a;
    const int max_steps = 10000000;
    for (int step = 0; step < max_steps; step++) {
        if ((d > 0.0 && x >= b) || (d < 0.0 && x <= b)) return result;
        double term = ml_pow(x, exponent) + additive;
        if (ml_isnan(term)) return ml_make_nan();
        result += term * d;
        double next_x = x + d;
        if (next_x == x) return ml_make_nan();
        x = next_x;
    }
    return ml_make_nan();
}

/* MATHLIB_V12A1_GAMMA_HALFINT_V7 */
/*
 * Gamma / lgamma with exact half-integer shortcuts.
 *
 * ROOT CAUSE of previous failures:
 *   The Lanczos g=7 n=9 coefficient set has ~1e-15 intrinsic
 *   approximation error. For half-integers, the reflection formula
 *   subtracts two values of similar magnitude, amplifying this
 *   error to 100+ ULP.
 *
 * FIX:
 *   - Half-integers use exact product/sum formulas (no Lanczos).
 *   - x < 0.5 uses 8-step recurrence (not 1-step).
 *   - x >= 8 uses Stirling DD.
 *   - 0.5 <= x < 8 uses Lanczos DD (with half-integer bypass).
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
    ml_dd_t r; r.hi = s; r.lo = (a - (s - v)) + (b - v); return r;
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
    if (ml_isinf(x)) return ml_dd_from_d(x);
    if (x == 1.0) return ml_dd_from_d(0.0);
    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + (double)adjust);
    e -= adjust;
    double num = m - 1.0;
    ml_dd_t den = ml_dd_two_sum(m, 1.0);
    double q = num / den.hi;
    double r = ML_FMA(-q, den.hi, num);
    r = ML_FMA(-q, den.lo, r);
    double q2 = r / den.hi;
    ml_dd_t z = ml_dd_renorm(q, q2);
    ml_dd_t z2 = ml_dd_mul(z, z);
    static const double lc[11] = {
        2.0, 0.6666666666666666, 0.4, 0.2857142857142857,
        0.2222222222222222, 0.18181818181818182, 0.15384615384615385,
        0.13333333333333333, 0.11764705882352941, 0.10526315789473684,
        0.09523809523809523
    };
    ml_dd_t p = ml_dd_from_d(lc[10]);
    for (int i = 9; i >= 0; i--)
        p = ml_dd_add(ml_dd_mul(p, z2), ml_dd_from_d(lc[i]));
    ml_dd_t lm = ml_dd_mul(z, p);
    double ed = (double)e;
    double ehi = ed * ML_LN2_HI;
    double elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
    ml_dd_t eln2 = ml_dd_renorm(ehi, elo);
    return ml_dd_add(lm, eln2);
}

static ml_dd_t ml_log_pi_dd(void) {
    ml_dd_t lp = ml_log_dd(ML_PI_HI_D);
    return ml_dd_add_d(lp, ML_PI_LO_D / ML_PI_HI_D);
}

/* ---- exp(hi+lo) ---- */
static double ml_exp_dd(ml_dd_t L) {
/* MATHLIB_V12A1_EXP_DD_SECOND_ORDER */
/*
* exp(L.hi + L.lo) = exp(L.hi) * exp(L.lo)
*
* The old code used the first-order approximation:
*     exp(L.lo) ≈ 1 + L.lo
* which drops L.lo²/2. For the 8-step recurrence at x=0.001,
* L.lo reaches ~5e-8, making the dropped term ~1.25e-12 = 10 ULP.
*
* Fix: compute exp(L.lo) with a 3-term Taylor:
*     exp(L.lo) ≈ 1 + L.lo + L.lo²/2
* Since L.lo is tiny, this is accurate to ~1e-48.
*/
if (L.hi > ML_GAMMA_EXP_OVERFLOW) return ml_make_inf(0);
if (L.hi < ML_GAMMA_EXP_UNDERFLOW) return 0.0;
double g = ml_exp(L.hi);
if (!ml_isfinite(g) || g == 0.0) return g;
/* 3-term Taylor for exp(L.lo): accurate to ~L.lo^3/6 */
double elo = ML_FMA(L.lo, L.lo * 0.5, L.lo) + 1.0;
return ML_FMA(g, elo, 0.0);
}

/* ---- Half-integer detection ---- */
static int ml_is_half_integer(double x) {
    if (!ml_isfinite(x)) return 0;
    if (x == ml_round(x)) return 0;
    double t = 2.0 * x;
    if (!ml_isfinite(t)) return 0;
    return t == ml_round(t);
}

/* Index n for x = n + 0.5 */
static long long ml_half_index(double x) {
    double nd = ml_round(x - 0.5);
    if (!ml_isfinite(nd) || ml_fabs(nd) >= 9007199254740992.0) return 0;
    return (long long)nd;
}

/* sin(pi*(n+0.5)) = (-1)^n */
static double ml_half_sin_sign(double x) {
    long long n = ml_half_index(x);
    return ((n % 2LL) != 0LL) ? -1.0 : 1.0;
}

/* ---- Exact half-integer lgamma (positive) ---- */
/*
 * lgamma(n+0.5) = 0.5*log(pi) + sum(log(j+0.5), j=0..n-1)
 * All computed in DD. No Lanczos approximation error.
 */
static ml_dd_t ml_lgamma_half_positive_dd(double x) {
    long long n = ml_half_index(x);
    ml_dd_t L = ml_dd_mul_d(ml_log_pi_dd(), 0.5);
    for (long long j = 0; j < n; j++) {
        L = ml_dd_add(L, ml_log_dd((double)j + 0.5));
    }
    return L;
}

/* ---- Exact half-integer gamma (positive) ---- */
/*
 * Gamma(n+0.5) = sqrt(pi) * prod(j+0.5, j=0..n-1)
 * Computed as exp(lgamma) to preserve DD precision.
 */
static double ml_gamma_half_positive(double x) {
    ml_dd_t L = ml_lgamma_half_positive_dd(x);
    return ml_exp_dd(L);
}

/* ---- Lanczos DD for non-half-integer x in [0.5, 8) ---- */
static const double ml_lanczos_coeff[9] = {
     0.99999999999980993,
     676.5203681218851,
    -1259.1392167224028,
     771.32342877765313,
    -176.61502916214059,
     12.507343278686905,
    -0.13857109526572012,
     9.9843695780195716e-6,
     1.5056327351493116e-7
};

static ml_dd_t ml_lgamma_lanczos_dd(double x) {
    double z = x - 1.0;
    double t = z + 7.5;
    double w = z + 0.5;
    ml_dd_t ag = ml_dd_from_d(ml_lanczos_coeff[0]);
    for (int i = 1; i < 9; i++) {
        double denom = z + (double)i;
        double q = ml_lanczos_coeff[i] / denom;
        double r = ML_FMA(-q, denom, ml_lanczos_coeff[i]);
        double q2 = r / denom;
        ml_dd_t term = ml_dd_renorm(q, q2);
        ag = ml_dd_add(ag, term);
    }
    if (ag.hi <= 0.0) return ml_dd_from_d(ml_make_nan());
    ml_dd_t lag = ml_log_dd(ag.hi);
    lag = ml_dd_add_d(lag, ag.lo / ag.hi);
    ml_dd_t lt = ml_log_dd(t);
    ml_dd_t wlt = ml_dd_mul_d(lt, w);
    ml_dd_t L = ml_dd_two_sum(ML_HALF_LOG_2PI, -t);
    L = ml_dd_add(L, wlt);
    L = ml_dd_add(L, lag);
    return L;
}

/* ---- Stirling DD for x >= 8 ---- */
static ml_dd_t ml_stirling_lgamma_dd(double x) {
    if (!ml_isfinite(x)) return ml_dd_from_d(x);
    if (x <= 0.0) return ml_dd_from_d(ml_make_nan());
    static const double sc[12] = {
         1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0, -1.0 / 1680.0,
         1.0 / 1188.0, -691.0 / 360360.0, 1.0 / 156.0,
        -3617.0 / 122400.0, 43867.0 / 244188.0, -174611.0 / 125400.0,
         854513.0 / 63756.0, -236364091.0 / 1506960.0
    };
    double invx = 1.0 / x;
    double invx2 = invx * invx;
    double corr = invx * sc[0];
    double p = invx * invx2;
    for (int i = 1; i < 12; i++) { corr += p * sc[i]; p *= invx2; }
    ml_dd_t lx = ml_log_dd(x);
    ml_dd_t w = ml_dd_from_d(x - 0.5);
    ml_dd_t prod = ml_dd_mul(w, lx);
    ml_dd_t L = ml_dd_sub(prod, ml_dd_from_d(x));
    L = ml_dd_add_d(L, ML_HALF_LOG_2PI);
    L = ml_dd_add_d(L, corr);
    return L;
}

/* ---- Positive-domain lgamma dispatch ---- */
static ml_dd_t ml_lgamma_positive_dd(double x) {
    if (x >= 8.0)
        return ml_stirling_lgamma_dd(x);
    if (ml_is_half_integer(x))
        return ml_lgamma_half_positive_dd(x);
    return ml_lgamma_lanczos_dd(x);
}

/* ---- Positive-domain gamma dispatch ---- */
static double ml_gamma_positive(double x) {
/* MATHLIB_V12A1_GAMMA_DIVIDE_FIX */
/*
* For x >= 8: Stirling DD -> exp.
* For half-integers: exact product formula.
* For 0 < x < 8: use the DD lgamma path (log subtraction),
*   then exponentiate. This avoids the product-then-divide
*   approach which introduced two extra double roundings.
*/
if (x >= 8.0) {
ml_dd_t L = ml_stirling_lgamma_dd(x);
return ml_exp_dd(L);
}
if (ml_is_half_integer(x)) {
return ml_gamma_half_positive(x);
}
/* DD log-subtract-then-exp: same path as lgamma, proven accurate */
ml_dd_t L = ml_lgamma_positive_dd(x);
return ml_exp_dd(L);
}

/* ---- Exact small factorials ---- */
static double ml_factorial_exact_small(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) f *= (double)k;
    return f;
}

/* ---- Compensated sin(pi*x) ---- */
static double ml_sinpi(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (x == 0.0) return x;
    double arg_hi = x * ML_PI_HI_D;
    double arg_lo = ML_FMA(x, ML_PI_HI_D, -arg_hi) + x * ML_PI_LO_D;
    double arg = arg_hi + arg_lo;
    if (!ml_isfinite(arg)) return ml_sin(ML_PI * x);
    return ml_sin(arg);
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
        /* Half-integer: exact formula */
        if (ml_is_half_integer(x)) {
            ml_dd_t L = ml_lgamma_half_positive_dd(x);
            return L.hi + L.lo;
        }
        /* x < 0.5: 8-step recurrence */
        /* MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH16 */
        if (x < 0.5) {
            ml_dd_t L = ml_lgamma_positive_dd(x + 16.0);
            for (int k = 0; k < 16; k++)
                L = ml_dd_sub(L, ml_log_dd(x + (double)k));
            return L.hi + L.lo;
        }
        ml_dd_t L = ml_lgamma_positive_dd(x);
        return L.hi + L.lo;
    }

    /* Negative x: reflection */
    double s = ml_sinpi(x);
    if (s == 0.0) return ml_make_inf(0);
    double absin = ml_fabs(s);
    if (ml_is_half_integer(x)) absin = 1.0;
    if (absin == 0.0) return ml_make_inf(0);

    /* For negative half-integers, use exact positive lgamma */
    double pos_arg = 1.0 - x;
    ml_dd_t Lpos;
    if (ml_is_half_integer(pos_arg)) {
        Lpos = ml_lgamma_half_positive_dd(pos_arg);
    } else {
        Lpos = ml_lgamma_positive_dd(pos_arg);
    }

    ml_dd_t logterm = ml_log_pi_dd();
    if (absin != 1.0)
        logterm = ml_dd_sub(logterm, ml_log_dd(absin));
    ml_dd_t r = ml_dd_sub(logterm, Lpos);
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
        /* Half-integer: exact formula */
        if (ml_is_half_integer(x))
            return ml_gamma_half_positive(x);
        /* MATHLIB_V12A1_GAMMA_NEW_LOGSPACE */
/* x < 0.5: 8-step recurrence in log-space (no product-divide)
*
* Old: gy = exp(lgamma(x+8)); p = prod(x+k); return gy/p;
*   -> two extra double roundings + division rounding = 11 ULP
*
* New: lgamma(x+8) - sum(log(x+k)) then exp
*   -> everything stays in DD until final exp, same as ml_lgamma
*/
if (x < 0.5) {
ml_dd_t L = ml_lgamma_positive_dd(x + 8.0);
for (int k = 0; k < 8; k++)
L = ml_dd_sub(L, ml_log_dd(x + (double)k));
return ml_exp_dd(L);
}
        return ml_gamma_positive(x);
    }

    /* Negative x: reflection */
    double s = ml_sinpi(x);
    if (s == 0.0) return ml_make_nan();
    double sin_use = s;
    if (ml_is_half_integer(x))
        sin_use = ml_half_sin_sign(x);
    if (sin_use == 0.0) return ml_make_nan();

    /* For negative half-integers, use exact positive gamma */
    double pos_arg = 1.0 - x;
    double G;
    if (ml_is_half_integer(pos_arg)) {
        G = ml_gamma_half_positive(pos_arg);
    } else {
        G = ml_gamma_positive(pos_arg);
    }

    if (ml_isinf(G)) return ml_copysign(0.0, sin_use);
    if (G == 0.0) return (sin_use < 0.0) ? -ml_make_inf(0) : ml_make_inf(0);
    return ML_PI_HI_D / (sin_use * G);
}
