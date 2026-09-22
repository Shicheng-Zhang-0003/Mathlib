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
        if (!ml_isfinite(term)) return ml_make_nan();
        result += term * d;
        double next_x = x + d;
        if (next_x == x) return ml_make_nan();
        x = next_x;
    }
    return ml_make_nan();
}

/* MATHLIB_V12A1_INTEGRAL_COMMENT_CLEANUP */
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
 *   - x < 0.5 uses 1-step recurrence: lgamma(x) = lgamma(x+1) - log(x).
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

/* ---- Exact (n-1)! for integer gamma: Gamma(n) = (n-1)! ----
 * Named ml_gamma_int_fact to avoid the off-by-one trap: for integer
 * n>=1, Gamma(n)=(n-1)!. Call sites pass the gamma argument n directly. */
static double ml_gamma_int_fact(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) f *= (double)k;
    return f;
}

/* ---- Compensated sin(pi*x): use public ml_sinpi from trig.c
 * (mod-2 reduction + HI/LO compensation). Static removed to avoid
 * clash with ML_API ml_sinpi. ---- */

/* ---- Public APIs ---- */
ML_API double ml_lgamma(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x <= 0.0 && x == ml_round(x)) return ml_make_inf(0);
    if (x == 1.0 || x == 2.0) return 0.0;

    if (x > 0.0) {
        if (x == ml_round(x) && x <= 23.0)
            return ml_log(ml_gamma_int_fact((int)x));
        /* Half-integer: exact formula is O(n); for x>=8 Stirling is
         * 1e-15-accurate and avoids hanging on 1e15+0.5. */
        if (ml_is_half_integer(x)) {
            if (x >= 8.0) {
                ml_dd_t L = ml_stirling_lgamma_dd(x);
                return L.hi + L.lo;
            }
            ml_dd_t L = ml_lgamma_half_positive_dd(x);
            return L.hi + L.lo;
        }
        /* x < 0.5: 8-step recurrence */
        /* MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH16 */
        /* MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE */
/* x < 0.5: 1-step recurrence in log-space.
 * lgamma(x) = lgamma(x+1) - log(x)
 * For x=0.1: lgamma(1.1)~-0.05, log(0.1)~-2.30.
 * This is addition (both positive), not subtraction.
 * Cancellation factor ~1.04 (almost zero).
 */
if (x < 0.5) {
    ml_dd_t L = ml_lgamma_positive_dd(x + 1.0);
    L = ml_dd_sub(L, ml_log_dd(x));
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

    /* For negative half-integers, use exact positive lgamma unless the
     * positive argument is huge (hang guard: Stirling instead). */
    double pos_arg = 1.0 - x;
    ml_dd_t Lpos;
    if (ml_is_half_integer(pos_arg)) {
        if (pos_arg >= 8.0) Lpos = ml_stirling_lgamma_dd(pos_arg);
        else Lpos = ml_lgamma_half_positive_dd(pos_arg);
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
            return ml_gamma_int_fact((int)x);
        /* Half-integer: exact formula is O(n); Stirling for x>=8. */
        if (ml_is_half_integer(x)) {
            if (x >= 8.0) {
                ml_dd_t L = ml_stirling_lgamma_dd(x);
                return ml_exp_dd(L);
            }
            return ml_gamma_half_positive(x);
        }
        /* MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS */
/* x < 0.5: use Lanczos directly via ml_gamma_positive.
 * Same cancellation fix as ml_lgamma above.
 */
/* MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE */
/* x < 0.5: 1-step recurrence in log-space, then exp.
 * Same cancellation fix as ml_lgamma above.
 * Uses ml_exp_dd(L) instead of g/x to avoid division rounding.
 */
if (x < 0.5) {
    ml_dd_t L = ml_lgamma_positive_dd(x + 1.0);
    L = ml_dd_sub(L, ml_log_dd(x));
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

    /* For negative half-integers, use exact positive gamma unless huge. */
    double pos_arg = 1.0 - x;
    double G;
    if (ml_is_half_integer(pos_arg)) {
        if (pos_arg >= 8.0) {
            ml_dd_t L = ml_stirling_lgamma_dd(pos_arg);
            G = ml_exp_dd(L);
        } else {
            G = ml_gamma_half_positive(pos_arg);
        }
    } else {
        G = ml_gamma_positive(pos_arg);
    }

    if (ml_isinf(G)) return ml_copysign(0.0, sin_use);
    if (G == 0.0) return (sin_use < 0.0) ? -ml_make_inf(0) : ml_make_inf(0);
    /* Full pi (HI+LO) for ~0.5 ULP gain over HI-only. Double-double
     * division via reciprocal refinement: q=HI/(s*G), then one
     * correction with LO. */
    {
        double denom = sin_use * G;
        double q = ML_PI_HI_D / denom;
        double r = ML_FMA(-q, denom, ML_PI_HI_D) + ML_PI_LO_D;
        q += r / denom;
        return q;
    }
}

ML_API double ml_digamma(double x) {
    /* psi = Gamma'/Gamma. Recurrence psi(x+1)=psi(x)+1/x to x>=32,
     * then Stirling: ln x - 1/2x - 1/12x^2 + 1/120x^4 - 1/252x^6.
     * Shift to 32 (not 8) so truncation ~1e-15 even without x^-8 term;
     * critical for Y1/K1 series which call digamma ~240x. */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x;
    if (x <= 0.0 && x == ml_round(x)) return ml_make_nan();
    if (x == 0.0) return -ml_make_inf(0);
    if (x < 0.0) {
        /* Reflection psi(1-x)-psi(x)=pi*cot(pi x) avoids ~1e15-step
         * recurrence hang for digamma(-1e15). */
        double y = 1.0 - x;
        double sn = ml_sinpi(x);
        if (sn == 0.0) return ml_make_nan();
        double dy = ml_digamma(y);
        if (!ml_isfinite(dy)) return dy;
        {
            double cot = ml_cospi(x) / sn;
            return dy - ML_PI * cot;
        }
    }
    {
        double r = 0.0;
        double xx = x;
        while (xx < 32.0) {
            r -= 1.0 / xx;
            xx += 1.0;
        }
        {
            double inv = 1.0 / xx;
            double inv2 = inv * inv;
            double s = ml_log(xx) - 0.5 * inv
                - inv2 * (1.0/12.0 - inv2 * (1.0/120.0 - inv2 * (1.0/252.0)));
            return r + s;
        }
    }
}

static double ml_gamma_series_p(double a, double x) {
    /* P(a,x) series for x < a+1: sum x^n/n! normalization. */
    double gln = ml_lgamma(a);
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    double n = 1.0;
    for (int i = 0; i < 1000; i++) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        n += 1.0;
        (void)n;
        if (ml_fabs(del) < ml_fabs(sum) * 1e-16) break;
    }
    return sum * ml_exp(-x + a * ml_log(x) - gln);
}

static double ml_gamma_cf_q(double a, double x) {
    /* Q(a,x) Lentz continued fraction for x >= a+1. */
    static const double FPMIN = 1e-300;
    double gln = ml_lgamma(a);
    double b = x + 1.0 - a;
    double c = 1.0 / FPMIN;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i < 1000; i++) {
        double an = -(double)i * ((double)i - a);
        b += 2.0;
        d = an * d + b;
        if (ml_fabs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (ml_fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        {
            double del = d * c;
            h *= del;
            if (ml_fabs(del - 1.0) < 1e-15) break;
        }
    }
    return ml_exp(-x + a * ml_log(x) - gln) * h;
}

ML_API double ml_gamma_p(double a, double x) {
    if (ml_isnan(a) || ml_isnan(x)) return ml_make_nan();
    if (!(a > 0.0) || !ml_isfinite(a)) return ml_make_nan();
    if (x < 0.0 || !ml_isfinite(x)) return ml_make_nan();
    if (x == 0.0) return 0.0;
    if (ml_isinf(x)) return 1.0;
    if (x < a + 1.0) return ml_gamma_series_p(a, x);
    return 1.0 - ml_gamma_cf_q(a, x);
}

ML_API double ml_gamma_q(double a, double x) {
    if (ml_isnan(a) || ml_isnan(x)) return ml_make_nan();
    if (!(a > 0.0) || !ml_isfinite(a)) return ml_make_nan();
    if (x < 0.0 || !ml_isfinite(x)) return ml_make_nan();
    if (x == 0.0) return 1.0;
    if (ml_isinf(x)) return 0.0;
    if (x < a + 1.0) return 1.0 - ml_gamma_series_p(a, x);
    return ml_gamma_cf_q(a, x);
}

ML_API double ml_beta(double a, double b) {
    if (ml_isnan(a) || ml_isnan(b)) return ml_make_nan();
    if (!(a > 0.0) || !(b > 0.0) || !ml_isfinite(a) || !ml_isfinite(b)) {
        return ml_make_nan();
    }
    /* Beta = exp(lgamma(a)+lgamma(b)-lgamma(a+b)), DD-stable. */
    {
        double L = ml_lgamma(a) + ml_lgamma(b) - ml_lgamma(a + b);
        return ml_exp(L);
    }
}

ML_API double ml_bessel_j0(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return 0.0;
    {
        double ax = ml_fabs(x);
        if (ax == 0.0) return 1.0;
        if (ax <= 16.0) {
            /* J0 = sum (-1)^m (x^2/4)^m/(m!^2), recurrence, Kahan. */
            double y = (ax * ax) * 0.25;
            double t = 1.0, sum = 1.0, comp = 0.0;
            for (int m = 1; m < 120; m++) {
                t *= -y / ((double)m * (double)m);
                {
                    double w = t - comp;
                    double tt = sum + w;
                    comp = (tt - sum) - w;
                    sum = tt;
                }
                if (ml_fabs(t) < 1e-18 * ml_fabs(sum)) break;
            }
            return sum;
        }
        /* Asymptotic: sqrt(2/pi x)[cos(x-pi/4)P - sin(x-pi/4)Q],
         * P=1-9/128z^2+..., Q=1/8z-..., z=x^2. 3 terms: ~1e-9 at x=8. */
        {
            double z = ax * ax;
            double p = 1.0 - 9.0 / (128.0 * z) + 3675.0 / (32768.0 * z * z);
            double q = 1.0 / (8.0 * ax) - 75.0 / (1024.0 * ax * z);
            double chi = ax - 0.78539816339744830962;
            double amp = 0.79788456080286535588 / ml_sqrt(ax);
            double r = amp * (ml_cos(chi) * p - ml_sin(chi) * q);
            return r;
        }
    }
}

ML_API double ml_bessel_j1(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return 0.0;
    if (x == 0.0) return x;
    {
        double ax = ml_fabs(x);
        double sgn = (x < 0.0) ? -1.0 : 1.0;
        if (ax <= 16.0) {
            double y = (ax * ax) * 0.25;
            double t = 0.5 * ax;
            double sum = t, comp = 0.0;
            for (int m = 1; m < 120; m++) {
                t *= -y / ((double)m * (double)(m + 1));
                {
                    double w = t - comp;
                    double tt = sum + w;
                    comp = (tt - sum) - w;
                    sum = tt;
                }
                if (ml_fabs(t) < 1e-18 * ml_fabs(sum)) break;
            }
            return sgn * sum;
        }
        {
            double z = ax * ax;
            double p = 1.0 + 15.0 / (128.0 * z) - 4725.0 / (32768.0 * z * z);
            double q = 3.0 / (8.0 * ax) - 105.0 / (1024.0 * ax * z);
            double chi = ax - 2.35619449019234492885;
            double amp = 0.79788456080286535588 / ml_sqrt(ax);
            double r = amp * (ml_cos(chi) * p - ml_sin(chi) * q);
            return sgn * r;
        }
    }
}

/* ---- Carlson symmetric forms (DLMF 19.16, 19.36): clean-room
 * duplication with 5th/7th-order Taylor. No printf, NaN on domain error.
 * errtol 1e-3 -> truncation < 1e-16 (TOMS577 tables). */

static double ml_rf_poly(double e2, double e3) {
    /* DLMF 19.36.1 degree 7. */
    double e22 = e2 * e2, e33 = e3 * e3, e23 = e2 * e3;
    return 1.0 - e2*(1.0/10.0) + e3*(1.0/14.0) + e22*(1.0/24.0)
        - e23*(3.0/44.0) - e22*e2*(5.0/208.0) + e33*(3.0/104.0)
        + e22*e3*(1.0/16.0);
}

ML_API double ml_carlson_rf(double x, double y, double z) {
    if (ml_isnan(x) || ml_isnan(y) || ml_isnan(z)) return ml_make_nan();
    if (x < 0.0 || y < 0.0 || z < 0.0) return ml_make_nan();
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z)) return ml_make_nan();
    if (x + y + z == 0.0) return ml_make_nan();
    {
        double xn = x, yn = y, zn = z;
        for (int it = 0; it < 100; it++) {
            double mu = (xn + yn + zn) / 3.0;
            if (mu == 0.0) return ml_make_nan();
            double dx = 2.0 - (mu + xn) / mu;
            double dy = 2.0 - (mu + yn) / mu;
            double dz = 2.0 - (mu + zn) / mu;
            double ax = ml_fabs(dx), ay = ml_fabs(dy), az = ml_fabs(dz);
            double ep = ax > ay ? ax : ay;
            ep = ep > az ? ep : az;
            if (ep < 1e-3) {
                double e2 = dx * dy - dz * dz;
                /* Carlson: e2 = dx*dy+dy*dz+dz*dx but dx+dy+dz=0 so
                 * dx*dy+dy*dz+dz*dx = dx*dy - dz^2 when dz=-(dx+dy)? Use full. */
                e2 = dx * dy + dy * dz + dz * dx;
                {
                    double e3 = dx * dy * dz;
                    double s = ml_rf_poly(e2, e3);
                    return s / ml_sqrt(mu);
                }
            }
            {
                double sx = ml_sqrt(xn), sy = ml_sqrt(yn), sz = ml_sqrt(zn);
                double lam = sx * (sy + sz) + sy * sz;
                xn = (xn + lam) * 0.25;
                yn = (yn + lam) * 0.25;
                zn = (zn + lam) * 0.25;
            }
        }
        return ml_make_nan();
    }
}

ML_API double ml_carlson_rc(double x, double y) {
    if (ml_isnan(x) || ml_isnan(y)) return ml_make_nan();
    if (x < 0.0 || y <= 0.0) return ml_make_nan();
    if (!ml_isfinite(x) || !ml_isfinite(y)) return ml_make_nan();
    /* RC(x,y) = RF(y,y,x) by symmetry of the integral. */
    return ml_carlson_rf(y, y, x);
}

ML_API double ml_carlson_rd(double x, double y, double z) {
    if (ml_isnan(x) || ml_isnan(y) || ml_isnan(z)) return ml_make_nan();
    if (x < 0.0 || y < 0.0 || z <= 0.0) return ml_make_nan();
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z)) return ml_make_nan();
    {
        double xn = x, yn = y, zn = z;
        double sigma = 0.0, fac = 1.0;
        for (int it = 0; it < 100; it++) {
            double mu = (xn + yn + 3.0 * zn) * 0.2;
            if (mu == 0.0) return ml_make_nan();
            double dx = (mu - xn) / mu, dy = (mu - yn) / mu, dz = (mu - zn) / mu;
            double ax = ml_fabs(dx), ay = ml_fabs(dy), az = ml_fabs(dz);
            double ep = ax > ay ? ax : ay;
            ep = ep > az ? ep : az;
            if (ep < 1e-3) {
                double ea = dx * dy, eb = dz * dz;
                double ec = ea - eb, ed = ea - 6.0 * eb;
                double ef = ed + ec + ec;
                double s1 = ed * (-3.0/14.0 + 0.25 * (9.0/22.0) * ed - 1.5 * (3.0/26.0) * dz * ef);
                double s2 = dz * ((1.0/6.0) * ef + dz * (-(9.0/22.0) * ec + dz * (3.0/26.0) * ea));
                return 3.0 * sigma + fac * (1.0 + s1 + s2) / (mu * ml_sqrt(mu));
            }
            {
                double sx = ml_sqrt(xn), sy = ml_sqrt(yn), sz = ml_sqrt(zn);
                double lam = sx * (sy + sz) + sy * sz;
                sigma += fac / (sz * (zn + lam));
                fac *= 0.25;
                xn = (xn + lam) * 0.25;
                yn = (yn + lam) * 0.25;
                zn = (zn + lam) * 0.25;
            }
        }
        return ml_make_nan();
    }
}

ML_API double ml_carlson_rj(double x, double y, double z, double p) {
    if (ml_isnan(x) || ml_isnan(y) || ml_isnan(z) || ml_isnan(p)) return ml_make_nan();
    if (x < 0.0 || y < 0.0 || z < 0.0 || p <= 0.0) return ml_make_nan();
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z) || !ml_isfinite(p)) {
        return ml_make_nan();
    }
    {
        double xn = x, yn = y, zn = z, pn = p;
        double sigma = 0.0, fac = 1.0;
        for (int it = 0; it < 100; it++) {
            double mu = (xn + yn + zn + pn + pn) * 0.2;
            if (mu == 0.0) return ml_make_nan();
            double dx = (mu - xn) / mu, dy = (mu - yn) / mu;
            double dz = (mu - zn) / mu, dp = (mu - pn) / mu;
            double ep = ml_fabs(dx);
            if (ml_fabs(dy) > ep) ep = ml_fabs(dy);
            if (ml_fabs(dz) > ep) ep = ml_fabs(dz);
            if (ml_fabs(dp) > ep) ep = ml_fabs(dp);
            if (ep < 1e-3) {
                double ea = dx * (dy + dz) + dy * dz;
                double eb = dx * dy * dz;
                double ec = dp * dp;
                double e2 = ea - 3.0 * ec;
                double e3 = eb + 2.0 * dp * (ea - ec);
                double s1 = 1.0 + e2 * (-3.0/14.0 + 0.75 * (3.0/22.0) * e2 - 1.5 * (3.0/26.0) * e3);
                /* -3/11 = -3/22-3/22 (DLMF 19.36.2 eb*dp coefficient). */
                double s2 = eb * (0.5 * (1.0/3.0) + dp * (-(3.0/11.0) + dp * (3.0/26.0)));
                double s3 = dp * ea * ((1.0/3.0) - dp * (3.0/22.0)) - (1.0/3.0) * dp * ec;
                return 3.0 * sigma + fac * (s1 + s2 + s3) / (mu * ml_sqrt(mu));
            }
            {
                double sx = ml_sqrt(xn), sy = ml_sqrt(yn), sz = ml_sqrt(zn);
                double lam = sx * (sy + sz) + sy * sz;
                double alfa = pn * (sx + sy + sz) + sx * sy * sz;
                alfa *= alfa;
                {
                    double beta = pn * (pn + lam) * (pn + lam);
                    double rcv = ml_carlson_rc(alfa, beta);
                    if (!ml_isfinite(rcv)) return ml_make_nan();
                    sigma += fac * rcv;
                }
                fac *= 0.25;
                xn = (xn + lam) * 0.25;
                yn = (yn + lam) * 0.25;
                zn = (zn + lam) * 0.25;
                pn = (pn + lam) * 0.25;
            }
        }
        return ml_make_nan();
    }
}

ML_API double ml_ellip_k(double k) {
    if (ml_isnan(k)) return k;
    {
        double ak = ml_fabs(k);
        if (ak > 1.0) return ml_make_nan();
        if (ak == 1.0) return ml_make_inf(0);
        if (ak == 0.0) return ML_PI * 0.5;
        /* AGM: K = pi/(2*agm(1,sqrt(1-k^2))). Quadratic, ~4 iters. */
        {
            double a = 1.0, b = ml_sqrt(1.0 - ak * ak);
            for (int i = 0; i < 20; i++) {
                double an = (a + b) * 0.5;
                double bn = ml_sqrt(a * b);
                if (!ml_isfinite(an) || !ml_isfinite(bn)) return ml_make_nan();
                if (an == a) break;
                a = an; b = bn;
            }
            return ML_PI / (2.0 * a);
        }
    }
}

ML_API double ml_ellip_e(double k) {
    if (ml_isnan(k)) return k;
    {
        double ak = ml_fabs(k);
        if (ak > 1.0) return ml_make_nan();
        if (ak == 0.0) return ML_PI * 0.5;
        if (ak == 1.0) return 1.0;
        /* Carlson: E = RF(0,1-k^2,1) - k^2*RD(0,1-k^2,1)/3. */
        {
            double m = ak * ak;
            double rf = ml_carlson_rf(0.0, 1.0 - m, 1.0);
            double rd = ml_carlson_rd(0.0, 1.0 - m, 1.0);
            if (!ml_isfinite(rf) || !ml_isfinite(rd)) return ml_make_nan();
            return rf - m * rd / 3.0;
        }
    }
}

ML_API double ml_ellip_f(double phi, double m) {
    if (ml_isnan(phi) || ml_isnan(m)) return ml_make_nan();
    if (!ml_isfinite(phi) || !ml_isfinite(m)) return ml_make_nan();
    if (m < 0.0 || m > 1.0) return ml_make_nan();
    /* Reduce mod pi: F(n*pi+phi0) = 2nK + F(phi0). */
    {
        double K = ml_ellip_k(ml_sqrt(m));
        double n = ml_round(phi / ML_PI);
        double r = phi - n * ML_PI;
        /* Map r to [-pi/2,pi/2] via symmetry F(-x)=-F(x). */
        double s = ml_sin(r);
        double c = ml_cos(r);
        double rf = ml_carlson_rf(c * c, 1.0 - m * s * s, 1.0);
        if (!ml_isfinite(rf)) return ml_make_nan();
        return 2.0 * n * K + s * rf;
    }
}

ML_API double ml_ellip_e_inc(double phi, double m) {
    if (ml_isnan(phi) || ml_isnan(m)) return ml_make_nan();
    if (!ml_isfinite(phi) || !ml_isfinite(m)) return ml_make_nan();
    if (m < 0.0 || m > 1.0) return ml_make_nan();
    {
        double E = ml_ellip_e(ml_sqrt(m));
        double n = ml_round(phi / ML_PI);
        double r = phi - n * ML_PI;
        double s = ml_sin(r);
        double c = ml_cos(r);
        double rf = ml_carlson_rf(c * c, 1.0 - m * s * s, 1.0);
        double rd = ml_carlson_rd(c * c, 1.0 - m * s * s, 1.0);
        if (!ml_isfinite(rf) || !ml_isfinite(rd)) return ml_make_nan();
        {
            double base = s * rf - m * s * s * s * rd / 3.0;
            return 2.0 * n * E + base;
        }
    }
}

ML_API double ml_hyp2f1(double a, double b, double c, double z) {
    /* Gauss 2F1 for |z|<1 via Taylor t_{n+1}=t_n(a+n)(b+n)z/((c+n)(n+1)),
     * Kahan. c not nonpos-int. Outside disk -> NaN (continuation deferred,
     * documented; prevents wrong answers). */
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(c) || ml_isnan(z)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c)) return ml_make_nan();
    if (!ml_isfinite(z)) return ml_make_nan();
    {
        double rc = ml_round(c);
        if (c == rc && c <= 0.0) return ml_make_nan();
    }
    if (z == 0.0) return 1.0;
    if (z == 1.0) {
        /* Gauss sum c-a-b>0: Gamma(c)Gamma(c-a-b)/Gamma(c-a)Gamma(c-b). */
        if (c - a - b > 0.0) {
            double L = ml_lgamma(c) + ml_lgamma(c - a - b) - ml_lgamma(c - a) - ml_lgamma(c - b);
            return ml_exp(L);
        }
        return ml_make_nan();
    }
    if (ml_fabs(z) >= 1.0) return ml_make_nan();
    {
        double t = 1.0, sum = 1.0, comp = 0.0;
        for (int n = 0; n < 5000; n++) {
            double nn = (double)n;
            t *= (a + nn) * (b + nn) / ((c + nn) * (nn + 1.0)) * z;
            if (!ml_isfinite(t)) return ml_make_nan();
            {
                double w = t - comp;
                double tt = sum + w;
                comp = (tt - sum) - w;
                sum = tt;
            }
            if (ml_fabs(t) < 1e-17 * ml_fabs(sum)) break;
            if (n == 4999) return ml_make_nan();
        }
        return sum;
    }
}

ML_API double ml_hyp1f1(double a, double b, double z) {
    /* Kummer 1F1 via series + Kummer transform for z<0:
     * 1F1(a;b;z)=e^z 1F1(b-a;b;-z) to avoid alternating cancellation. */
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(z)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b)) return ml_make_nan();
    if (!ml_isfinite(z)) {
        if (ml_isinf(z)) {
            if (z > 0.0) return ml_make_inf(0);
            return 0.0;
        }
        return ml_make_nan();
    }
    {
        double rb = ml_round(b);
        if (b == rb && b <= 0.0) return ml_make_nan();
    }
    if (z == 0.0) return 1.0;
    if (z < 0.0) {
        double f = ml_hyp1f1(b - a, b, -z);
        if (!ml_isfinite(f)) return ml_make_nan();
        return ml_exp(z) * f;
    }
    {
        double t = 1.0, sum = 1.0, comp = 0.0;
        for (int n = 0; n < 5000; n++) {
            double nn = (double)n;
            t *= (a + nn) / ((b + nn) * (nn + 1.0)) * z;
            if (!ml_isfinite(t)) {
                if (z > 0.0) return ml_make_inf(0);
                return ml_make_nan();
            }
            {
                double w = t - comp;
                double tt = sum + w;
                comp = (tt - sum) - w;
                sum = tt;
            }
            if (ml_fabs(t) < 1e-17 * ml_fabs(sum)) break;
            if (n == 4999) return ml_make_nan();
        }
        return sum;
    }
}

ML_API double ml_zeta(double s) {
    static const double B[7] = {0.0, 1.0/6.0, -1.0/30.0, 1.0/42.0, -1.0/30.0, 5.0/66.0, -691.0/2730.0};
    if (ml_isnan(s)) return s;
    if (ml_isinf(s)) return (s > 0.0) ? 1.0 : ml_make_nan();
    if (s == 1.0) return ml_make_inf(0);
    if (s > 1.0) {
        const int N = 50;
        double sum = 0.0, comp = 0.0;
        for (int n = 1; n < N; n++) {
            double t = ml_pow((double)n, -s);
            double w = t - comp;
            double tt = sum + w;
            comp = (tt - sum) - w;
            sum = tt;
        }
        {
            double nN = (double)N;
            double tail = ml_pow(nN, -s) * 0.5 + ml_pow(nN, 1.0 - s) / (s - 1.0);
            double poch = s;
            double pw = ml_pow(nN, -s - 1.0);
            for (int k = 1; k <= 6; k++) {
                /* term = B2k/(2k)! * (s)_{2k-1} * N^{-s-2k+1} */
                double fact = 1.0;
                for (int j = 1; j <= 2 * k; j++) fact *= (double)j;
                double rising = 1.0;
                for (int j = 0; j < 2 * k - 1; j++) rising *= s + (double)j;
                (void)poch;
                tail += B[k] / fact * rising * pw;
                pw /= nN * nN;
                if (!ml_isfinite(tail)) break;
            }
            return sum + tail;
        }
    }
    if (s > 0.0) {
        /* 0<s<1: Dirichlet eta via Euler transform (Van Wijngaarden):
         * eta = sum_{k>=0} d_k/2^{k+1}, d_k = k-th forward difference of
         * a_n=(n+1)^{-s}. 60 terms -> ~1e-14. zeta = eta/(1-2^{1-s}). */
        const int K = 60;
        double a[61];
        for (int n = 0; n <= K; n++) {
            a[n] = ml_pow((double)(n + 1), -s);
            if (!ml_isfinite(a[n])) return ml_make_nan();
        }
        {
            double eta = 0.0;
            double inv2 = 0.5;
            for (int k = 0; k <= K; k++) {
                eta += a[0] * inv2;
                inv2 *= 0.5;
                for (int j = 0; j < K - k; j++) {
                    a[j] = a[j] - a[j + 1];
                }
            }
            {
                double d = 1.0 - ml_pow(2.0, 1.0 - s);
                if (d == 0.0 || !ml_isfinite(d)) return ml_make_nan();
                return eta / d;
            }
        }
    }
    /* s <= 0: reflection (1-s>1, EM path above, no recursion loop). */
    {
        /* Trivial zeros at negative even integers: sin(pi s/2)=0.
         * The product p2*ppi*sn*g*z1 would be 0*Inf=NaN for large
         * |s| (ppi underflows, gamma overflows); return 0 directly. */
        double rs = ml_round(s);
        if (s == rs && s <= 0.0) {
            long long n = (long long)rs;
            /* n even -> zeta(n)=0; guard cast range. */
            if (ml_fabs(rs) < 9007199254740992.0 && (n % 2LL) == 0LL) return 0.0;
        }
        double t = 1.0 - s;
        double z1 = ml_zeta(t);
        if (!ml_isfinite(z1)) return ml_make_nan();
        {
            /* Long-double product avoids 0*Inf when ppi underflows
             * while gamma overflows (both finite in extended range). */
            long double p2 = __builtin_powl(2.0L, (long double)s);
            long double ppi = __builtin_powl(3.14159265358979323846264338327950288L, (long double)(s - 1.0));
            long double sn = __builtin_sinl((long double)s * 3.14159265358979323846264338327950288L * 0.5L);
            long double g = (long double)ml_gamma_new(1.0 - s);
            long double r = p2 * ppi * sn * g * (long double)z1;
            if (!(r == r)) return ml_make_nan();
            if (r > (long double)1.7976931348623157e308L) return ml_make_inf(0);
            if (r < -(long double)1.7976931348623157e308L) return ml_make_inf(1);
            return (double)r;
        }
    }
}

ML_API double ml_bessel_y0(double x) {
    /* Weber: Y0 = 2/pi[(ln(x/2)+gamma)J0 + series with harmonic numbers].
     * x>0 only; x<=8 series, else Hankel asymptotic. */
    static const double GAM = 0.57721566490153286061;
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return 0.0;
    if (x <= 16.0) {
        double j0 = ml_bessel_j0(x);
        double y = (x * x) * 0.25;
        double t = 1.0, s = 0.0, comp = 0.0;
        double hn = 0.0;
        for (int m = 1; m < 60; m++) {
            hn += 1.0 / (double)m;
            t *= -y / ((double)m * (double)m);
            /* Series needs (-1)^{m+1} H_m y^m/(m!^2) = -t*H_m. */
            {
                double w = (-t * hn - comp);
                double tt = s + w;
                comp = (tt - s) - w;
                s = tt;
            }
            if (ml_fabs(t * hn) < 1e-18 * (ml_fabs(s) + 1.0)) break;
        }
        return (2.0 / ML_PI) * ((ml_log(x * 0.5) + GAM) * j0 + s);
    }
    {
        double p = 1.0 - 9.0 / (128.0 * x * x);
        double q = 1.0 / (8.0 * x);
        double chi = x - 0.78539816339744830962;
        double amp = 0.79788456080286535588 / ml_sqrt(x);
        return amp * (ml_sin(chi) * p + ml_cos(chi) * q);
    }
}

ML_API double ml_bessel_y1(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return 0.0;
    if (x <= 16.0) {
        /* DLMF 10.8.1 (n=1): -2/(pi z) + (2/pi)ln(z/2)J1
         *                    - (z/2)/pi * S,
         * S = sum_{k>=0} (psi(k+1)+psi(k+2)) (-y)^k/(k!(k+1)!). */
        double j1 = ml_bessel_j1(x);
        double y = (x * x) * 0.25;
        double t = 1.0, s = 0.0, comp = 0.0;
        for (int k = 0; k < 120; k++) {
            double ck = ml_digamma((double)k + 1.0) + ml_digamma((double)k + 2.0);
            {
                double w = (t * ck - comp);
                double tt = s + w;
                comp = (tt - s) - w;
                s = tt;
            }
            if (k > 3 && ml_fabs(t * ck) < 1e-17 * ml_fabs(s)) break;
            t *= -y / ((double)(k + 1) * (double)(k + 2));
        }
        {
            double lead = -2.0 / (ML_PI * x);
            double lterm = (2.0 / ML_PI) * ml_log(x * 0.5) * j1;
            double sterm = -(x * 0.5) / ML_PI * s;
            return lead + lterm + sterm;
        }
    }
    {
        double p = 1.0 + 15.0 / (128.0 * x * x);
        double q = 3.0 / (8.0 * x);
        double chi = x - 2.35619449019234492885;
        double amp = 0.79788456080286535588 / ml_sqrt(x);
        return amp * (ml_sin(chi) * p + ml_cos(chi) * q);
    }
}

ML_API double ml_bessel_i0(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);
    {
        double ax = ml_fabs(x);
        if (ax <= 16.0) {
            double y = (ax * ax) * 0.25;
            double t = 1.0, sum = 1.0, comp = 0.0;
            for (int m = 1; m < 120; m++) {
                t *= y / ((double)m * (double)m);
                double w = t - comp;
                double tt = sum + w;
                comp = (tt - sum) - w;
                sum = tt;
                if (ml_fabs(t) < 1e-18 * sum) break;
            }
            return sum;
        }
        {
            double e = ml_exp(ax);
            double amp = e / ml_sqrt(2.0 * ML_PI * ax);
            double corr = 1.0 + 1.0 / (8.0 * ax) + 9.0 / (128.0 * ax * ax);
            return amp * corr;
        }
    }
}

ML_API double ml_bessel_i1(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return x;
    if (ml_isinf(x)) return ml_copysign(ml_make_inf(0), x);
    {
        double ax = ml_fabs(x);
        double sgn = (x < 0.0) ? -1.0 : 1.0;
        if (ax <= 16.0) {
            double y = (ax * ax) * 0.25;
            double t = 0.5 * ax, sum = t, comp = 0.0;
            for (int m = 1; m < 120; m++) {
                t *= y / ((double)m * (double)(m + 1));
                double w = t - comp;
                double tt = sum + w;
                comp = (tt - sum) - w;
                sum = tt;
                if (ml_fabs(t) < 1e-18 * sum) break;
            }
            return sgn * sum;
        }
        {
            double e = ml_exp(ax);
            double amp = e / ml_sqrt(2.0 * ML_PI * ax);
            double corr = 1.0 - 3.0 / (8.0 * ax) - 15.0 / (128.0 * ax * ax);
            return sgn * amp * corr;
        }
    }
}

ML_API double ml_bessel_k0(double x) {
    static const double GAM = 0.57721566490153286061;
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return 0.0;
    if (x <= 16.0) {
        double i0 = ml_bessel_i0(x);
        double y = (x * x) * 0.25;
        double t = 1.0, s = 0.0, comp = 0.0;
        double hn = 0.0;
        for (int m = 1; m < 60; m++) {
            hn += 1.0 / (double)m;
            t *= y / ((double)m * (double)m);
            {
                double w = (t * hn - comp);
                double tt = s + w;
                comp = (tt - s) - w;
                s = tt;
            }
            if (ml_fabs(t * hn) < 1e-18 * (ml_fabs(s) + 1.0)) break;
        }
        return -(ml_log(x * 0.5) + GAM) * i0 + s;
    }
    {
        double e = ml_exp(-x);
        double amp = e * ml_sqrt(ML_PI / (2.0 * x));
        return amp * (1.0 - 1.0 / (8.0 * x));
    }
}

ML_API double ml_bessel_k1(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return 0.0;
    if (x <= 16.0) {
        /* DLMF 10.31.1 (n=1): 1/z + ln(z/2)I1 - (z/2) S,
         * S = sum (psi(k+1)+psi(k+2)) y^k/(k!(k+1)!). */
        double i1 = ml_bessel_i1(x);
        double y = (x * x) * 0.25;
        double t = 1.0, s = 0.0, comp = 0.0;
        for (int k = 0; k < 120; k++) {
            double ck = ml_digamma((double)k + 1.0) + ml_digamma((double)k + 2.0);
            {
                double w = (t * ck - comp);
                double tt = s + w;
                comp = (tt - s) - w;
                s = tt;
            }
            if (k > 3 && ml_fabs(t * ck) < 1e-17 * ml_fabs(s)) break;
            t *= y / ((double)(k + 1) * (double)(k + 2));
        }
        return 1.0 / x + ml_log(x * 0.5) * i1 - (x * 0.5) * s;
    }
    {
        double e = ml_exp(-x);
        double amp = e * ml_sqrt(ML_PI / (2.0 * x));
        return amp * (1.0 + 3.0 / (8.0 * x));
    }
}

ML_API double ml_airy_ai(double x) {
    /* Ai''=x Ai with Ai(0)=0.355028053887817, Ai'(0)=-0.258819403792807.
     * Taylor a[0],a[1] given, (n+2)(n+1)a_{n+2}=a_{n-1}. Kahan to n=60:
     * ~1 ULP for |x|<=5. Outside: exponential/oscillatory asymptotics. */
    static const double A0 = 0.35502805388781723943;
    static const double A1 = -0.25881940379280679841;
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? 0.0 : ml_make_nan();
    if (ml_fabs(x) <= 5.0) {
        /* a[0]=A0, a[1]=A1, a[2]=0, a[n+2]=a[n-1]/((n+2)(n+1)). Kahan. */
        double a[102];
        a[0] = A0; a[1] = A1;
        {
            double s = A0 + A1 * x, cc = 0.0;
            double xpow = x;
            for (int n = 0; n < 100; n++) {
                double anp2;
                if (n == 0) anp2 = 0.0;
                else anp2 = a[n - 1] / ((double)(n + 2) * (double)(n + 1));
                a[n + 2] = anp2;
                xpow *= x;
                if (n + 2 >= 2) {
                    double term = anp2 * xpow;
                    double w = term - cc;
                    double tt = s + w;
                    cc = (tt - s) - w;
                    s = tt;
                    /* Every 3rd coefficient is exactly 0 (a2=0 cascade);
                     * never break on an exact zero — it is not convergence. */
                    if (term != 0.0 && (n + 2) > 3 && ml_fabs(term) < 1e-18 * ml_fabs(s)) break;
                }
            }
            return s;
        }
    }
    if (x > 5.0) {
        /* Ai(x)~e^{-z}/(2 sqrt(pi) x^{1/4}) u(z), u=1-5/48z+385/4608z^2. */
        double z = 2.0 * ml_pow(x, 1.5) / 3.0;
        double e = ml_exp(-z);
        double u = 1.0 - 5.0 / (48.0 * z) + 385.0 / (4608.0 * z * z);
        return e * u / (2.0 * ml_sqrt(ML_PI) * ml_pow(x, 0.25));
    }
    /* x < -5: Ai(-z)~[sin f - cos g]/... f=1-..., g=5/48z-... */
    {
        double z = -x;
        double zt = 2.0 * ml_pow(z, 1.5) / 3.0;
        double ph = zt + ML_PI / 4.0;
        double f = 1.0 - 385.0 / (4608.0 * zt * zt);
        double g = 5.0 / (48.0 * zt) - 85085.0 / (663552.0 * zt * zt * zt);
        return (ml_sin(ph) * f - ml_cos(ph) * g) / (ml_sqrt(ML_PI) * ml_pow(z, 0.25));
    }
}
