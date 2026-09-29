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
    if (x == 0.0) return -ml_make_inf(0);
    if (x < 0.0 && x == ml_round(x)) return ml_make_nan();
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
        double r = 0.0, cr = 0.0;
        double xx = x;
        while (xx < 32.0) {
            double w = -1.0 / xx - cr;
            double t = r + w;
            cr = (t - r) - w;
            r = t;
            xx += 1.0;
        }
        {
            double inv = 1.0 / xx;
            double inv2 = inv * inv;
            double s = ml_log(xx) - 0.5 * inv
                - inv2 * (1.0/12.0 - inv2 * (1.0/120.0 - inv2 * (1.0/252.0
                  - inv2 * (1.0/240.0 - inv2 * (1.0/132.0 - inv2 * (691.0/32760.0))))));
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

/* ---- Bessel helpers -------------------------------------------------
 * Hankel asymptotic polynomials, DLMF 10.17.3, for integer order nu.
 *   a_k(nu) = prod_{j=1..k} (4 nu^2 - (2j-1)^2) / (k! 8^k)
 *   P = sum_k (-1)^k a_{2k}   / x^{2k}
 *   Q = sum_k (-1)^k a_{2k+1} / x^{2k+1}
 * Each sum is truncated at its least term, so no term count needs tuning.
 * Accumulated in long double: the optimal-truncation series still cancels
 * for moderate x, and the ascending series below cancels harder. */
/* Airy: Taylor series handles [-ML_AIRY_XN, ML_AIRY_XP]; outside, the
 * DLMF 9.7.5/9.7.6 asymptotics are summed to the least term.  The branches
 * cross near |x| ~ 5-7, where the Taylor series starts losing digits to
 * cancellation and the asymptotic starts losing them to its least term. */
#define ML_AIRY_XP 5.5
#define ML_AIRY_XN 7.0
#define ML_AIRY_TERMS 40

#define ML_BESSEL_ASYM_MAX 40

/* Ascending-series / asymptotic crossovers, chosen where the two branches
 * have comparable error.  The ascending series loses digits to cancellation
 * that grows like exp(x^2/4); the asymptotic series bottoms out at its least
 * term, whose error falls like exp(-2x).  The optimum sits where they meet. */
#define ML_BESSEL_XJ 14.0
#define ML_BESSEL_XY 14.0
#define ML_BESSEL_XK 9.0
#define ML_BESSEL_XI 20.0

static void ml_hankel_pq(long double nu, long double x,
                         long double *pp, long double *pq) {
    long double a = 1.0L;                       /* a_0 */
    long double inv = 1.0L;                     /* becomes 1/x on first pass */
    long double n2 = 4.0L * nu * nu;
    long double sP = 1.0L, sQ = 0.0L;
    long double pP = 1e300L, pQ = 1e300L;   /* sentinel: force first term in */
    int n;
    for (n = 1; n <= ML_BESSEL_ASYM_MAX; n++) {
        long double f = (long double)(2 * n - 1);
        long double t;
        a = a * (n2 - f * f) / ((long double)(8 * n));
        inv /= x;
        t = a * inv;
        if (n % 2 == 0) {            /* P: (-1)^(n/2) a_n / x^n */
            if (((n / 2) % 2) != 0) t = -t;
            if (__builtin_fabsl(t) > __builtin_fabsl(pP)) break;
            sP += t;
            pP = t;
        } else {                     /* Q: (-1)^((n-1)/2) a_n / x^n */
            if ((((n - 1) / 2) % 2) != 0) t = -t;
            if (__builtin_fabsl(t) > __builtin_fabsl(pQ)) break;
            sQ += t;
            pQ = t;
        }
    }
    *pp = sP;
    *pq = sQ;
}

/* Modified-Bessel single series, DLMF 10.40.2/10.40.3.
 *   sgn < 0: I_nu = e^x / sqrt(2 pi x) * sum_k (-1)^k a_k / x^k
 *   sgn > 0: K_nu = sqrt(pi / 2x) e^-x * sum_k a_k / x^k          */
static long double ml_bessel_sk(long double nu, long double x, int sgn) {
    long double a = 1.0L, inv = 1.0L;
    long double n2 = 4.0L * nu * nu;
    long double s = 1.0L, prev = 1.0L;
    int k;
    for (k = 1; k <= ML_BESSEL_ASYM_MAX; k++) {
        long double f = (long double)(2 * k - 1);
        long double t;
        a = a * (n2 - f * f) / ((long double)(8 * k));
        inv /= x;
        t = a * inv;
        if (sgn < 0 && (k % 2)) t = -t;
        if (__builtin_fabsl(t) > __builtin_fabsl(prev)) break;
        s += t;
        prev = t;
    }
    return s;
}

/* Oscillatory part of J_nu/Y_nu, evaluated without ever forming
 * x - (2nu+1)pi/4 in double: that subtraction loses every digit of x for
 * large x.  Instead rotate the unit phasor (cos x, sin x) -- which
 * ml_sin/ml_cos already reduce accurately -- by the exact angle.
 * 1/sqrt(2) is exact to 64 bits, so cos/sin of pi/4 and 3pi/4 are the
 * same constant with a sign. */
#define ML_RSQRT2 0.70710678118654752440

static void ml_bessel_phase(int n, double x, double *cwp, double *swp) {
    double s = ml_sin(x), c = ml_cos(x);
    double cc = (n == 0) ? ML_RSQRT2 : -ML_RSQRT2;
    double ss = ML_RSQRT2;
    *cwp = c * cc + s * ss;   /* cos(x - (2n+1)pi/4) */
    *swp = s * cc - c * ss;   /* sin(x - (2n+1)pi/4) */
}

ML_API double ml_bessel_j0(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return 0.0;
    {
        double ax = ml_fabs(x);
        if (ax == 0.0) return 1.0;
        if (ax < ML_BESSEL_XJ) {
            /* J0 = sum (-1)^m (x^2/4)^m/(m!^2), long double + Kahan. */
            long double y = (long double)ax * (long double)ax * 0.25L;
            long double t = 1.0L, s = 1.0L, c = 0.0L;
            int m;
            for (m = 1; m < 200; m++) {
                long double w, tt;
                t *= -y / ((long double)m * (long double)m);
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl(t) < 1e-22L * __builtin_fabsl(s)) break;
            }
            return (double)s;
        }
        {
            long double P, Q;
            double cw, sw;
            ml_hankel_pq(0.0L, (long double)ax, &P, &Q);
            ml_bessel_phase(0, ax, &cw, &sw);
            return (double)((0.79788456080286535588L / __builtin_sqrtl((long double)ax))
                   * ((long double)cw * P - (long double)sw * Q));
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
        if (ax < ML_BESSEL_XJ) {
            long double y = (long double)ax * (long double)ax * 0.25L;
            long double t = 0.5L * (long double)ax, s = t, c = 0.0L;
            int m;
            for (m = 1; m < 200; m++) {
                long double w, tt;
                t *= -y / ((long double)m * (long double)(m + 1));
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl(t) < 1e-22L * __builtin_fabsl(s)) break;
            }
            return sgn * (double)s;
        }
        {
            long double P, Q;
            double cw, sw;
            ml_hankel_pq(1.0L, (long double)ax, &P, &Q);
            ml_bessel_phase(1, ax, &cw, &sw);
            return sgn * (double)((0.79788456080286535588L / __builtin_sqrtl((long double)ax))
                   * ((long double)cw * P - (long double)sw * Q));
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
    if (x < ML_BESSEL_XY) {
        long double j0 = (long double)ml_bessel_j0(x);
        long double y = (long double)x * (long double)x * 0.25L;
        long double t = 1.0L, s = 0.0L, c = 0.0L, hn = 0.0L;
        int m;
        for (m = 1; m < 200; m++) {
            long double w, tt;
            hn += 1.0L / (long double)m;
            t *= -y / ((long double)m * (long double)m);
            w = -t * hn - c;
            tt = s + w;
            c = (tt - s) - w;
            s = tt;
            if (__builtin_fabsl(t * hn) < 1e-22L * (__builtin_fabsl(s) + 1.0L)) break;
        }
        return (double)(2.0L / 3.14159265358979323846264338327950288L
                        * ((ml_log(x * 0.5) + GAM) * (double)j0 + (double)s));
    }
    {
        long double P, Q;
        double cw, sw;
        ml_hankel_pq(0.0L, (long double)x, &P, &Q);
        ml_bessel_phase(0, x, &cw, &sw);
        return (double)((0.79788456080286535588L / __builtin_sqrtl((long double)x))
               * ((long double)sw * P + (long double)cw * Q));
    }
}

ML_API double ml_bessel_y1(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return 0.0;
    if (x < ML_BESSEL_XY) {
        /* DLMF 10.8.1 (n=1): -2/(pi z) + (2/pi)ln(z/2)J1
         *                    - (z/2)/pi * S,
         * S = sum_{k>=0} (psi(k+1)+psi(k+2)) (-y)^k/(k!(k+1)!). */
        long double j1 = (long double)ml_bessel_j1(x);
        long double y = (long double)x * (long double)x * 0.25L;
        long double t = 1.0L, s = 0.0L, c = 0.0L;
        int k;
        for (k = 0; k < 200; k++) {
            long double ck = (long double)ml_digamma((double)k + 1.0)
                           + (long double)ml_digamma((double)k + 2.0);
            long double w, tt;
            w = t * ck - c;
            tt = s + w;
            c = (tt - s) - w;
            s = tt;
            if (k > 3 && __builtin_fabsl(t * ck) < 1e-22L * __builtin_fabsl(s)) break;
            t *= -y / ((long double)(k + 1) * (long double)(k + 2));
        }
        {
            long double pi = 3.14159265358979323846264338327950288L;
            return (double)(-2.0L / (pi * (long double)x)
                   + (2.0L / pi) * (long double)ml_log(x * 0.5) * j1
                   - ((long double)x * 0.5L) / pi * s);
        }
    }
    {
        long double P, Q;
        double cw, sw;
        ml_hankel_pq(1.0L, (long double)x, &P, &Q);
        ml_bessel_phase(1, x, &cw, &sw);
        return (double)((0.79788456080286535588L / __builtin_sqrtl((long double)x))
               * ((long double)sw * P + (long double)cw * Q));
    }
}

ML_API double ml_bessel_i0(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);
    {
        double ax = ml_fabs(x);
        if (ax < ML_BESSEL_XI) {
            /* All terms positive, so this series never cancels; long double
             * keeps it exact to well past the asymptotic crossover. */
            long double y = (long double)ax * (long double)ax * 0.25L;
            long double t = 1.0L, s = 1.0L, c = 0.0L;
            int m;
            for (m = 1; m < 400; m++) {
                long double w, tt;
                t *= y / ((long double)m * (long double)m);
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl(t) < 1e-22L * s) break;
            }
            return (double)s;
        }
        {
            long double pi = 3.14159265358979323846264338327950288L;
            long double sk = ml_bessel_sk(0.0L, (long double)ax, -1);
            return (double)(__builtin_expl((long double)ax)
                   / __builtin_sqrtl(2.0L * pi * (long double)ax) * sk);
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
        if (ax < ML_BESSEL_XI) {
            long double y = (long double)ax * (long double)ax * 0.25L;
            long double t = 0.5L * (long double)ax, s = t, c = 0.0L;
            int m;
            for (m = 1; m < 400; m++) {
                long double w, tt;
                t *= y / ((long double)m * (long double)(m + 1));
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl(t) < 1e-22L * s) break;
            }
            return sgn * (double)s;
        }
        {
            long double pi = 3.14159265358979323846264338327950288L;
            long double sk = ml_bessel_sk(1.0L, (long double)ax, -1);
            return sgn * (double)(__builtin_expl((long double)ax)
                   / __builtin_sqrtl(2.0L * pi * (long double)ax) * sk);
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
    if (x < ML_BESSEL_XK) {
        long double i0 = (long double)ml_bessel_i0(x);
        long double y = (long double)x * (long double)x * 0.25L;
        long double t = 1.0L, s = 0.0L, c = 0.0L, hn = 0.0L;
        int m;
        for (m = 1; m < 200; m++) {
            long double w, tt;
            hn += 1.0L / (long double)m;
            t *= y / ((long double)m * (long double)m);
            w = t * hn - c;
            tt = s + w;
            c = (tt - s) - w;
            s = tt;
            if (__builtin_fabsl(t * hn) < 1e-22L * (__builtin_fabsl(s) + 1.0L)) break;
        }
        return -(ml_log(x * 0.5) + GAM) * (double)i0 + (double)s;
    }
    {
        long double pi = 3.14159265358979323846264338327950288L;
        long double sk = ml_bessel_sk(0.0L, (long double)x, 1);
        return (double)(__builtin_sqrtl(pi / (2.0L * (long double)x))
               * __builtin_expl(-(long double)x) * sk);
    }
}

ML_API double ml_bessel_k1(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return 0.0;
    if (x < ML_BESSEL_XK) {
        /* DLMF 10.31.1 (n=1): 1/z + ln(z/2)I1 - (z/2) S,
         * S = sum (psi(k+1)+psi(k+2)) y^k/(k!(k+1)!). */
        long double i1 = (long double)ml_bessel_i1(x);
        long double y = (long double)x * (long double)x * 0.25L;
        long double t = 1.0L, s = 0.0L, c = 0.0L;
        int k;
        for (k = 0; k < 200; k++) {
            long double ck = (long double)ml_digamma((double)k + 1.0)
                           + (long double)ml_digamma((double)k + 2.0);
            long double w, tt;
            w = t * ck - c;
            tt = s + w;
            c = (tt - s) - w;
            s = tt;
            if (k > 3 && __builtin_fabsl(t * ck) < 1e-22L * __builtin_fabsl(s)) break;
            t *= y / ((long double)(k + 1) * (long double)(k + 2));
        }
        return (double)(1.0L / (long double)x
               + (long double)ml_log(x * 0.5) * i1
               - ((long double)x * 0.25L) * s);
    }
    {
        long double pi = 3.14159265358979323846264338327950288L;
        long double sk = ml_bessel_sk(1.0L, (long double)x, 1);
        return (double)(__builtin_sqrtl(pi / (2.0L * (long double)x))
               * __builtin_expl(-(long double)x) * sk);
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
    if (x >= -ML_AIRY_XN && x <= ML_AIRY_XP) {
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
    if (x > ML_AIRY_XP) {
        /* Ai(x) ~ e^-z/(2 sqrt(pi) x^(1/4)) * U,  U = sum_k (-1)^k c_k,
         * z = 2 x^(3/2)/3,  c_0 = 1,  c_k = c_{k-1} (6k-5)(6k-1)/(72 k z). */
        long double xx = (long double)x;
        long double z = (2.0L / 3.0L) * xx * __builtin_sqrtl(xx);
        long double c = 1.0L, u = 1.0L, prev = 1.0L;
        int k;
        for (k = 1; k <= ML_AIRY_TERMS; k++) {
            c = c * (long double)(6 * k - 5) * (long double)(6 * k - 1)
                / ((long double)(72 * k) * z);
            if (__builtin_fabsl(c) > prev) break;
            u += (k & 1) ? -c : c;
            prev = __builtin_fabsl(c);
        }
        return (double)(__builtin_expl(-z) * u
               / (2.0L * __builtin_sqrtl(3.14159265358979323846264338327950288L)
                  * __builtin_sqrtl(__builtin_sqrtl(xx))));
    }
    if (x < -ML_AIRY_XN) {
        /* Ai(-z) ~ 1/(sqrt(2 pi) z^(1/4)) * [ (sin w + cos w) F
         *                                       - (cos w - sin w) G ],
         * w = zeta,  F = sum_j (-1)^j c_{2j},  G = sum_j (-1)^j c_{2j+1}. */
        long double xx = -(long double)x;
        long double z = (2.0L / 3.0L) * xx * __builtin_sqrtl(xx);
        long double c = 1.0L, f = 1.0L, g = 0.0L, pf = 1.0L, pg = 1.0L;
        long double sz, cz, r2;
        int k;
        for (k = 1; k <= ML_AIRY_TERMS; k++) {
            c = c * (long double)(6 * k - 5) * (long double)(6 * k - 1)
                / ((long double)(72 * k) * z);
            if (k & 1) {
                if (c > pg) break;
                g += ((k / 2) & 1) ? -c : c;
                pg = c;
            } else {
                if (c > pf) break;
                f += ((k / 2) & 1) ? -c : c;
                pf = c;
            }
        }
        /* sin(w+pi/4) = (sin w + cos w)/sqrt2 without forming w+pi/4,
         * which would round away every digit of w once w is large. */
        sz = __builtin_sinl(z);
        cz = __builtin_cosl(z);
        r2 = __builtin_sqrtl(2.0L * 3.14159265358979323846264338327950288L);
        return (double)(((sz + cz) * f - (cz - sz) * g)
               / (r2 * __builtin_sqrtl(__builtin_sqrtl(xx))));
    }
    return 0.0;   /* unreachable; silences -Wreturn-type on float compares */
}
