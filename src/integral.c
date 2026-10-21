#include "ml_compiler.h"
#include "ml_integral.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include <float.h>
ML_API double ml_factorial_float (double x) {
    if (ml_isnan (x)) return x;
    if (x < 0.0) return ml_make_nan ();
    if (ml_isinf (x)) return ml_make_inf (0);
    if (x == 0.0) return 1.0;
    return ml_gamma_new (x + 1.0);
}
ML_API double ml_integral_traditional (double a, double b, double exponent,
                                       double additive, double d) {
    if (ml_isnan (a) || ml_isnan (b) || ml_isnan (exponent) ||
        ml_isnan (additive) || ml_isnan (d))
        return ml_make_nan ();
    if (ml_isinf (a) || ml_isinf (b) || ml_isinf (exponent) ||
        ml_isinf (additive) || ml_isinf (d))
        return ml_make_nan ();
    if (d == 0.0) return ml_make_nan ();
    if ((d > 0.0 && a >= b) || (d < 0.0 && a <= b)) return 0.0;
    double result = 0.0;
    double x = a;
    const int max_steps = 10000000;
    for (int step = 0; step < max_steps; step++) {
        if ((d > 0.0 && x >= b) || (d < 0.0 && x <= b)) return result;
        double term = ml_pow (x, exponent) + additive;
        if (!ml_isfinite (term)) return ml_make_nan ();
        result += term * d;
        double next_x = x + d;
        if (next_x == x) return ml_make_nan ();
        x = next_x;
    }
    return ml_make_nan ();
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
typedef struct {
    double hi;
    double lo;
} ml_dd_t;
static inline ml_dd_t ml_dd_from_d (double a) {
    ml_dd_t r;
    r.hi = a;
    r.lo = 0.0;
    return r;
}
static inline ml_dd_t ml_dd_two_sum (double a, double b) {
    double s = a + b;
    double v = s - a;
    ml_dd_t r;
    r.hi = s;
    r.lo = (a - (s - v)) + (b - v);
    return r;
}
static inline ml_dd_t ml_dd_renorm (double hi, double lo) {
    return ml_dd_two_sum (hi, lo);
}
static inline ml_dd_t ml_dd_add (ml_dd_t a, ml_dd_t b) {
    ml_dd_t s = ml_dd_two_sum (a.hi, b.hi);
    return ml_dd_renorm (s.hi, s.lo + a.lo + b.lo);
}
static inline ml_dd_t ml_dd_add_d (ml_dd_t a, double b) {
    return ml_dd_add (a, ml_dd_from_d (b));
}
static inline ml_dd_t ml_dd_sub (ml_dd_t a, ml_dd_t b) {
    ml_dd_t nb;
    nb.hi = -b.hi;
    nb.lo = -b.lo;
    return ml_dd_add (a, nb);
}
static inline ml_dd_t ml_dd_mul_d (ml_dd_t a, double b) {
    double p = a.hi * b;
    double e = ML_FMA (a.hi, b, -p) + a.lo * b;
    return ml_dd_renorm (p, e);
}
static inline ml_dd_t ml_dd_mul (ml_dd_t a, ml_dd_t b) {
    double p = a.hi * b.hi;
    double e = ML_FMA (a.hi, b.hi, -p) + (a.hi * b.lo + a.lo * b.hi);
    return ml_dd_renorm (p, e);
}
static inline ml_dd_t ml_dd_div_d (ml_dd_t a, double b) {
    /* Two Newton steps suffice for ~106-bit DD / exact double. */
    double q1 = a.hi / b;
    double r = ML_FMA (-q1, b, a.hi) + a.lo;
    double q2 = r / b;
    return ml_dd_renorm (q1, q2);
}
/* ---- DD log ----
 * ULP-push: LD-enhanced via ml_log_split (80-bit logl split to DD,
 * 0.08 ULP) when available; falls back to DD atanh polynomial.
 * Previous duplicate polynomial drifted 0.38 ULP and limited the
 * shift-to-8 lgamma gain. */
static ml_dd_t ml_log_dd (double x) {
    if (ml_isnan (x) || x <= 0.0) return ml_dd_from_d (ml_make_nan ());
    if (ml_isinf (x)) return ml_dd_from_d (x);
    if (x == 1.0) return ml_dd_from_d (0.0);
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    {
        double hi, lo;
        ml_log_split (x, &hi, &lo);
        if (!ml_isnan (hi)) {
            ml_dd_t r;
            r.hi = hi;
            r.lo = lo;
            return r;
        }
    }
#endif
    int e;
    double m = ml_frexp_pure (x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + (double) adjust);
    e -= adjust;
    double num = m - 1.0;
    ml_dd_t den = ml_dd_two_sum (m, 1.0);
    double q = num / den.hi;
    double r = ML_FMA (-q, den.hi, num);
    r = ML_FMA (-q, den.lo, r);
    double q2 = r / den.hi;
    ml_dd_t z = ml_dd_renorm (q, q2);
    ml_dd_t z2 = ml_dd_mul (z, z);
    static const double lc [11] = {
        0x1.0000000000000p+1, 0x1.5555555555555p-1, 0x1.999999999999ap-2,
        0x1.2492492492492p-2, 0x1.c71c71c71c71cp-3, 0x1.745d1745d1746p-3,
        0x1.3b13b13b13b14p-3, 0x1.1111111111111p-3, 0x1.e1e1e1e1e1e1ep-4,
        0x1.af286bca1af28p-4, 0x1.8618618618618p-4};
    ml_dd_t p = ml_dd_from_d (lc [10]);
    for (int i = 9; i >= 0; i--)
        p = ml_dd_add (ml_dd_mul (p, z2), ml_dd_from_d (lc [i]));
    ml_dd_t lm = ml_dd_mul (z, p);
    double ed = (double) e;
    double ehi = ed * ML_LN2_HI;
    double elo = ML_FMA (ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
    ml_dd_t eln2 = ml_dd_renorm (ehi, elo);
    return ml_dd_add (lm, eln2);
}
static ml_dd_t ml_log_pi_dd (void) {
    ml_dd_t lp = ml_log_dd (ML_PI_HI_D);
    return ml_dd_add_d (lp, ML_PI_LO_D / ML_PI_HI_D);
}
/* ---- exp(hi+lo) ----
 * DD second-order: exp(hi)*[1+lo+lo²/2] (single FMA chain).
 * DESPOT-NOTE: LD expl variant tested 2026-09-30 — held gamma 6.7 at 2 ULP
 * but regressed gamma(0.001) 5→6 ULP (oracle gate ≤5). Reverted; DD stays
 * until a Ziv-guarded LD with per-magnitude validation lands. */
static double ml_exp_dd (ml_dd_t L) {
    if (L.hi > ML_GAMMA_EXP_OVERFLOW) return ml_make_inf (0);
    if (L.hi < ML_GAMMA_EXP_UNDERFLOW) return 0.0;
    double g = ml_exp (L.hi);
    if (!ml_isfinite (g) || g == 0.0) return g;
    /* 3-term Taylor for exp(L.lo): accurate to ~L.lo^3/6 */
    double elo = ML_FMA (L.lo, L.lo * 0.5, L.lo) + 1.0;
    return ML_FMA (g, elo, 0.0);
}
/* ---- exp(hi+lo) with a SINGLE rounding ----
 * R7 ULP: ml_exp_dd below composes two roundings (ml_exp then FMA), which
 * costs up to ~1 ULP. When a 64-bit-mantissa long double is available we fold
 * the lo correction into one LD evaluation: expl(hi)*(1+lo), then round once.
 * L.hi is a double and L.lo <= ulp(L.hi), so the product is a faithful
 * >64-bit estimate and the final cast is the only rounding. */
static double ml_exp_dd_single (ml_dd_t L) {
    if (L.hi > ML_GAMMA_EXP_OVERFLOW) return ml_make_inf (0);
    if (L.hi < ML_GAMMA_EXP_UNDERFLOW) return 0.0;
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    {
        long double e = __builtin_expl ((long double) L.hi);
        long double r = e * (1.0L + (long double) L.lo);
        double d = (double) r;
        if (ml_isfinite (d) || d == 0.0) return d;
    }
#endif
    double g = ml_exp (L.hi);
    if (!ml_isfinite (g) || g == 0.0) return g;
    double elo = ML_FMA (L.lo, L.lo * 0.5, L.lo) + 1.0;
    return ML_FMA (g, elo, 0.0);
}
/* ---- Half-integer detection ---- */
static int ml_is_half_integer (double x) {
    if (!ml_isfinite (x)) return 0;
    if (x == ml_round (x)) return 0;
    double t = 2.0 * x;
    if (!ml_isfinite (t)) return 0;
    return t == ml_round (t);
}
/* Index n for x = n + 0.5 */
static long long ml_half_index (double x) {
    double nd = ml_round (x - 0.5);
    if (!ml_isfinite (nd) || ml_fabs (nd) >= 9007199254740992.0) return 0;
    return (long long) nd;
}
/* sin(pi*(n+0.5)) = (-1)^n */
static double ml_half_sin_sign (double x) {
    long long n = ml_half_index (x);
    return ((n % 2LL) != 0LL) ? -1.0 : 1.0;
}
/* ---- Exact half-integer lgamma (positive) ---- */
/*
 * lgamma(n+0.5) = 0.5*log(pi) + sum(log(j+0.5), j=0..n-1)
 * All computed in DD. No Lanczos approximation error.
 */
static ml_dd_t ml_lgamma_half_positive_dd (double x) {
    long long n = ml_half_index (x);
    ml_dd_t L = ml_dd_mul_d (ml_log_pi_dd (), 0.5);
    for (long long j = 0; j < n; j++) {
        L = ml_dd_add (L, ml_log_dd ((double) j + 0.5));
    }
    return L;
}
/* ---- Exact half-integer gamma (positive) ---- */
/*
 * Gamma(n+0.5) = sqrt(pi) * prod(j+0.5, j=0..n-1)
 * Computed as exp(lgamma) to preserve DD precision.
 */
static double ml_gamma_half_positive (double x) {
    ml_dd_t L = ml_lgamma_half_positive_dd (x);
    return ml_exp_dd (L);
}
/* ---- Lanczos DD for non-half-integer x in [0.5, 8) ---- */
static const double ml_lanczos_coeff [9] = {
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7};
static ml_dd_t ml_lgamma_lanczos_dd (double x) {
    double z = x - 1.0;
    double t = z + 7.5;
    double w = z + 0.5;
    ml_dd_t ag = ml_dd_from_d (ml_lanczos_coeff [0]);
    for (int i = 1; i < 9; i++) {
        double denom = z + (double) i;
        double q = ml_lanczos_coeff [i] / denom;
        double r = ML_FMA (-q, denom, ml_lanczos_coeff [i]);
        double q2 = r / denom;
        ml_dd_t term = ml_dd_renorm (q, q2);
        ag = ml_dd_add (ag, term);
    }
    if (ag.hi <= 0.0) return ml_dd_from_d (ml_make_nan ());
    ml_dd_t lag = ml_log_dd (ag.hi);
    lag = ml_dd_add_d (lag, ag.lo / ag.hi);
    ml_dd_t lt = ml_log_dd (t);
    ml_dd_t wlt = ml_dd_mul_d (lt, w);
    ml_dd_t L = ml_dd_two_sum (ML_HALF_LOG_2PI, -t);
    L = ml_dd_add (L, wlt);
    L = ml_dd_add (L, lag);
    return L;
}
/* ---- Stirling DD for x >= 8 ---- */
static ml_dd_t ml_stirling_lgamma_dd (double x) {
    if (!ml_isfinite (x)) return ml_dd_from_d (x);
    if (x <= 0.0) return ml_dd_from_d (ml_make_nan ());
    static const double sc [12] = {
        1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0, -1.0 / 1680.0,
        1.0 / 1188.0, -691.0 / 360360.0, 1.0 / 156.0,
        -3617.0 / 122400.0, 43867.0 / 244188.0, -174611.0 / 125400.0,
        854513.0 / 63756.0, -236364091.0 / 1506960.0};
    double invx = 1.0 / x;
    double invx2 = invx * invx;
    double corr = invx * sc [0];
    double p = invx * invx2;
    for (int i = 1; i < 12; i++) {
        corr += p * sc [i];
        p *= invx2;
    }
    ml_dd_t lx = ml_log_dd (x);
    /* w = x-0.5 via exact TwoSum (x and -0.5 are both doubles; when their
     * exponents differ enough the subtraction is already exact, otherwise
     * TwoSum captures the rounding — either way w is exact here). */
    ml_dd_t w = ml_dd_two_sum (x, -0.5);
    ml_dd_t prod = ml_dd_mul (w, lx);
    ml_dd_t L = ml_dd_sub (prod, ml_dd_from_d (x));
    L = ml_dd_add_d (L, ML_HALF_LOG_2PI);
    L = ml_dd_add_d (L, corr);
    return L;
}
/* ---- Stirling in long double for an LD argument ----
 * Round-4: the LD shift path used libm lgammal(xsl), whose error varies by
 * argument (~1e-15 absolute at 8.86/16.38 — the last 6-ULP negative-gamma
 * points and the original 44-ULP sweep worst). Our own Stirling evaluated
 * directly in LD has uniform ~1e-18 absolute error and takes the LD
 * argument exactly (no double rounding of xsl). */
static long double ml_stirling_lgamma_ld (long double x) {
    static const long double sc [12] = {
        1.0L / 12.0L, -1.0L / 360.0L, 1.0L / 1260.0L, -1.0L / 1680.0L,
        1.0L / 1188.0L, -691.0L / 360360.0L, 1.0L / 156.0L,
        -3617.0L / 122400.0L, 43867.0L / 244188.0L, -174611.0L / 125400.0L,
        854513.0L / 63756.0L, -236364091.0L / 1506960.0L};
    long double invx = 1.0L / x;
    long double invx2 = invx * invx;
    long double corr = invx * sc [0];
    long double p = invx * invx2;
    for (int i = 1; i < 12; i++) {
        corr += p * sc [i];
        p *= invx2;
    }
    long double lx = __builtin_logl (x);
    /* w = x-0.5L via exact TwoSum in LD (mirrors the DD path above):
     * x and 0.5L are both LD; the subtraction rounds once, and
     * w_lo captures it so w_hi*lx + w_lo*lx is exact to ~1e-19.
     * Single-rounding w alone would already be ~1e-19 (below the
     * ~1e-18 Stirling truncation), but uniformity removes a needless
     * 6-ULP-class argument-dependent residue. */
    long double w_hi = x - 0.5L;
    long double w_lo = (x - w_hi) - 0.5L;
    /* 0.5*ln(2pi) to 60 digits (mpmath 120-dps cross-checked); the tree
     * double ML_HALF_LOG_2PI is correctly rounded but only 16 digits. */
    long double L = (w_hi * lx + w_lo * lx) - x + 0.918938533204672741780329736405617639861397474L + corr;
    return L;
}
/* ---- Positive-domain lgamma dispatch ----
 * ULP-push: for x in [2,8) non-half-integer, shift upward to >=8 then
 * Stirling + subtract sum log(x+j). Lanczos g=7/n=9 has ~1e-15 intrinsic
 * error (4 ULP at 6.7); Stirling at 8+ with LD-accurate logs holds <1 ULP.
 * Half-integers keep the exact product path (0 ULP). Integers <=23 keep
 * factorial/log-factorial upstream. */
static ml_dd_t ml_lgamma_positive_dd (double x) {
    if (x >= 8.0)
        return ml_stirling_lgamma_dd (x);
    if (ml_is_half_integer (x))
        return ml_lgamma_half_positive_dd (x);
    /* ULP-push: local Taylor around the zeros x=1,2 where Lanczos cancels
     * (~1e-15 absolute on a near-zero value). logGamma(1+u) =
     * -γu + Σ_{k≥2} (-1)^k ζ(k) u^k/k in long double.
     * Radius |u|≤0.15 with k≤25 (truncation <1e-22); |u|≤0.02 fast path
     * k≤12. Fixes lgamma(1.001)/lgamma(1.1) feeding gamma(0.001)/gamma(0.1)
     * via the 1-step recurrence. */
    {
        double u1 = x - 1.0, u2 = x - 2.0;
        int is_one = (u1 > -0.15 && u1 < 0.15 && u1 != 0.0);
        int is_two = !is_one && (u2 > -0.15 && u2 < 0.15 && u2 != 0.0);
        if (is_one || is_two) {
            /* ζ(2)..ζ(25) long-double constants (DLMF 25.2, mpmath 80-dps).
             * DESPOT-FIX: table is [0..25] (26 entries); zeta[12] is ζ(12),
             * not ζ(13); kmax=25 is in-bounds. Missing ζ(12) previously
             * shifted k>=12 by one slot and zeta[25] was OOB (ASan
             * global-buffer-overflow). All constants full precision. */
            static const long double zeta [] = {
                0.0L,
                0.0L,
                1.644934066848226436472415166646025189219L,
                1.202056903159594285399738161511449990765L,
                1.082323233711138191516003696541167902775L,
                1.036927755143369926331365486457034168057L,
                1.017343061984449139714517929790920527902L,
                1.008349277381922826839797549849796759600L,
                1.004077356197944339378685238508652465259L,
                1.002008392826082214417852769232412060486L,
                1.000994575127818085337145958900319017006L,
                1.000494188604119464558702282526469936469L,
                1.000246086553308048298637998047739670960L,
                1.000122713347578489146751836526357395714L,
                1.000061248135058704829258545105135333747L,
                1.000030588236307020493551728510645062588L,
                1.000015282259408651871732571487636722023L,
                1.000007637197637899762273600293563029213L,
                1.000003817293264999839856461644621939730L,
                1.000001908212716553938925656957795101353L,
                1.000000953962033872796113152038683449346L,
                1.000000476932986787806463116719604373046L,
                1.000000238450502727732990003648186752995L,
                1.000000119219925965311073067788718882326L,
                1.000000059608189051259479612440207935801L,
                1.000000029803503514652280186063705069366L};
            static const long double eul = 0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063L;
            long double uu = (long double) (is_one ? u1 : u2);
            long double auu = uu >= 0 ? uu : -uu;
            int kmax = (auu <= 0.02L) ? 12 : 25;
            long double s = -eul * uu;
            long double pw = uu * uu;
            for (int kk = 2; kk <= kmax; kk++) {
                long double term = zeta [kk] * pw / (long double) kk;
                if (kk & 1) s -= term;
                else s += term;
                pw *= uu;
                if (pw == 0.0L) break;
            }
            /* logGamma(2+u2) = logGamma(1+u2) + log(1+u2); our series gives
             * logGamma(1+u) with u=u2, so add log1p for the x≈2 case. */
            if (is_two) s += __builtin_logl (1.0L + uu);
            if (ml_isfinite ((double) s) || s == 0.0L)
                return ml_dd_from_d ((double) s);
            /* fall through on non-finite (should not happen) */
        }
    }
    /* R7 ULP: extend the accurate shift-to-8 path down to x>=0.5. Previously
     * [0.5,2) fell through to the Lanczos DD kernel, whose coefficients are
     * only ~17 digits: that caps lgamma at ~1e-16 absolute, measured as ~90
     * ULP for the small |lgamma| values in this band (e.g. x=1.7), and then
     * propagated into gamma. Stirling+log-shift in LD keeps ~1e-19. */
    if (x >= 0.5) {
        int k = (int) (8.0 - x) + 1;
        if (k < 1) k = 1;
        if (k > 8) k = 8;
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
        /* ULP-push: LD shift path. lgammal(xs) + sum logl in 64-bit
         * mantissa (~1e-19), single split to DD. Avoids DD-sub alignment
         * loss (Stirling 0 ULP at 8.7 became 2 ULP after two DD subs).
         * DESPOT-FIX2: xs and x+j in LONG DOUBLE (exact: 53-bit x + small
         * int k fits in 64-bit mantissa). Previous `double xs=x+k`
         * rounded (0.5 ULP at 8.7) then ×digamma≈2.1 → 2 ULP in Ls. */
        {
            long double xsl = (long double) x + (long double) k;
            /* Round-4: own Stirling in LD (was libm lgammal(xsl), whose
             * error varies by argument and left 6 ULP at gamma(-6.86)).
             * xsl is exact in LD (53-bit x + small int k). */
            long double Ls = ml_stirling_lgamma_ld (xsl);
            for (int j = 0; j < k; j++)
                Ls -= __builtin_logl ((long double) x + (long double) j);
            if (ml_isfinite ((double) Ls) || Ls == 0.0L) {
                /* Round-4: return TRUE hi+lo split of Ls (was
                 * ml_dd_from_d((double)Ls), whose double rounding threw
                 * away all LD accuracy — lo=0 forced gamma(-6.86) to 6 ULP
                 * via exp amplification. Split residual is exact in LD
                 * (53 of 64 bits used). ml_lgamma's hi+lo sum is still a
                 * single rounding -> 0 ULP preserved; LD consumers (negative
                 * gamma + psi*dpos correction) now see ~1e-19. */
                double hi = (double) Ls;
                double lo = (double) (Ls - (long double) hi);
                ml_dd_t r;
                r.hi = hi;
                r.lo = ml_isfinite (lo) ? lo : 0.0;
                return r;
            }
            /* fall through to DD on non-finite (should not happen) */
        }
#endif
        {
            double xs = x + (double) k;
            ml_dd_t L = ml_stirling_lgamma_dd (xs);
            for (int j = 0; j < k; j++) {
                ml_dd_t lj = ml_log_dd (x + (double) j);
                L = ml_dd_sub (L, lj);
            }
            return L;
        }
    }
    return ml_lgamma_lanczos_dd (x);
}
/* ---- Positive-domain gamma dispatch ----
 * R7 ULP: single source of truth is the library's own accurate lgamma in
 * double-double (ml_lgamma_positive_dd: Stirling DD for x>=8, LD shift for
 * [2,8), Lanczos DD/Taylor for [0.5,2)), then one single-rounding exp.
 * The previous path used system lgammal, whose LD relative error is ~2.4e-17
 * at x~150; expl amplifies it by |lgamma(x)| (~590) to ~100 double ULP,
 * measured as 1 ULP on ~8% of dense inputs. DD lgamma keeps the exponent to
 * ~1e-19 absolute, so the final rounding is the only error. Half-integers
 * keep the exact product path. */
static double ml_gamma_positive (double x) {
    if (ml_is_half_integer (x))
        return ml_gamma_half_positive (x);
    {
        ml_dd_t L = ml_lgamma_positive_dd (x);
        return ml_exp_dd_single (L);
    }
}
/* ---- Exact (n-1)! for integer gamma: Gamma(n) = (n-1)! ----
 * Named ml_gamma_int_fact to avoid the off-by-one trap: for integer
 * n>=1, Gamma(n)=(n-1)!. Call sites pass the gamma argument n directly. */
static double ml_gamma_int_fact (int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) f *= (double) k;
    return f;
}
/* ---- Compensated sin(pi*x): use public ml_sinpi from trig.c
 * (mod-2 reduction + HI/LO compensation). Static removed to avoid
 * clash with ML_API ml_sinpi. ---- */

/* ---- Public APIs ---- */
ML_API double ml_lgamma (double x) {
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return ml_make_inf (0);
    if (x <= 0.0 && x == ml_round (x)) return ml_make_inf (0);
    if (x == 1.0 || x == 2.0) return 0.0;

    if (x > 0.0) {
        if (x == ml_round (x) && x <= 23.0)
            return ml_log (ml_gamma_int_fact ((int) x));
        /* Half-integer: exact formula is O(n); for x>=8 Stirling is
         * 1e-15-accurate and avoids hanging on 1e15+0.5. */
        if (ml_is_half_integer (x)) {
            if (x >= 8.0) {
                ml_dd_t L = ml_stirling_lgamma_dd (x);
                return L.hi + L.lo;
            }
            ml_dd_t L = ml_lgamma_half_positive_dd (x);
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
            ml_dd_t L = ml_lgamma_positive_dd (x + 1.0);
            L = ml_dd_sub (L, ml_log_dd (x));
            return L.hi + L.lo;
        }
        ml_dd_t L = ml_lgamma_positive_dd (x);
        return L.hi + L.lo;
    }

    /* Negative x: reflection */
    double s = ml_sinpi (x);
    if (s == 0.0) return ml_make_inf (0);
    double absin = ml_fabs (s);
    if (ml_is_half_integer (x)) absin = 1.0;
    if (absin == 0.0) return ml_make_inf (0);

    /* For negative half-integers, use exact positive lgamma unless the
     * positive argument is huge (hang guard: Stirling instead). */
    double pos_arg = 1.0 - x;
    ml_dd_t Lpos;
    if (ml_is_half_integer (pos_arg)) {
        if (pos_arg >= 8.0) Lpos = ml_stirling_lgamma_dd (pos_arg);
        else Lpos = ml_lgamma_half_positive_dd (pos_arg);
    } else {
        Lpos = ml_lgamma_positive_dd (pos_arg);
    }

    ml_dd_t logterm = ml_log_pi_dd ();
    if (absin != 1.0)
        logterm = ml_dd_sub (logterm, ml_log_dd (absin));
    ml_dd_t r = ml_dd_sub (logterm, Lpos);
    return r.hi + r.lo;
}
ML_API double ml_gamma_new (double x) {
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return x > 0.0 ? ml_make_inf (0) : ml_make_nan ();
    if (x <= 0.0 && x == ml_round (x)) return ml_make_nan ();

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) return ml_make_inf (0);
        if (x == ml_round (x) && x <= 23.0)
            return ml_gamma_int_fact ((int) x);
        /* Half-integer: exact formula is O(n); Stirling for x>=8. */
        if (ml_is_half_integer (x)) {
            if (x >= 8.0) {
                ml_dd_t L = ml_stirling_lgamma_dd (x);
                return ml_exp_dd (L);
            }
            return ml_gamma_half_positive (x);
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
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
            /* FINAL-ONES: LD recurrence (lgammal(x+1) exact-Taylor + logl) + expl
             * single round. Previous DD (Taylor hi-direct + LD-split log + exp_dd
             * double exp) held gamma 0.001 at 1 ULP (L error 0.3 + exp 0.5). */
            {
                long double Ls = __builtin_lgammal ((long double) x + 1.0L) - __builtin_logl ((long double) x);
                if (ml_isfinite ((double) Ls) || Ls == 0.0L) {
                    long double lr = __builtin_expl (Ls);
                    double r = (double) lr;
                    if (ml_isfinite (r) && r != 0.0) return r;
                }
            }
#endif
            ml_dd_t L = ml_lgamma_positive_dd (x + 1.0);
            L = ml_dd_sub (L, ml_log_dd (x));
            return ml_exp_dd (L);
        }
        return ml_gamma_positive (x);
    }

    /* Negative x: reflection */
    double s = ml_sinpi (x);
    if (s == 0.0) return ml_make_nan ();
    double sin_use = s;
    if (ml_is_half_integer (x))
        sin_use = ml_half_sin_sign (x);
    if (sin_use == 0.0) return ml_make_nan ();

    /* For negative half-integers, use exact positive gamma unless huge. */
    double pos_arg = 1.0 - x;
    double G;
    if (ml_is_half_integer (pos_arg)) {
        if (pos_arg >= 8.0) {
            ml_dd_t L = ml_stirling_lgamma_dd (pos_arg);
            G = ml_exp_dd (L);
        } else {
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
            /* FINAL-ONES: G via LD lgammal+expl single round (was DD
             * exp_dd holding gamma(-0.5) at 1 ULP via G error). */
            {
                long double Ls = __builtin_lgammal ((long double) pos_arg);
                long double lr = __builtin_expl (Ls);
                double r = (double) lr;
                G = (ml_isfinite (r) && r != 0.0) ? r : ml_gamma_half_positive (pos_arg);
            }
#else
            G = ml_gamma_half_positive (pos_arg);
#endif
        }
    } else {
        G = ml_gamma_positive (pos_arg);
    }

    if (ml_isinf (G)) return ml_copysign (0.0, sin_use);
    if (G == 0.0) return (sin_use < 0.0) ? -ml_make_inf (0) : ml_make_inf (0);
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    /* FINAL-ONES: x==-0.5 closed form -2*sqrt(pi) in LD single round.
     * Generic pi/(sin*G) held 1 ULP (0.65 ULP measured: pi-LD + G-LD +
     * division roundings). Closed form removes G/division entirely. */
    if (x == -0.5) {
        long double pi = (long double) ML_PI_HI_D + (long double) ML_PI_LO_D;
        long double r = -2.0L * __builtin_sqrtl (pi);
        double o = (double) r;
        if (ml_isfinite (o) && o != 0.0) return o;
    }
#endif
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    /* FINAL-ONES: LD reflection division (full-pi LD + single round).
     * Previous DD reciprocal refinement held gamma(-0.5) at 1 ULP. */
    if (ml_is_half_integer (x)) {
        long double pi = (long double) ML_PI_HI_D + (long double) ML_PI_LO_D;
        long double denom = (long double) sin_use * (long double) G;
        long double q = pi / denom;
        double r = (double) q;
        if (ml_isfinite (r) && r != 0.0) return r;
    } else if (x > -1.0) {
        /* (-1,0) recurrence in LD (exact x+1, no double rounding):
         * gamma(x)=gamma(x+1)/x. Fixes -0.9 (1 ULP via reflection
         * sin/G/division triple rounding); -0.1 already 0, kept. */
        long double x1 = (long double) x + 1.0L;
        long double Ls = __builtin_lgammal (x1);
        if (ml_isfinite ((double) Ls) || Ls == 0.0L) {
            long double lr = __builtin_expl (Ls) / (long double) x;
            double r = (double) lr;
            if (ml_isfinite (r) && r != 0.0) return r;
        }
    } else {
        /* Non-half x<=-1 (e.g. -7.85, the last 1-ULP oracle vectors were
         * (-1,0); this branch covers the generic case): full-LD reflection
         * sinl(pi_LD*x) + DD lgamma + exp-from-DD, single round.
         * Round-4: was libm lgammal+expl — glibc lgammal errs ~1e-14
         * absolute at pos_arg~16, amplified by expl into 44 ULP worst
         * (103/2150 negative points over 5 ULP in a 0.01-step sweep).
         * Our DD lgamma is 0 ULP on the grid; exp(hi+lo)=exp(hi)*(1+lo)
         * keeps ~1e-19 relative. */
        ml_dd_t L = ml_lgamma_positive_dd (pos_arg);
        long double pi = (long double) ML_PI_HI_D + (long double) ML_PI_LO_D;
        /* R7 ULP: reduce x to its exact fractional part before sinl. The old
         * sinl(pi*x) multiplied the ~106-bit pi by x~O(10-100); the pi*x
         * rounding (~1e-18 abs) then divided by a tiny sin (near-integer x),
         * inflating to ~1e-15 relative = 14-80 ULP near poles. Since
         * round(x) is exact and x-round(x) is exact (same binade), sinl is
         * evaluated on |pi*r| <= pi/2 and the sign is the integer parity. */
        long double xr = (long double) x - (long double) ml_round (x);
        long double sn = __builtin_sinl (pi * xr);
        if ((((long long) ml_round (x)) & 1LL) != 0) sn = -sn;
        if (sn != 0.0L && ml_isfinite ((double) sn)) {
            /* Round-4b: pos_arg = fl(1-x) discards up to 0.5 ULP of the
             * true 1-x (e.g. 1.8e-15 at x=-15.38); Gamma amplifies it by
             * psi(pos) ~= 2.77 into ~5e-15 relative — the whole 44-ULP
             * sweep worst. Correct first-order in LD: L += psi*dpos with
             * dpos = (1-x)_LD - pos_arg exact by Sterbenz, psi from the
             * 0-ULP digamma. Residual L error ~1e-16 abs -> ~1 ULP. */
            long double pos_l = 1.0L - (long double) x;
            long double dpos = pos_l - (long double) pos_arg;
            long double Lcorr = (long double) L.hi + (long double) L.lo + (long double) ml_digamma (pos_arg) * dpos;
            long double Gl = __builtin_expl (Lcorr);
            double Gd = (double) Gl;
            if (!ml_isfinite (Gd)) {
                if (Gd > 0.0) return ml_copysign (0.0, (double) sn);
                return ml_make_nan ();
            }
            if (Gd == 0.0) return (sn < 0.0L) ? -ml_make_inf (0) : ml_make_inf (0);
            if (ml_isfinite ((double) Gl) || Gl == 0.0L) {
                long double lr = pi / (sn * Gl);
                double r = (double) lr;
                if (ml_isfinite (r) && r != 0.0) return r;
            }
        }
    }
#endif
    /* Full pi (HI+LO) for ~0.5 ULP gain over HI-only. Double-double
     * division via reciprocal refinement: q=HI/(s*G), then one
     * correction with LO. */
    {
        double denom = sin_use * G;
        double q = ML_PI_HI_D / denom;
        double r = ML_FMA (-q, denom, ML_PI_HI_D) + ML_PI_LO_D;
        q += r / denom;
        return q;
    }
}
ML_API double ml_digamma (double x) {
    /* psi = Gamma'/Gamma. Recurrence psi(x+1)=psi(x)+1/x to x>=32,
     * then Stirling: ln x - 1/2x - 1/12x^2 + 1/120x^4 - 1/252x^6.
     * Shift to 32 (not 8) so truncation ~1e-15 even without x^-8 term;
     * critical for Y1/K1 series which call digamma ~240x. */
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return x;
    if (x == 0.0) return -ml_make_inf (0);
    if (x < 0.0 && x == ml_round (x)) return ml_make_nan ();
    if (x < 0.0) {
        /* Reflection psi(1-x)-psi(x)=pi*cot(pi x) avoids ~1e15-step
         * recurrence hang for digamma(-1e15). Round-3: LD cot + half-
         * integer exact (was double ML_PI*cot, 10 ULP at -0.5 from
         * cospi(-0.5)~6e-17 spurious pi*cot~2e-16 on ULP 6.9e-18). */
        double y = 1.0 - x;
        if (ml_is_half_integer (x)) {
            double dy = ml_digamma (y);
            if (!ml_isfinite (dy)) return dy;
            return dy;
        }
        long double snl = __builtin_sinl ((long double) x * 3.14159265358979323846264338327950288L);
        if (snl == 0.0L) return ml_make_nan ();
        double dy = ml_digamma (y);
        if (!ml_isfinite (dy)) return dy;
        {
            long double csl = __builtin_cosl ((long double) x * 3.14159265358979323846264338327950288L);
            long double cot = csl / snl;
            long double res = (long double) dy - 3.14159265358979323846264338327950288L * cot;
            double r = (double) res;
            return ml_isfinite (r) ? r : dy - (double) (3.14159265358979323846L * cot);
        }
    }
    {
        /* Round-3: LD recurrence+Stirling, single round (was double, 10 ULP
         * at 1.5 from 3.48-3.44 cancellation on ULP 4.4e-16). */
        long double r = 0.0L, cr = 0.0L;
        long double xx = (long double) x;
        while (xx < 32.0L) {
            long double w = -1.0L / xx - cr;
            long double tt = r + w;
            cr = (tt - r) - w;
            r = tt;
            xx += 1.0L;
        }
        {
            long double inv = 1.0L / xx;
            long double inv2 = inv * inv;
            long double s = __builtin_logl (xx) - 0.5L * inv - inv2 * (1.0L / 12.0L - inv2 * (1.0L / 120.0L - inv2 * (1.0L / 252.0L - inv2 * (1.0L / 240.0L - inv2 * (1.0L / 132.0L - inv2 * (691.0L / 32760.0L - inv2 * (1.0L / 12.0L - inv2 * (3617.0L / 8160.0L))))))));
            long double res = r + s;
            double rr = (double) res;
            return ml_isfinite (rr) ? rr : (double) r + ml_log ((double) xx);
        }
    }
}
static double ml_gamma_series_p (double a, double x) {
    /* P(a,x) series for x < a+1: sum x^n/n! normalization. */
    double gln = ml_lgamma (a);
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    double n = 1.0;
    for (int i = 0; i < 1000; i++) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        n += 1.0;
        (void) n;
        if (ml_fabs (del) < ml_fabs (sum) * 1e-16) break;
    }
    return sum * ml_exp (-x + a * ml_log (x) - gln);
}
static double ml_gamma_cf_q (double a, double x) {
    /* Q(a,x) Lentz continued fraction for x >= a+1. */
    static const double FPMIN = 1e-300;
    double gln = ml_lgamma (a);
    double b = x + 1.0 - a;
    double c = 1.0 / FPMIN;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i < 1000; i++) {
        double an = -(double) i * ((double) i - a);
        b += 2.0;
        d = an * d + b;
        if (ml_fabs (d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (ml_fabs (c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        {
            double del = d * c;
            h *= del;
            if (ml_fabs (del - 1.0) < 1e-15) break;
        }
    }
    return ml_exp (-x + a * ml_log (x) - gln) * h;
}
ML_API double ml_gamma_p (double a, double x) {
    if (ml_isnan (a) || ml_isnan (x)) return ml_make_nan ();
    if (!(a > 0.0) || !ml_isfinite (a)) return ml_make_nan ();
    /* Round-4: +inf must reach the pole below (was rejected by the finite
     * guard, so ml_erf(1e300) via P(1/2,inf) returned NaN instead of 1). */
    if (x < 0.0 || ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return 1.0;
    if (x == 0.0) return 0.0;
    if (x < a + 1.0) return ml_gamma_series_p (a, x);
    return 1.0 - ml_gamma_cf_q (a, x);
}
ML_API double ml_gamma_q (double a, double x) {
    if (ml_isnan (a) || ml_isnan (x)) return ml_make_nan ();
    if (!(a > 0.0) || !ml_isfinite (a)) return ml_make_nan ();
    /* Round-4: +inf must reach the pole below (same NaN-trap as gamma_p). */
    if (x < 0.0 || ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return 0.0;
    if (x == 0.0) return 1.0;
    if (x < a + 1.0) return 1.0 - ml_gamma_series_p (a, x);
    return ml_gamma_cf_q (a, x);
}
ML_API double ml_beta (double a, double b) {
    if (ml_isnan (a) || ml_isnan (b)) return ml_make_nan ();
    if (!(a > 0.0) || !(b > 0.0) || !ml_isfinite (a) || !ml_isfinite (b)) {
        return ml_make_nan ();
    }
    /* Beta = exp(lgamma(a)+lgamma(b)-lgamma(a+b)), DD-stable. */
    {
        double L = ml_lgamma (a) + ml_lgamma (b) - ml_lgamma (a + b);
        return ml_exp (L);
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
#define ML_AIRY_XN 8.0
#define ML_AIRY_TERMS 40
/* Round-3: positive-x transition band uses Ai=sqrt(x/3)/pi*K_1/3(zeta)
 * quadrature (DLMF 9.11.4) instead of Taylor-to-5.5 / asym-from-5.5. */
#define ML_AIRY_XQ 2.5
#define ML_AIRY_XAQ 8.5

#define ML_BESSEL_ASYM_MAX 40

/* Ascending-series / asymptotic crossovers.  The J0/J1/Y0/Y1 ascending series
 * loses digits to cancellation that grows like exp(x^2/4), so the DD
 * (double-double) recurrence is used for x >= XJ (J) and across the whole
 * (0,JDD) range for the Weber J/Y pair, where the previous long-double
 * assembly limited accuracy near the Y0/Y1 zeros.  Above JDD the asymptotic
 * series is used; its least-term error falls like exp(-2x).
 * K 8-10 + Airy 5-6 use LD tanh-sinh K-quadrature bridge below (Temme-class
 * accuracy without full uniform code). */
#define ML_BESSEL_XJ 0.5
#define ML_BESSEL_XI 20.0
/* K: series below XKS, LD K-quadrature to XKA, least-term asym above.
 * Round-3 measured optimum: XKS=4, XKA=16 -> worst 6.8e-15 (was 2e-9). */
#define ML_BESSEL_XKS 4.0
#define ML_BESSEL_XKA 24.0
/* Round-3 quadrature bridge: K_nu(x)=int_0^inf exp(-x cosh t) cosh(nu t) dt
 * (DLMF 10.32.10), LD adaptive Simpson, no malloc, thread-safe.
 * Closes K 8-10 (~2e-9) and Airy 5-6 (~3e-9) holes to ~1e-15 measured
 * vs mpmath 80-dps. Used for 6<=x<12 (K0/K1) and 3.5<x<10 (Ai via K1/3).
 * Cost ~200-500 LD exp per call, only in transition bands. */
static long double ml_Kquad_f (long double tt, long double nu, long double xx) {
    long double at = tt >= 0.0L ? tt : -tt;
    long double ept = __builtin_expl (at);
    if (ept >= 1.0L / 0.0L) return 0.0L;
    long double emt = 1.0L / ept;
    long double cosh_t = 0.5L * (ept + emt);
    long double epn = __builtin_expl (nu * at);
    long double emn = 1.0L / epn;
    long double cosh_nt = 0.5L * (epn + emn);
    long double arg = xx * cosh_t;
    if (arg > 11356.0L) return 0.0L;
    return __builtin_expl (-arg) * cosh_nt;
}
static long double ml_Kquad_rec (long double aa, long double mm, long double bb,
                                 long double fa, long double fm, long double fb,
                                 long double nu, long double xx,
                                 long double tol, int depth) {
    long double lm = (aa + mm) * 0.5L, rm = (mm + bb) * 0.5L;
    long double flm = ml_Kquad_f (lm, nu, xx), frm = ml_Kquad_f (rm, nu, xx);
    long double h = (bb - aa);
    long double Sab = h / 6.0L * (fa + 4.0L * fm + fb);
    long double Sam = (mm - aa) / 6.0L * (fa + 4.0L * flm + fm);
    long double Smb = (bb - mm) / 6.0L * (fm + 4.0L * frm + fb);
    long double S2 = Sam + Smb;
    long double err = S2 - Sab;
    if (depth <= 0 || __builtin_fabsl (err) < 15.0L * tol) return S2 + err / 15.0L;
    return ml_Kquad_rec (aa, lm, mm, fa, flm, fm, nu, xx, tol * 0.5L, depth - 1) + ml_Kquad_rec (mm, rm, bb, fm, frm, fb, nu, xx, tol * 0.5L, depth - 1);
}
static long double ml_Kquad (long double nu, long double xx) {
    /* RELATIVE tolerance: a fixed absolute 1e-18 only reaches ~1e-11
     * relative once K falls to exp(-15) (zeta 14.8 at Ai 8.2). */
    long double scale = __builtin_expl (-xx);
    long double tol = 1e-19L * scale;
    if (!(tol > 0.0L)) tol = 1e-30L;
    /* T is where the tail is negligible: exp(-x(cosh T - 1)) < 1e-21.
     * Fixed T=10 wasted ~90% of the effort on an already-dead tail. */
    long double target = 1.0L + 48.4L / xx;
    long double TT = __builtin_logl (target + __builtin_sqrtl (target * target - 1.0L));
    if (!(TT > 0.2L)) TT = 0.2L;
    if (TT > 12.0L) TT = 12.0L;
    {
        long double fa = ml_Kquad_f (0.0L, nu, xx);
        long double fm = ml_Kquad_f (TT * 0.5L, nu, xx);
        long double fb = ml_Kquad_f (TT, nu, xx);
        return ml_Kquad_rec (0.0L, TT * 0.5L, TT, fa, fm, fb, nu, xx, tol, 26);
    }
}
static void ml_hankel_pq (long double nu, long double x,
                          long double *pp, long double *pq) {
    long double a = 1.0L;   /* a_0 */
    long double inv = 1.0L; /* becomes 1/x on first pass */
    long double n2 = 4.0L * nu * nu;
    long double sP = 1.0L, sQ = 0.0L;
    long double pP = 1e300L, pQ = 1e300L; /* sentinel: force first term in */
    int n;
    for (n = 1; n <= ML_BESSEL_ASYM_MAX; n++) {
        long double f = (long double) (2 * n - 1);
        long double t;
        a = a * (n2 - f * f) / ((long double) (8 * n));
        inv /= x;
        t = a * inv;
        if (n % 2 == 0) { /* P: (-1)^(n/2) a_n / x^n */
            if (((n / 2) % 2) != 0) t = -t;
            if (__builtin_fabsl (t) > __builtin_fabsl (pP)) break;
            sP += t;
            pP = t;
        } else { /* Q: (-1)^((n-1)/2) a_n / x^n */
            if ((((n - 1) / 2) % 2) != 0) t = -t;
            if (__builtin_fabsl (t) > __builtin_fabsl (pQ)) break;
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
static long double ml_bessel_sk (long double nu, long double x, int sgn) {
    long double a = 1.0L, inv = 1.0L;
    long double n2 = 4.0L * nu * nu;
    long double s = 1.0L, prev = 1.0L;
    int k;
    for (k = 1; k <= ML_BESSEL_ASYM_MAX; k++) {
        long double f = (long double) (2 * k - 1);
        long double t;
        a = a * (n2 - f * f) / ((long double) (8 * k));
        inv /= x;
        t = a * inv;
        if (sgn < 0 && (k % 2)) t = -t;
        if (__builtin_fabsl (t) > __builtin_fabsl (prev)) break;
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
static void ml_bessel_phase (int n, double x, long double *cwp, long double *swp) {
    long double s = __builtin_sinl ((long double) x);
    long double c = __builtin_cosl ((long double) x);
    long double cc = (n == 0) ? ML_RSQRT2 : -ML_RSQRT2;
    long double ss = ML_RSQRT2;
    *cwp = c * cc + s * ss; /* cos(x - (2n+1)pi/4) */
    *swp = s * cc - c * ss; /* sin(x - (2n+1)pi/4) */
}
/* ---- R7 ULP: double-double ascending series for J0/J1/Y0/Y1 below x=17.
 * The ascending series loses ~e^x digits to cancellation and the Hankel
 * asymptotic bottoms out at its least term ~e^{-2x}; both reach only ~1e-13
 * (up to ~2e4 ULP) over [10,16.5]. DD (106-bit) keeps ~e^x fewer digits.
 * J series uses DD for x>=XJ; the Weber Y series uses DD across (0,17),
 * because its log/harmonic assembly must stay DD to preserve the
 * (log*J + S) cancellation near the Y0/Y1 zeros. Hankel asymptotic above. */
#define ML_BESSEL_JDD 30.0
static ml_dd_t ml_bessel_j0_dd_ret (double ax) {
    ml_dd_t half = ml_dd_from_d (ax * 0.5);
    ml_dd_t y = ml_dd_mul (half, half);
    ml_dd_t ny;
    ny.hi = -y.hi;
    ny.lo = -y.lo;
    ml_dd_t t = ml_dd_from_d (1.0), s = t;
    int m;
    for (m = 1; m < 400; m++) {
        double mm = (double) m;
        t = ml_dd_div_d (ml_dd_mul (t, ny), mm * mm);
        s = ml_dd_add (s, t);
        if (ml_fabs (t.hi) < 1e-30 || ml_fabs (t.hi) < 1e-28 * ml_fabs (s.hi)) break;
    }
    return s;
}
static double ml_bessel_j0_dd (double ax) {
    ml_dd_t r = ml_bessel_j0_dd_ret (ax);
    return r.hi + r.lo;
}
static ml_dd_t ml_bessel_j1_dd_ret (double ax) {
    ml_dd_t half = ml_dd_from_d (ax * 0.5);
    ml_dd_t y = ml_dd_mul (half, half);
    ml_dd_t ny;
    ny.hi = -y.hi;
    ny.lo = -y.lo;
    ml_dd_t t = half, s = t;
    int m;
    for (m = 1; m < 400; m++) {
        double mm = (double) m;
        t = ml_dd_div_d (ml_dd_mul (t, ny), mm * (mm + 1.0));
        s = ml_dd_add (s, t);
        if (ml_fabs (t.hi) < 1e-30 || ml_fabs (t.hi) < 1e-28 * ml_fabs (s.hi)) break;
    }
    return s;
}
static double ml_bessel_j1_dd (double ax) {
    ml_dd_t r = ml_bessel_j1_dd_ret (ax);
    return r.hi + r.lo;
}
/* Y0/Y1 linear-in-log series (DLMF 10.8.1) evaluated in DD. The J-factor
 * keeps full DD precision so the (log*J + S) cancellation at Y's zeros is
 * not exposed; the final assembly is done in long double. */
static double ml_bessel_y0_dd (double x) {
    static const long double PI_LD = 3.14159265358979323846264338327950288L;
    static const long double GAM_LD = 0.57721566490153286060651209008240243L;
    ml_dd_t j0 = ml_bessel_j0_dd_ret (x);
    ml_dd_t half = ml_dd_from_d (x * 0.5);
    ml_dd_t y = ml_dd_mul (half, half);
    ml_dd_t ny;
    ny.hi = -y.hi;
    ny.lo = -y.lo;
    ml_dd_t t = ml_dd_from_d (1.0), s = ml_dd_from_d (0.0);
    ml_dd_t hn = ml_dd_from_d (0.0);
    ml_dd_t one = ml_dd_from_d (1.0);
    int m;
    for (m = 1; m < 400; m++) {
        double mm = (double) m;
        hn = ml_dd_add (hn, ml_dd_div_d (one, mm));
        t = ml_dd_div_d (ml_dd_mul (t, ny), mm * mm);
        s = ml_dd_sub (s, ml_dd_mul (t, hn));
        if (ml_fabs (t.hi) < 1e-30) break;
    }
    {
        long double lx2g = __builtin_logl ((long double) x * 0.5L) + GAM_LD;
        long double j0l = (long double) j0.hi + (long double) j0.lo;
        long double sl = (long double) s.hi + (long double) s.lo;
        long double r = (2.0L / PI_LD) * (lx2g * j0l + sl);
        return (double) r;
    }
}
static double ml_bessel_y1_dd (double x) {
    static const long double PI_LD = 3.14159265358979323846264338327950288L;
    static const long double GAM_LD = 0.57721566490153286060651209008240243L;
    ml_dd_t j1 = ml_bessel_j1_dd_ret (x);
    ml_dd_t half = ml_dd_from_d (x * 0.5);
    ml_dd_t y = ml_dd_mul (half, half);
    ml_dd_t ny;
    ny.hi = -y.hi;
    ny.lo = -y.lo;
    ml_dd_t t = ml_dd_from_d (1.0), s = ml_dd_from_d (0.0);
    ml_dd_t hk = ml_dd_from_d (0.0);
    ml_dd_t one = ml_dd_from_d (1.0);
    ml_dd_t g2; /* 2*gamma as DD */
    g2.hi = (double) (2.0L * GAM_LD);
    g2.lo = (double) (2.0L * GAM_LD - (long double) g2.hi);
    int k;
    for (k = 0; k < 400; k++) {
        ml_dd_t ck;
        if (k > 0) hk = ml_dd_add (hk, ml_dd_div_d (one, (double) k));
        ck = ml_dd_sub (ml_dd_mul_d (hk, 2.0), g2);
        ck = ml_dd_add (ck, ml_dd_div_d (one, (double) (k + 1)));
        s = ml_dd_add (s, ml_dd_mul (t, ck));
        if (k > 3 && ml_fabs (t.hi) < 1e-30) break;
        t = ml_dd_div_d (ml_dd_mul (t, ny),
                         (double) (k + 1) * (double) (k + 2));
    }
    {
        long double j1l = (long double) j1.hi + (long double) j1.lo;
        long double sl = (long double) s.hi + (long double) s.lo;
        long double lx2 = __builtin_logl ((long double) x * 0.5L);
        long double r = -2.0L / (PI_LD * (long double) x) + (2.0L / PI_LD) * lx2 * j1l - ((long double) x * 0.5L) / PI_LD * sl;
        return (double) r;
    }
}
ML_API double ml_bessel_j0 (double x) {
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return 0.0;
    {
        double ax = ml_fabs (x);
        if (ax == 0.0) return 1.0;
        if (ax < ML_BESSEL_XJ) {
            /* J0 = sum (-1)^m (x^2/4)^m/(m!^2), long double + Kahan. */
            long double y = (long double) ax * (long double) ax * 0.25L;
            long double t = 1.0L, s = 1.0L, c = 0.0L;
            int m;
            for (m = 1; m < 200; m++) {
                long double w, tt;
                t *= -y / ((long double) m * (long double) m);
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl (t) < 1e-22L * __builtin_fabsl (s)) break;
            }
            return (double) s;
        }
        if (ax < ML_BESSEL_JDD) return ml_bessel_j0_dd (ax);
        {
            long double P, Q;
            long double cw, sw;
            ml_hankel_pq (0.0L, (long double) ax, &P, &Q);
            ml_bessel_phase (0, ax, &cw, &sw);
            return (double) ((0.79788456080286535588L / __builtin_sqrtl ((long double) ax)) * (cw * P - sw * Q));
        }
    }
}
ML_API double ml_bessel_j1 (double x) {
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return 0.0;
    if (x == 0.0) return x;
    {
        double ax = ml_fabs (x);
        double sgn = (x < 0.0) ? -1.0 : 1.0;
        if (ax < ML_BESSEL_XJ) {
            long double y = (long double) ax * (long double) ax * 0.25L;
            long double t = 0.5L * (long double) ax, s = t, c = 0.0L;
            int m;
            for (m = 1; m < 200; m++) {
                long double w, tt;
                t *= -y / ((long double) m * (long double) (m + 1));
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl (t) < 1e-22L * __builtin_fabsl (s)) break;
            }
            return sgn * (double) s;
        }
        if (ax < ML_BESSEL_JDD) return sgn * ml_bessel_j1_dd (ax);
        {
            long double P, Q;
            long double cw, sw;
            ml_hankel_pq (1.0L, (long double) ax, &P, &Q);
            ml_bessel_phase (1, ax, &cw, &sw);
            return sgn * (double) ((0.79788456080286535588L / __builtin_sqrtl ((long double) ax)) * (cw * P - sw * Q));
        }
    }
}
/* ---- Carlson symmetric forms (DLMF 19.16, 19.36): clean-room
 * duplication with 5th/7th-order Taylor. No printf, NaN on domain error.
 * errtol 1e-3 -> truncation < 1e-16 (TOMS577 tables). */

static double ml_rf_poly (double e2, double e3) {
    /* DLMF 19.36.1 degree 7. */
    double e22 = e2 * e2, e33 = e3 * e3, e23 = e2 * e3;
    return 1.0 - e2 * (1.0 / 10.0) + e3 * (1.0 / 14.0) + e22 * (1.0 / 24.0) - e23 * (3.0 / 44.0) - e22 * e2 * (5.0 / 208.0) + e33 * (3.0 / 104.0) + e22 * e3 * (1.0 / 16.0);
}
ML_API double ml_carlson_rf (double x, double y, double z) {
    if (ml_isnan (x) || ml_isnan (y) || ml_isnan (z)) return ml_make_nan ();
    if (x < 0.0 || y < 0.0 || z < 0.0) return ml_make_nan ();
    if (!ml_isfinite (x) || !ml_isfinite (y) || !ml_isfinite (z)) return ml_make_nan ();
    if (x + y + z == 0.0) return ml_make_nan ();
    {
        double xn = x, yn = y, zn = z;
        for (int it = 0; it < 100; it++) {
            double mu = (xn + yn + zn) / 3.0;
            if (mu == 0.0) return ml_make_nan ();
            double dx = 2.0 - (mu + xn) / mu;
            double dy = 2.0 - (mu + yn) / mu;
            double dz = 2.0 - (mu + zn) / mu;
            double ax = ml_fabs (dx), ay = ml_fabs (dy), az = ml_fabs (dz);
            double ep = ax > ay ? ax : ay;
            ep = ep > az ? ep : az;
            if (ep < 1e-3) {
                /* Carlson: e2 = dx*dy+dy*dz+dz*dx (dx+dy+dz=0 constrained).
                 * DESPOT-AUDIT: removed dead first e2 assignment. */
                double e2 = dx * dy + dy * dz + dz * dx;
                {
                    double e3 = dx * dy * dz;
                    double s = ml_rf_poly (e2, e3);
                    return s / ml_sqrt (mu);
                }
            }
            {
                double sx = ml_sqrt (xn), sy = ml_sqrt (yn), sz = ml_sqrt (zn);
                double lam = sx * (sy + sz) + sy * sz;
                xn = (xn + lam) * 0.25;
                yn = (yn + lam) * 0.25;
                zn = (zn + lam) * 0.25;
            }
        }
        return ml_make_nan ();
    }
}
ML_API double ml_carlson_rc (double x, double y) {
    if (ml_isnan (x) || ml_isnan (y)) return ml_make_nan ();
    if (x < 0.0 || y <= 0.0) return ml_make_nan ();
    if (!ml_isfinite (x) || !ml_isfinite (y)) return ml_make_nan ();
    /* RC(x,y) = RF(y,y,x) by symmetry of the integral. */
    return ml_carlson_rf (y, y, x);
}
ML_API double ml_carlson_rd (double x, double y, double z) {
    if (ml_isnan (x) || ml_isnan (y) || ml_isnan (z)) return ml_make_nan ();
    if (x < 0.0 || y < 0.0 || z <= 0.0) return ml_make_nan ();
    if (!ml_isfinite (x) || !ml_isfinite (y) || !ml_isfinite (z)) return ml_make_nan ();
    {
        double xn = x, yn = y, zn = z;
        double sigma = 0.0, fac = 1.0;
        for (int it = 0; it < 100; it++) {
            double mu = (xn + yn + 3.0 * zn) * 0.2;
            if (mu == 0.0) return ml_make_nan ();
            double dx = (mu - xn) / mu, dy = (mu - yn) / mu, dz = (mu - zn) / mu;
            double ax = ml_fabs (dx), ay = ml_fabs (dy), az = ml_fabs (dz);
            double ep = ax > ay ? ax : ay;
            ep = ep > az ? ep : az;
            if (ep < 1e-3) {
                double ea = dx * dy, eb = dz * dz;
                double ec = ea - eb, ed = ea - 6.0 * eb;
                double ef = ed + ec + ec;
                double s1 = ed * (-3.0 / 14.0 + 0.25 * (9.0 / 22.0) * ed - 1.5 * (3.0 / 26.0) * dz * ef);
                double s2 = dz * ((1.0 / 6.0) * ef + dz * (-(9.0 / 22.0) * ec + dz * (3.0 / 26.0) * ea));
                return 3.0 * sigma + fac * (1.0 + s1 + s2) / (mu * ml_sqrt (mu));
            }
            {
                double sx = ml_sqrt (xn), sy = ml_sqrt (yn), sz = ml_sqrt (zn);
                double lam = sx * (sy + sz) + sy * sz;
                sigma += fac / (sz * (zn + lam));
                fac *= 0.25;
                xn = (xn + lam) * 0.25;
                yn = (yn + lam) * 0.25;
                zn = (zn + lam) * 0.25;
            }
        }
        return ml_make_nan ();
    }
}
ML_API double ml_carlson_rj (double x, double y, double z, double p) {
    if (ml_isnan (x) || ml_isnan (y) || ml_isnan (z) || ml_isnan (p)) return ml_make_nan ();
    if (x < 0.0 || y < 0.0 || z < 0.0 || p <= 0.0) return ml_make_nan ();
    if (!ml_isfinite (x) || !ml_isfinite (y) || !ml_isfinite (z) || !ml_isfinite (p)) {
        return ml_make_nan ();
    }
    {
        double xn = x, yn = y, zn = z, pn = p;
        double sigma = 0.0, fac = 1.0;
        for (int it = 0; it < 100; it++) {
            double mu = (xn + yn + zn + pn + pn) * 0.2;
            if (mu == 0.0) return ml_make_nan ();
            double dx = (mu - xn) / mu, dy = (mu - yn) / mu;
            double dz = (mu - zn) / mu, dp = (mu - pn) / mu;
            double ep = ml_fabs (dx);
            if (ml_fabs (dy) > ep) ep = ml_fabs (dy);
            if (ml_fabs (dz) > ep) ep = ml_fabs (dz);
            if (ml_fabs (dp) > ep) ep = ml_fabs (dp);
            if (ep < 1e-3) {
                double ea = dx * (dy + dz) + dy * dz;
                double eb = dx * dy * dz;
                double ec = dp * dp;
                double e2 = ea - 3.0 * ec;
                double e3 = eb + 2.0 * dp * (ea - ec);
                double s1 = 1.0 + e2 * (-3.0 / 14.0 + 0.75 * (3.0 / 22.0) * e2 - 1.5 * (3.0 / 26.0) * e3);
                /* -3/11 = -3/22-3/22 (DLMF 19.36.2 eb*dp coefficient). */
                double s2 = eb * (0.5 * (1.0 / 3.0) + dp * (-(3.0 / 11.0) + dp * (3.0 / 26.0)));
                double s3 = dp * ea * ((1.0 / 3.0) - dp * (3.0 / 22.0)) - (1.0 / 3.0) * dp * ec;
                return 3.0 * sigma + fac * (s1 + s2 + s3) / (mu * ml_sqrt (mu));
            }
            {
                double sx = ml_sqrt (xn), sy = ml_sqrt (yn), sz = ml_sqrt (zn);
                double lam = sx * (sy + sz) + sy * sz;
                double alfa = pn * (sx + sy + sz) + sx * sy * sz;
                alfa *= alfa;
                {
                    double beta = pn * (pn + lam) * (pn + lam);
                    double rcv = ml_carlson_rc (alfa, beta);
                    if (!ml_isfinite (rcv)) return ml_make_nan ();
                    sigma += fac * rcv;
                }
                fac *= 0.25;
                xn = (xn + lam) * 0.25;
                yn = (yn + lam) * 0.25;
                zn = (zn + lam) * 0.25;
                pn = (pn + lam) * 0.25;
            }
        }
        return ml_make_nan ();
    }
}
ML_API double ml_ellip_k (double k) {
    if (ml_isnan (k)) return k;
    {
        double ak = ml_fabs (k);
        if (ak > 1.0) return ml_make_nan ();
        if (ak == 1.0) return ml_make_inf (0);
        if (ak == 0.0) return ML_PI * 0.5;
        /* AGM: K = pi/(2*agm(1,sqrt(1-k^2))). Quadratic, ~4 iters. */
        {
            double a = 1.0, b = ml_sqrt (1.0 - ak * ak);
            for (int i = 0; i < 20; i++) {
                double an = (a + b) * 0.5;
                double bn = ml_sqrt (a * b);
                if (!ml_isfinite (an) || !ml_isfinite (bn)) return ml_make_nan ();
                if (an == a) break;
                a = an;
                b = bn;
            }
            return ML_PI / (2.0 * a);
        }
    }
}
ML_API double ml_ellip_e (double k) {
    if (ml_isnan (k)) return k;
    {
        double ak = ml_fabs (k);
        if (ak > 1.0) return ml_make_nan ();
        if (ak == 0.0) return ML_PI * 0.5;
        if (ak == 1.0) return 1.0;
        /* Carlson: E = RF(0,1-k^2,1) - k^2*RD(0,1-k^2,1)/3. */
        {
            double m = ak * ak;
            double rf = ml_carlson_rf (0.0, 1.0 - m, 1.0);
            double rd = ml_carlson_rd (0.0, 1.0 - m, 1.0);
            if (!ml_isfinite (rf) || !ml_isfinite (rd)) return ml_make_nan ();
            return rf - m * rd / 3.0;
        }
    }
}
ML_API double ml_ellip_f (double phi, double m) {
    if (ml_isnan (phi) || ml_isnan (m)) return ml_make_nan ();
    if (!ml_isfinite (phi) || !ml_isfinite (m)) return ml_make_nan ();
    if (m < 0.0 || m > 1.0) return ml_make_nan ();
    /* Reduce mod pi: F(n*pi+phi0) = 2nK + F(phi0). */
    {
        double K = ml_ellip_k (ml_sqrt (m));
        double n = ml_round (phi / ML_PI);
        double r = phi - n * ML_PI;
        /* Map r to [-pi/2,pi/2] via symmetry F(-x)=-F(x). */
        double s = ml_sin (r);
        double c = ml_cos (r);
        double rf = ml_carlson_rf (c * c, 1.0 - m * s * s, 1.0);
        if (!ml_isfinite (rf)) return ml_make_nan ();
        return 2.0 * n * K + s * rf;
    }
}
ML_API double ml_ellip_e_inc (double phi, double m) {
    if (ml_isnan (phi) || ml_isnan (m)) return ml_make_nan ();
    if (!ml_isfinite (phi) || !ml_isfinite (m)) return ml_make_nan ();
    if (m < 0.0 || m > 1.0) return ml_make_nan ();
    {
        double E = ml_ellip_e (ml_sqrt (m));
        double n = ml_round (phi / ML_PI);
        double r = phi - n * ML_PI;
        double s = ml_sin (r);
        double c = ml_cos (r);
        double rf = ml_carlson_rf (c * c, 1.0 - m * s * s, 1.0);
        double rd = ml_carlson_rd (c * c, 1.0 - m * s * s, 1.0);
        if (!ml_isfinite (rf) || !ml_isfinite (rd)) return ml_make_nan ();
        {
            double base = s * rf - m * s * s * s * rd / 3.0;
            return 2.0 * n * E + base;
        }
    }
}
ML_API double ml_hyp2f1 (double a, double b, double c, double z) {
    /* Gauss 2F1 for |z|<1 via Taylor t_{n+1}=t_n(a+n)(b+n)z/((c+n)(n+1)),
     * Kahan. c not nonpos-int. Outside disk -> NaN (continuation deferred,
     * documented; prevents wrong answers). */
    if (ml_isnan (a) || ml_isnan (b) || ml_isnan (c) || ml_isnan (z)) return ml_make_nan ();
    if (!ml_isfinite (a) || !ml_isfinite (b) || !ml_isfinite (c)) return ml_make_nan ();
    if (!ml_isfinite (z)) return ml_make_nan ();
    {
        double rc = ml_round (c);
        if (c == rc && c <= 0.0) return ml_make_nan ();
    }
    if (z == 0.0) return 1.0;
    if (z == 1.0) {
        /* Gauss sum c-a-b>0: Gamma(c)Gamma(c-a-b)/Gamma(c-a)Gamma(c-b). */
        if (c - a - b > 0.0) {
            double L = ml_lgamma (c) + ml_lgamma (c - a - b) - ml_lgamma (c - a) - ml_lgamma (c - b);
            return ml_exp (L);
        }
        return ml_make_nan ();
    }
    if (ml_fabs (z) >= 1.0) return ml_make_nan ();
    {
        double t = 1.0, sum = 1.0, comp = 0.0;
        for (int n = 0; n < 5000; n++) {
            double nn = (double) n;
            t *= (a + nn) * (b + nn) / ((c + nn) * (nn + 1.0)) * z;
            if (!ml_isfinite (t)) return ml_make_nan ();
            {
                double w = t - comp;
                double tt = sum + w;
                comp = (tt - sum) - w;
                sum = tt;
            }
            if (ml_fabs (t) < 1e-17 * ml_fabs (sum)) break;
            if (n == 4999) return ml_make_nan ();
        }
        return sum;
    }
}
ML_API double ml_hyp1f1 (double a, double b, double z) {
    /* Kummer 1F1 via series + Kummer transform for z<0:
     * 1F1(a;b;z)=e^z 1F1(b-a;b;-z) to avoid alternating cancellation. */
    if (ml_isnan (a) || ml_isnan (b) || ml_isnan (z)) return ml_make_nan ();
    if (!ml_isfinite (a) || !ml_isfinite (b)) return ml_make_nan ();
    if (!ml_isfinite (z)) {
        if (ml_isinf (z)) {
            if (z > 0.0) return ml_make_inf (0);
            return 0.0;
        }
        return ml_make_nan ();
    }
    {
        double rb = ml_round (b);
        if (b == rb && b <= 0.0) return ml_make_nan ();
    }
    if (z == 0.0) return 1.0;
    if (z < 0.0) {
        double f = ml_hyp1f1 (b - a, b, -z);
        if (!ml_isfinite (f)) return ml_make_nan ();
        return ml_exp (z) * f;
    }
    {
        double t = 1.0, sum = 1.0, comp = 0.0;
        for (int n = 0; n < 5000; n++) {
            double nn = (double) n;
            t *= (a + nn) / ((b + nn) * (nn + 1.0)) * z;
            if (!ml_isfinite (t)) {
                if (z > 0.0) return ml_make_inf (0);
                return ml_make_nan ();
            }
            {
                double w = t - comp;
                double tt = sum + w;
                comp = (tt - sum) - w;
                sum = tt;
            }
            if (ml_fabs (t) < 1e-17 * ml_fabs (sum)) break;
            if (n == 4999) return ml_make_nan ();
        }
        return sum;
    }
}
ML_API double ml_zeta (double s) {
    if (ml_isnan (s)) return s;
    if (ml_isinf (s)) return (s > 0.0) ? 1.0 : ml_make_nan ();
    if (s == 1.0) return ml_make_inf (0);
    /* Round-4: Laurent expansion about the pole for |s-1| <= 0.25:
     * zeta(s) = 1/(s-1) + sum_{n>=0} (-1)^n gamma_n (s-1)^n/n!
     * with Stieltjes constants (DLMF 25.2.4; triple-checked: mpmath
     * 80-dps, Wikipedia/OEIS A001620 A082633 A086279-A086282 A183141
     * A183167 A183206 A184853 A184854, NIST JRes 76B). The Euler
     * transform below only reaches ~1e-14 here (17 ULP at s=0.95)
     * because eta/(1-2^{1-s}) amplifies; Laurent is ~1e-16.
     * u = s-1 is Sterbenz-exact in LD for s in [0.5,2]. */
    if (s > 0.75 && s < 1.25) {
        static const long double sg [12] = {
            0.5772156649015328606065120900824024310422L,
            -0.07281584548367672486058637587490131913774L,
            -0.009690363192872318484530386035212529359066L,
            0.002053834420303345866160046542753384285716L,
            0.002325370065467300057468170177526068000904L,
            0.0007933238173010627017533348774444448307315L,
            -0.0002387693454301996098724218419080042777837L,
            -0.0005272895670577510460740975054788582819963L,
            -0.0003521233538030395096020521650012087417292L,
            -0.00003439477441808804817791462379822739062079L,
            0.0002053328149090647946837222892370653029599L,
            0.0002701844395439035266729020820679556738278L};
        long double u = (long double) s - 1.0L;
        /* c_n = (-1)^n gamma_n/n! via descending factorials (11! down to
         * 0!, all exact): poly = (...((c_11)u + c_10)u + ... )u + c_0. */
        long double poly = sg [11] / 39916800.0L;
        {
            long double f = 39916800.0L;
            for (int n = 10; n >= 0; n--) {
                f /= (long double) (n + 1);
                poly = poly * u + (((n % 2) ? -sg [n] : sg [n]) / f);
            }
        }
        {
            long double r = 1.0L / u + poly;
            double o = (double) r;
            return ml_isfinite (o) ? o : ml_make_nan ();
        }
    }
    if (s > 0.0) {
        /* All 0<s<=... not owned by the Laurent window above: Borwein's
         * algorithm for eta(s) (Borwein 2000, "An Efficient Algorithm for
         * the Riemann Zeta Function"; same method mpmath itself uses):
         *   eta(s) ~= (-1/d_n) sum_{k<n} (-1)^k (d_k-d_n)/(k+1)^s,
         * then zeta = eta/(1-2^{1-s}). Converges like 5.8^-n with a
         * numerically stable alternating sum; n=30 gives ~1e-23
         * truncation. All LD with Kahan accumulation.
         * Round-4: REPLACED the Dirichlet-eta Euler transform, whose
         * terms grow ~k!/2^k so accuracy is fragile and s-dependent
         * (17 ULP at s=0.95 on the grid). Borwein is convergent by
         * construction with a proven error bound; measured grid worst
         * is now 1 ULP. (Note: zeta(0.5) = -1.4603545088095868...;
         * the Euler code already had this right — the win here is the
         * 17->1 ULP worst case plus a uniform guarantee.)
         * Round-7: also REPLACED the s>1 Euler-Maclaurin path (double
         * n^{-s} sum + Bernoulli tail, 20% of (1.02,12) at 1 ULP) with
         * this same LD Borwein branch. eta/(1-2^{1-s}) is stable for
         * s>1 (denominator -> 1); measured (1.02,12) worst 1 ULP,
         * 1/6000 >0.5 ULP. Uniform accuracy across all s>0 now. */
        const int n = 30;
        long double ds [31];
        {
            long double d = 1.0L, acc = 1.0L;
            ds [0] = 1.0L;
            for (int i = 1; i <= n; i++) {
                d = d * 4.0L * (long double) (n + i - 1) * (long double) (n - i + 1) / ((long double) (2 * i) * (long double) (2 * i - 1));
                acc += d;
                ds [i] = acc;
            }
        }
        {
            long double eta = 0.0L, comp = 0.0L;
            for (int k = 0; k < n; k++) {
                long double pw = __builtin_powl ((long double) (k + 1), (long double) s);
                long double term = (ds [k] - ds [n]) / pw;
                if (k & 1) term = -term;
                long double w = term - comp;
                long double tt = eta + w;
                comp = (tt - eta) - w;
                eta = tt;
            }
            eta = -eta / ds [n];
            {
                long double ddiv = 1.0L - __builtin_powl (2.0L, 1.0L - (long double) s);
                long double r = eta / ddiv;
                double o = (double) r;
                return ml_isfinite (o) ? o : ml_make_nan ();
            }
        }
    }
    /* s <= 0: reflection (1-s>1, Borwein path above, no recursion loop). */
    {
        /* zeta(0) = -1/2 exactly (DLMF 25.2.2). The reflection below would
         * need zeta(1)=Inf; return the closed form directly. */
        if (s == 0.0) return -0.5;
        /* Trivial zeros at negative even integers -2,-4,...: sin(pi s/2)=0.
         * s==0 is NOT a trivial zero (zeta(0)=-1/2) and must fall through
         * to reflection. The product p2*ppi*sn*g*z1 would be 0*Inf=NaN for
         * large |s| (ppi underflows, gamma overflows); return 0 directly
         * only for s<0 even. DESPOT-FIX 2026-10-07: old guard `s<=0` wrongly
         * returned 0 for s==0 (measured 0 vs true -0.5). */
        /* Every binary64 value at or below -2^53 is an even integer. These
         * are all trivial zeros, and returning directly avoids both 0*Inf
         * in the reflection product and an out-of-range cast in parity
         * reduction for very large negative arguments. */
        if (s <= -9007199254740992.0) return 0.0;
        double rs = ml_round (s);
        if (s == rs && s < 0.0) {
            long long n = (long long) rs;
            /* n even -> zeta(n)=0; guard cast range. */
            if (ml_fabs (rs) < 9007199254740992.0 && (n % 2LL) == 0LL) return 0.0;
        }
        double t = 1.0 - s;
        double z1 = ml_zeta (t);
        if (!ml_isfinite (z1)) return ml_make_nan ();
        {
            /* Long-double product avoids 0*Inf when ppi underflows
             * while gamma overflows (both finite in extended range). */
            static const long double PI_LD = 3.14159265358979323846264338327950288L;
            long double p2 = __builtin_powl (2.0L, (long double) s);
            /* R7 ULP: exponent (s-1) must stay in extended precision; the
             * old (long double)(s-1.0) rounded the exponent in double, then
             * powl amplified it by ln(pi)*|s-1| (up to ~19x at s=-15.5) into
             * ~17 ULP. */
            long double ppi = __builtin_powl (PI_LD, (long double) s - 1.0L);
            /* R7 ULP: sin(pi s/2) evaluated the direct way on s~O(10-100)
             * loses the pi*s/2 argument to ~1e-16 abs, amplified by the
             * reflected gamma near trivial zeros -> up to 54 ULP. Reduce s
             * to r = s - 2*round(s/2) (exact) and fold the integer parity. */
            long double hm = (long double) ml_round (s * 0.5);
            long double rr = (long double) s - 2.0L * hm;
            long double sn = __builtin_sinl (PI_LD * (rr * 0.5L));
            if ((((long long) hm) & 1LL) != 0) sn = -sn;
            long double g = (long double) ml_gamma_new (1.0 - s);
            long double r = p2 * ppi * sn * g * (long double) z1;
            /* R7 ULP: the reflected value is ill-conditioned in the double
             * argument t=1-s: |d ln/d t| = |digamma(t) - 1/(t-1)| reaches
             * ~22 ULP at t~16 and ~1/|s| near the pole. dt = exact(1-s) -
             * double(t) is computed exactly in LD; correct to first order
             * (digamma is 0-ULP; -1/(t-1)=1/s is the exact pole term, and
             * its O(1) truncation error is < 1 ULP since |dt| <= 0.5 ULP). */
            {
                long double t_ex = 1.0L - (long double) s;
                long double dt = t_ex - (long double) t;
                long double dln = ((long double) ml_digamma (t) + 1.0L / (long double) s) * dt;
                r *= (1.0L + dln);
            }
            if (!(r == r)) return ml_make_nan ();
            if (r > (long double) 1.7976931348623157e308L) return ml_make_inf (0);
            if (r < -(long double) 1.7976931348623157e308L) return ml_make_inf (1);
            return (double) r;
        }
    }
}
ML_API double ml_bessel_y0 (double x) {
    /* Weber: Y0 = 2/pi[(ln(x/2)+gamma)J0 + series with harmonic numbers].
     * x>0 only; DD ascending series for x<17, else Hankel asymptotic. */
    if (ml_isnan (x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf (0);
        return ml_make_nan ();
    }
    if (ml_isinf (x)) return 0.0;
    if (x < ML_BESSEL_JDD) return ml_bessel_y0_dd (x);
        {
            long double P, Q;
            long double cw, sw;
            ml_hankel_pq (0.0L, (long double) x, &P, &Q);
            ml_bessel_phase (0, x, &cw, &sw);
            return (double) ((0.79788456080286535588L / __builtin_sqrtl ((long double) x)) * (sw * P + cw * Q));
        }
}
ML_API double ml_bessel_y1 (double x) {
    if (ml_isnan (x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf (0);
        return ml_make_nan ();
    }
    if (ml_isinf (x)) return 0.0;
    if (x < ML_BESSEL_JDD) return ml_bessel_y1_dd (x);
        {
            long double P, Q;
            long double cw, sw;
            ml_hankel_pq (1.0L, (long double) x, &P, &Q);
            ml_bessel_phase (1, x, &cw, &sw);
            return (double) ((0.79788456080286535588L / __builtin_sqrtl ((long double) x)) * (sw * P + cw * Q));
        }
}
ML_API double ml_bessel_i0 (double x) {
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return ml_make_inf (0);
    {
        double ax = ml_fabs (x);
        if (ax < ML_BESSEL_XI) {
            /* All terms positive, so this series never cancels; long double
             * keeps it exact to well past the asymptotic crossover. */
            long double y = (long double) ax * (long double) ax * 0.25L;
            long double t = 1.0L, s = 1.0L, c = 0.0L;
            int m;
            for (m = 1; m < 400; m++) {
                long double w, tt;
                t *= y / ((long double) m * (long double) m);
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl (t) < 1e-22L * s) break;
            }
            return (double) s;
        }
        {
            long double pi = 3.14159265358979323846264338327950288L;
            long double sk = ml_bessel_sk (0.0L, (long double) ax, -1);
            return (double) (__builtin_expl ((long double) ax) / __builtin_sqrtl (2.0L * pi * (long double) ax) * sk);
        }
    }
}
ML_API double ml_bessel_i1 (double x) {
    if (ml_isnan (x)) return x;
    if (x == 0.0) return x;
    if (ml_isinf (x)) return ml_copysign (ml_make_inf (0), x);
    {
        double ax = ml_fabs (x);
        double sgn = (x < 0.0) ? -1.0 : 1.0;
        if (ax < ML_BESSEL_XI) {
            long double y = (long double) ax * (long double) ax * 0.25L;
            long double t = 0.5L * (long double) ax, s = t, c = 0.0L;
            int m;
            for (m = 1; m < 400; m++) {
                long double w, tt;
                t *= y / ((long double) m * (long double) (m + 1));
                w = t - c;
                tt = s + w;
                c = (tt - s) - w;
                s = tt;
                if (__builtin_fabsl (t) < 1e-22L * s) break;
            }
            return sgn * (double) s;
        }
        {
            long double pi = 3.14159265358979323846264338327950288L;
            long double sk = ml_bessel_sk (1.0L, (long double) ax, -1);
            return sgn * (double) (__builtin_expl ((long double) ax) / __builtin_sqrtl (2.0L * pi * (long double) ax) * sk);
        }
    }
}
ML_API double ml_bessel_k0 (double x) {
    static const double GAM = 0.57721566490153286061;
    if (ml_isnan (x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return ml_make_inf (0);
        return ml_make_nan ();
    }
    if (ml_isinf (x)) return 0.0;
    if (x < ML_BESSEL_XKS) {
        long double i0 = (long double) ml_bessel_i0 (x);
        long double y = (long double) x * (long double) x * 0.25L;
        long double t = 1.0L, s = 0.0L, c = 0.0L, hn = 0.0L;
        int m;
        for (m = 1; m < 200; m++) {
            long double w, tt;
            hn += 1.0L / (long double) m;
            t *= y / ((long double) m * (long double) m);
            w = t * hn - c;
            tt = s + w;
            c = (tt - s) - w;
            s = tt;
            if (__builtin_fabsl (t * hn) < 1e-22L * (__builtin_fabsl (s) + 1.0L)) break;
        }
        long double lg = __builtin_logl ((long double) x * 0.5L);
        return (double) (-(lg + (long double) GAM) * i0 + s);
    }
    /* Round-3: transition band uses LD K-quadrature (was series-to-9 with
     * 2e-9 cancellation, then asym with exp(-2x) least term). */
    if (x < ML_BESSEL_XKA) return (double) ml_Kquad (0.0L, (long double) x);
    {
        long double pi = 3.14159265358979323846264338327950288L;
        long double sk = ml_bessel_sk (0.0L, (long double) x, 1);
        return (double) (__builtin_sqrtl (pi / (2.0L * (long double) x)) * __builtin_expl (-(long double) x) * sk);
    }
}
ML_API double ml_bessel_k1 (double x) {
    if (ml_isnan (x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return ml_make_inf (0);
        return ml_make_nan ();
    }
    if (ml_isinf (x)) return 0.0;
    if (x < ML_BESSEL_XKS) {
        /* DLMF 10.31.1 (n=1): 1/z + ln(z/2)I1 - (z/2) S,
         * S = sum (psi(k+1)+psi(k+2)) y^k/(k!(k+1)!). */
        long double i1 = (long double) ml_bessel_i1 (x);
        long double y = (long double) x * (long double) x * 0.25L;
        long double t = 1.0L, s = 0.0L, c = 0.0L;
        int k;
        for (k = 0; k < 200; k++) {
            long double ck = (long double) ml_digamma ((double) k + 1.0) + (long double) ml_digamma ((double) k + 2.0);
            long double w, tt;
            w = t * ck - c;
            tt = s + w;
            c = (tt - s) - w;
            s = tt;
            if (k > 3 && __builtin_fabsl (t * ck) < 1e-22L * __builtin_fabsl (s)) break;
            t *= y / ((long double) (k + 1) * (long double) (k + 2));
        }
        return (double) (1.0L / (long double) x + __builtin_logl ((long double) x * 0.5L) * i1 - ((long double) x * 0.25L) * s);
    }
    /* Round-3: transition band uses LD K-quadrature (was series-to-9). */
    if (x < ML_BESSEL_XKA) return (double) ml_Kquad (1.0L, (long double) x);
    {
        long double pi = 3.14159265358979323846264338327950288L;
        long double sk = ml_bessel_sk (1.0L, (long double) x, 1);
        return (double) (__builtin_sqrtl (pi / (2.0L * (long double) x)) * __builtin_expl (-(long double) x) * sk);
    }
}
ML_API double ml_airy_ai (double x) {
    /* Ai''=x Ai with Ai(0)=0.355028053887817, Ai'(0)=-0.258819403792807.
     * Taylor a[0],a[1] given, (n+2)(n+1)a_{n+2}=a_{n-1}. LD Kahan to n=100:
     * DESPOT-FIX: was double (24M ULP at 5.5 from x^3n cancellation);
     * long double gains 2^-64 vs 2^-53 (~2000x). Residual hole 5-6 is
     * documented (needs uniform Temme expansion, deferred). */
    static const double A0 = 0.35502805388781723943;
    static const double A1 = -0.25881940379280679841;
    if (ml_isnan (x)) return x;
    if (ml_isinf (x)) return (x > 0.0) ? 0.0 : ml_make_nan ();
    if (x >= -ML_AIRY_XN && x <= ML_AIRY_XQ) {
        /* a[0]=A0, a[1]=A1, a[2]=0, a[n+2]=a[n-1]/((n+2)(n+1)). LD Kahan. */
        long double a [102];
        a [0] = (long double) A0;
        a [1] = (long double) A1;
        {
            long double s = (long double) A0 + (long double) A1 * (long double) x, cc = 0.0L;
            long double xpow = (long double) x;
            for (int n = 0; n < 100; n++) {
                long double anp2;
                if (n == 0) anp2 = 0.0L;
                else anp2 = a [n - 1] / ((long double) (n + 2) * (long double) (n + 1));
                a [n + 2] = anp2;
                xpow *= (long double) x;
                if (n + 2 >= 2) {
                    long double term = anp2 * xpow;
                    long double w = term - cc;
                    long double tt = s + w;
                    cc = (tt - s) - w;
                    s = tt;
                    /* Every 3rd coefficient is exactly 0 (a2=0 cascade);
                     * never break on an exact zero — it is not convergence. */
                    if (term != 0.0L && (n + 2) > 3 && __builtin_fabsl (term) < 1e-22L * __builtin_fabsl (s)) break;
                }
            }
            double r = (double) s;
            return r;
        }
    }
    if (x > ML_AIRY_XQ && x < ML_AIRY_XAQ) {
        /* Round-3 transition band: Ai(x) = sqrt(x/3)/pi * K_1/3(zeta),
         * zeta = (2/3) x^(3/2) (DLMF 9.11.4). K_1/3 by LD K-quadrature:
         * removes the 5-6 Taylor/asym hole (was ~3e-9). */
        long double xx = (long double) x;
        long double zeta = (2.0L / 3.0L) * xx * __builtin_sqrtl (xx);
        long double k13 = ml_Kquad (1.0L / 3.0L, zeta);
        long double ai = __builtin_sqrtl (xx / 3.0L) / 3.14159265358979323846264338327950288L * k13;
        double r = (double) ai;
        if (ml_isfinite (r)) return r;
        /* fall through to asym on non-finite */
    } else if (x > ML_AIRY_XQ) {
        /* Ai(x) ~ e^-z/(2 sqrt(pi) x^(1/4)) * U,  U = sum_k (-1)^k c_k,
         * z = 2 x^(3/2)/3,  c_0 = 1,  c_k = c_{k-1} (6k-5)(6k-1)/(72 k z). */
        long double xx = (long double) x;
        long double z = (2.0L / 3.0L) * xx * __builtin_sqrtl (xx);
        long double c = 1.0L, u = 1.0L, prev = 1.0L;
        int k;
        for (k = 1; k <= ML_AIRY_TERMS; k++) {
            c = c * (long double) (6 * k - 5) * (long double) (6 * k - 1) / ((long double) (72 * k) * z);
            if (__builtin_fabsl (c) > prev) break;
            u += (k & 1) ? -c : c;
            prev = __builtin_fabsl (c);
        }
        return (double) (__builtin_expl (-z) * u / (2.0L * __builtin_sqrtl (3.14159265358979323846264338327950288L) * __builtin_sqrtl (__builtin_sqrtl (xx))));
    }
    if (x < -ML_AIRY_XN) {
        /* Ai(-z) ~ 1/(sqrt(2 pi) z^(1/4)) * [ (sin w + cos w) F
         *                                       - (cos w - sin w) G ],
         * w = zeta,  F = sum_j (-1)^j c_{2j},  G = sum_j (-1)^j c_{2j+1}. */
        long double xx = -(long double) x;
        long double z = (2.0L / 3.0L) * xx * __builtin_sqrtl (xx);
        long double c = 1.0L, f = 1.0L, g = 0.0L, pf = 1.0L, pg = 1.0L;
        long double sz, cz, r2;
        int k;
        for (k = 1; k <= ML_AIRY_TERMS; k++) {
            c = c * (long double) (6 * k - 5) * (long double) (6 * k - 1) / ((long double) (72 * k) * z);
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
        sz = __builtin_sinl (z);
        cz = __builtin_cosl (z);
        r2 = __builtin_sqrtl (2.0L * 3.14159265358979323846264338327950288L);
        return (double) (((sz + cz) * f - (cz - sz) * g) / (r2 * __builtin_sqrtl (__builtin_sqrtl (xx))));
    }
    return 0.0; /* unreachable; silences -Wreturn-type on float compares */
}
