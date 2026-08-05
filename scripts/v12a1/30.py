#!/usr/bin/env python3
"""
30_a1_payne_hanek_v6_gamma_v5.py
Run from the folder that CONTAINS the v12A1 working folder.

Fixes:
1. Payne-Hanek: V5 fixed-array accumulation lost the fractional part
   because it never reduced mod 8. Restores V3 per-term integer/fractional
   separation with full-table sweep.
2. Gamma: recurrence depth 12 -> 16 for x < 1. Fixes B_22/B_24 Stirling
   coefficients (wrong denominators).

Targets:
  v12A1/src/internal/payne_hanek.h
  v12A1/src/integral.c

Usage:
  python3 30_a1_payne_hanek_v6_gamma_v5.py
  python3 30_a1_payne_hanek_v6_gamma_v5.py --force
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

MARKER_PH  = "MATHLIB_V12A1_PAYNE_HANEK_V6"
MARKER_GAM = "MATHLIB_V12A1_GAMMA_STIRLING_V5"

# ===================================================================
# Payne-Hanek V6: per-term separation + full-table sweep
# ===================================================================
NEW_PAYNE_HANEK = r"""#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H
/* MATHLIB_V12A1_PAYNE_HANEK_V6 */
/*
* RANGE REDUCTION FOR TRIGONOMETRIC FUNCTIONS
*
* Path 1  |x| <= 1e6 : 2-term Cody-Waite with error-free transforms.
* Path 2  |x| >  1e6 : Payne-Hanek with per-term integer/fractional
*                       separation and full-table sweep.
*
* V5 failed because it accumulated the full product x*(2/pi) in
* double precision without reducing mod 8, so the integer part
* (billions) swallowed the fractional part.
*
* V6 restores the V3 per-term approach:
*   - For each table term, extract integer contribution mod 4
*   - Kahan-accumulate only the fractional parts
*   - Full-table sweep (k=0..65) so no bits are dropped
*/
#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

/* ---- Cody-Waite constants (path 1) ---- */
static const double
    ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
    ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

/* ---- 2/pi table, 24-bit chunks (Cephes / musl) ---- */
static const int32_t ml_two_over_pi[66] = {
    0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62,
    0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
    0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
    0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41,
    0x3991D6, 0x398353, 0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8,
    0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
    0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
    0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08,
    0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
    0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880,
    0x4D7327, 0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B
};

/* pi/2 high / low for reconstruction */
static const double
    ML_PH_PI2_HI = 0x1.921fb54442d18p+0,
    ML_PH_PI2_LO = 0x1.1a62633145c07p-54;
static const double ML_PH_TWO_OVER_PI = 0.63661977236758134308;

/*
* Process one term: prod * 2^shift.
*
* prod is an exact integer < 2^53.
* Extracts integer contribution mod 4 into *n.
* Kahan-accumulates fractional part into *acc / *comp.
*/
static inline void ml_ph_process_term(
    double prod, int shift, int *n, double *acc, double *comp
) {
    if (shift >= 2) {
        /* prod * 2^shift is integer divisible by 4. No contribution. */
        return;
    }
    if (shift == 1) {
        /* (prod*2) mod 4 = 2*(prod mod 2) */
        double fl = (double)((long long)prod & 1LL);
        *n = (*n + (int)(2.0 * fl)) & 3;
        return;
    }
    if (shift == 0) {
        /* prod mod 4 */
        double fl = (double)((long long)prod & 3LL);
        *n = (*n + (int)fl) & 3;
        return;
    }
    /* shift < 0: fractional term */
    double t = ml_ldexp_pure(prod, shift);
    double int_part = 0.0;
    double frac_part;
    if (shift >= -52) {
        /* t may have an integer part */
        int_part = (double)(long long)t;
        frac_part = t - int_part;
        *n = (*n + ((int)(long long)int_part & 3)) & 3;
    } else {
        /* t < 1, purely fractional */
        frac_part = t;
    }
    /* Kahan summation of fractional part */
    double y_k = frac_part - *comp;
    double t_k = *acc + y_k;
    *comp = (t_k - *acc) - y_k;
    *acc = t_k;
}

/*
* Payne-Hanek reduction for large arguments (|x| > 1e6).
*
* Computes:
*   *y = x mod (pi/2), in [-pi/4, pi/4]
*   returns quadrant n in {0, 1, 2, 3}
*/
static inline int ml_rem_pio2_large(double x, double *y) {
    uint64_t bits;
    double   ax   = ml_fabs(x);
    int      sign = (x < 0.0);

    memcpy(&bits, &ax, sizeof(uint64_t));
    int biased_e = (int)((bits >> 52) & 0x7FF);
    /* E: exponent such that ax = m * 2^E (m is 53-bit integer) */
    int E = biased_e - 1075;
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);

    /*
    * Split m into high 28 bits and low 25 bits.
    * m_hi has bits 52..25 (28 bits), m_lo has bits 24..0 (25 bits).
    */
    double m_hi = (double)(m >> 25);
    double m_lo = (double)(m & 0x1FFFFFFULL);

    /*
    * Full-table sweep.
    *
    * The skip logic inside ml_ph_process_term() handles terms that
    * don't contribute (shift >= 2). Processing all 66 entries
    * ensures no fractional bits are dropped for any |x|.
    */
    int k_start = 0;
    int k_end   = 65;

    /* Accumulate quadrant and fractional part */
    int    n         = 0;
    double frac_acc  = 0.0;
    double frac_comp = 0.0;

    for (int k = k_start; k <= k_end; k++) {
        double tk = (double)ml_two_over_pi[k];
        int base_shift = E - 24 * k - 24;

        /* m_hi term: m_hi * tk * 2^(base_shift + 25) */
        double prod_hi = m_hi * tk;  /* exact: 28+24 = 52 bits */
        ml_ph_process_term(prod_hi, base_shift + 25,
                           &n, &frac_acc, &frac_comp);

        /* m_lo term: m_lo * tk * 2^base_shift */
        double prod_lo = m_lo * tk;  /* exact: 25+24 = 49 bits */
        ml_ph_process_term(prod_lo, base_shift,
                           &n, &frac_acc, &frac_comp);
    }

    /* Extract integer part from fractional accumulator */
    int extra = (int)(long long)frac_acc;
    n = (n + (extra & 3)) & 3;
    double frac = frac_acc - (double)extra;

    /* Center fractional part in [-0.5, 0.5] */
    if (frac >  0.5) { frac -= 1.0; n = (n + 1) & 3; }
    if (frac < -0.5) { frac += 1.0; n = (n + 3) & 3; }

    /* Reconstruct: reduced_arg = frac * pi/2 */
    double result = ML_FMA(frac, ML_PH_PI2_HI, frac * ML_PH_PI2_LO);

    if (sign) {
        result = -result;
        n = (4 - n) & 3;
    }
    *y = result;
    return n;
}

/* ==================================================================
* UNIFIED ENTRY POINT
* ================================================================== */
static inline int ml_rem_pio2(double x, double *y) {
    if (ml_isnan(x) || ml_isinf(x)) {
        *y = ml_make_nan();
        return 0;
    }
    double ax = ml_fabs(x);

    /*
    * Small/medium arguments: fast Cody-Waite.
    * Accurate to < 1 ULP for |x| <= 1e6.
    */
    if (ax <= 1.0e6) {
        double fn = ml_round(x * ML_PH_TWO_OVER_PI);
        long long n_ll = (long long)fn;
        int n = (int)(n_ll % 4);
        if (n < 0) n += 4;

        double p     = fn * ML_PH_PIO2_HI;
        double p_err = ML_FMA(fn, ML_PH_PIO2_HI, -p);
        double r1, r1_err;
        r1 = ml_two_sum(x, -p, &r1_err);
        double r2 = r1_err - p_err - (fn * ML_PH_PIO2_LO);
        *y = r1 + r2;
        return n;
    }

    /*
    * Large arguments: Payne-Hanek.
    * Works for the full double range up to ~1.8e308.
    */
    return ml_rem_pio2_large(x, y);
}

#endif /* LIBMATHC_PAYNE_HANEK_H */
"""

# ===================================================================
# Gamma V5: deeper recurrence + corrected Stirling coefficients
# ===================================================================
NEW_GAMMA_SECTION = r"""/* MATHLIB_V12A1_GAMMA_STIRLING_V5 */
/*
* Gamma / log-gamma via double-double Stirling expansion.
*
* Changes vs V4:
*   - Recurrence depth 16 for x < 1 (was 12), 8 for x >= 1.
*   - Corrected Stirling coefficients for B_22 and B_24:
*       B_22/(22*21) = 854513/63756   (was 854513/138600)
*       B_24/(24*23) = -236364091/1506960 (was -236364091/6547290)
*   - 12-term Stirling correction retained.
*   - Compensated sin(pi*x) via ml_sinpi().
*   - Half-integer reflection uses |sin|=1 exactly.
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

    /*
    * Corrected coefficients:
    *   B_22/(22*21) = 854513 / (138*462) = 854513/63756
    *   B_24/(24*23) = -236364091 / (2730*552) = -236364091/1506960
    */
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
         854513.0 / 63756.0,
        -236364091.0 / 1506960.0
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
    *   x < 1   -> 16 steps  (deeper than V4's 12, for tiny x)
    *   x >= 1  ->  8 steps  (sufficient for x in [1,8))
    */
    int depth = (x < 1.0) ? 16 : 8;
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
"""

# ===================================================================
# Helpers
# ===================================================================
def fail(msg: str) -> None:
    print("ERROR: " + msg)
    sys.exit(1)

def normalize(t: str) -> str:
    return t.replace("\r\n", "\n").replace("\r", "\n")

def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")

def locate_v12a1() -> tuple:
    root = Path.cwd()
    cand = root / "v12A1"
    if cand.is_dir():
        return root, cand
    if (root / "src" / "internal" / "payne_hanek.h").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/.")

def patch_payne_hanek(v12: Path, force: bool) -> None:
    path = v12 / "src" / "internal" / "payne_hanek.h"
    if not path.is_file():
        fail(f"Missing: {path}")
    if MARKER_PH in path.read_text(encoding="utf-8") and not force:
        print(f"  [skip] {path}: already at V6")
        return
    write_text(path, NEW_PAYNE_HANEK)

def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER_GAM in text and not force:
        print(f"  [skip] {path}: already at GAMMA_STIRLING_V5")
        return

    markers = [
        "MATHLIB_V12A1_GAMMA_STIRLING_V5",
        "MATHLIB_V12A1_GAMMA_STIRLING_V4",
        "MATHLIB_V12A1_GAMMA_STIRLING_V3",
        "MATHLIB_V12A1_GAMMA_STIRLING_V2",
        "MATHLIB_V12A1_GAMMA_STIRLING",
        "MATHLIB_V12A1_GAMMA_DD2",
        "MATHLIB_V12A1_GAMMA_LANCZOS",
    ]
    start = None
    for mk in markers:
        idx = text.find("/* " + mk)
        if idx == -1:
            idx = text.find("/*" + mk)
        if idx >= 0:
            start = idx
            break
    if start is None:
        fail(f"{path}: could not locate gamma section.")

    new_text = text[:start] + NEW_GAMMA_SECTION
    write_text(path, new_text)

def archive_self(v12: Path, force: bool) -> None:
    try:
        src = Path(__file__).resolve()
        dst = v12 / "scripts" / "v12a1" / src.name
        if src == dst:
            return
        if dst.exists() and not force:
            print(f"  [skip] {dst}: already archived")
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        print(f"  [archive] {dst}")
    except NameError:
        pass

def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()

    print("=========================================================")
    print("  MATHLIB v12A1: PAYNE-HANEK V6 + GAMMA STIRLING V5")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] payne_hanek.h — V6 per-term separation")
    patch_payne_hanek(v12, force)

    print("\n[2/2] integral.c — Stirling V5 (depth 16, fixed coefficients)")
    patch_integral(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Both fixes applied.")
    print("")
    print("  Payne-Hanek V6 changes:")
    print("    - Restored per-term integer/fractional separation (V3)")
    print("    - Full-table sweep k=0..65")
    print("    - Kahan accumulation of fractional parts only")
    print("    - Integer part extracted mod 4 per term")
    print("")
    print("  Gamma Stirling V5 changes:")
    print("    - Recurrence depth 16 for x < 1 (was 12)")
    print("    - B_22 coeff: 854513/63756 (was 854513/138600)")
    print("    - B_24 coeff: -236364091/1506960 (was -236364091/6547290)")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("=========================================================")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
