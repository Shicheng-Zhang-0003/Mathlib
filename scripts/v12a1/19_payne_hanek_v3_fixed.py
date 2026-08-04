#!/usr/bin/env python3
"""
19_payne_hanek_v3_fixed.py
Run from the folder that CONTAINS the v12A1 working folder.

Fixes the Payne-Hanek large-argument reduction. The V2 implementation
accumulated the full product x*(2/pi) in double precision, losing the
fractional part for |x| > ~1e10.

This V3 implementation:
- Extracts integer contribution (mod 4) per table term
- Kahan-accumulates only the fractional parts
- Preserves the fractional part regardless of |x| magnitude

Targets:
  v12A1/src/internal/payne_hanek.h
  v12A1/tests/test_edge_trig.c

Usage:
  python3 19_payne_hanek_v3_fixed.py
  python3 19_payne_hanek_v3_fixed.py --force
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_PAYNE_HANEK_V3_FIXED"

NEW_PAYNE_HANEK_H = r'''#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H

/* MATHLIB_V12A1_PAYNE_HANEK_V3_FIXED */

/*
 * RANGE REDUCTION FOR TRIGONOMETRIC FUNCTIONS
 *
 * Two paths:
 *
 * 1. |x| <= 1e6: Fast 2-term Cody-Waite with error-free transforms.
 *    Accurate to < 1 ULP.
 *
 * 2. |x| > 1e6: Payne-Hanek reduction with per-term integer/fractional
 *    separation. For each term in the 2/pi table multiplication:
 *      - shift >= 2: integer divisible by 4, skip (no contribution)
 *      - shift == 1: contributes 2*(prod mod 2) to quadrant
 *      - shift == 0: contributes (prod mod 4) to quadrant
 *      - shift <  0: fractional part accumulated via Kahan summation
 *
 *    This preserves the fractional part regardless of how large the
 *    integer part of x*(2/pi) is.
 */

#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

/* Cody-Waite constants for |x| <= 1e6 */
static const double
ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

/*
 * 2/pi stored as 24-bit chunks (Cephes/musl table).
 *   2/pi = sum_{k=0}^{65} two_over_pi[k] * 2^(-24*(k+1))
 */
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

/* pi/2 high and low parts for reconstruction */
static const double
ML_PH_PI2_HI = 0x1.921fb54442d18p+0,
ML_PH_PI2_LO = 0x1.1a62633145c07p-54;

static const double ML_PH_TWO_OVER_PI = 0.63661977236758134308;

/*
 * Process one term: prod * 2^shift.
 *
 * prod is an exact integer < 2^52.
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
    double ax = ml_fabs(x);
    int sign = (x < 0.0);
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
     * Relevant table range.
     * k_start: first chunk where the term could have a fractional part.
     * k_end: last chunk before terms become negligible.
     */
    int k_start = (E - 77) / 24;
    if (k_start < 0) k_start = 0;
    int k_end = (E + 53) / 24;
    if (k_end > 65) k_end = 65;

    /* Accumulate quadrant and fractional part */
    int n = 0;
    double frac_acc = 0.0;
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
    if (frac > 0.5) {
        frac -= 1.0;
        n = (n + 1) & 3;
    } else if (frac < -0.5) {
        frac += 1.0;
        n = (n + 3) & 3;
    }

    /* Reconstruct: reduced_arg = frac * pi/2 */
    double result = frac * ML_PH_PI2_HI + frac * ML_PH_PI2_LO;

    if (sign) {
        result = -result;
        n = (4 - n) & 3;
    }

    *y = result;
    return n;
}

/* ========================================================================
 * UNIFIED ENTRY POINT
 * ====================================================================== */
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

        double p = fn * ML_PH_PIO2_HI;
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
'''

def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")


def locate_v12a1() -> tuple[Path, Path]:
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "internal" / "payne_hanek.h").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_payne_hanek(v12: Path, force: bool) -> None:
    path = v12 / "src" / "internal" / "payne_hanek.h"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = path.read_text(encoding="utf-8")
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched (v3_fixed)")
        return
    write_text(path, NEW_PAYNE_HANEK_H)


def patch_edge_trig(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "test_edge_trig.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = path.read_text(encoding="utf-8")
    if "MATHLIB_V12A1_PAYNE_HANEK_V3_TEST" in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Replace the NaN assertion for sin(1e50) with a finite-result assertion
    old_assert = '    ASSERT_TRUE(&ctx, ml_isnan(ml_sin(1e50)), "sin(1e50) safely NaN");'
    new_assert = """    /* MATHLIB_V12A1_PAYNE_HANEK_V3_TEST */
    /* Payne-Hanek v3: large arguments produce finite results */
    double s1e50 = ml_sin(1e50);
    ASSERT_TRUE(&ctx, !ml_isnan(s1e50) && !ml_isinf(s1e50), "sin(1e50) is finite");
    ASSERT_TRUE(&ctx, s1e50 >= -1.0 && s1e50 <= 1.0, "sin(1e50) in [-1,1]");
    double c1e50 = ml_cos(1e50);
    ASSERT_TRUE(&ctx, !ml_isnan(c1e50) && !ml_isinf(c1e50), "cos(1e50) is finite");
    ASSERT_TRUE(&ctx, c1e50 >= -1.0 && c1e50 <= 1.0, "cos(1e50) in [-1,1]");
    /* Pythagorean identity must hold even for huge arguments */
    ASSERT_NEAR(&ctx, s1e50*s1e50 + c1e50*c1e50, 1.0, 1e-12, "sin^2+cos^2 at 1e50");
    /* Even larger */
    double s1e300 = ml_sin(1e300);
    ASSERT_TRUE(&ctx, !ml_isnan(s1e300) && !ml_isinf(s1e300), "sin(1e300) is finite");
    ASSERT_TRUE(&ctx, s1e300 >= -1.0 && s1e300 <= 1.0, "sin(1e300) in [-1,1]");"""

    if old_assert in text:
        text = text.replace(old_assert, new_assert, 1)
        write_text(path, text)
    else:
        print(f"  [warn] {path}: could not find sin(1e50) NaN assertion. Manual check needed.")


def archive_self(v12: Path, force: bool) -> None:
    try:
        source = Path(__file__).resolve()
        dest = v12 / "scripts" / "v12a1" / source.name
        if source == dest:
            return
        if dest.exists() and not force:
            print(f"  [skip] {dest}: already archived")
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        print(f"  [archive] {dest}")
    except NameError:
        pass


def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()
    print("=========================================================")
    print("  MATHLIB v12A1: PAYNE-HANEK V3 FIXED")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] payne_hanek.h — v3 fixed implementation")
    patch_payne_hanek(v12, force)

    print("\n[2/2] test_edge_trig.c — updated assertions")
    patch_edge_trig(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Payne-Hanek V3 fixed applied.")
    print("")
    print("  What changed:")
    print("    - Replaced double-precision accumulation with per-term")
    print("      integer/fractional separation")
    print("    - Quadrant extracted mod 4 from integer parts")
    print("    - Fractional parts Kahan-accumulated")
    print("    - Works for full double range up to ~1.8e308")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC")
    print("    cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("    # Then re-run oracle:")
    print("    gcc -std=c99 -O3 -fPIE \\")
    print("      -Iinclude/mathlib -Isrc \\")
    print("      -DMATHLIB_HAS_ORACLE_DATA \\")
    print("      -o build/oracle_check \\")
    print("      tests/test_oracle.c \\")
    print("      -Lbuild -lmathc -lm")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: sin/cos large-argument ULP errors drop from")
    print("  10^15+ down to < 5 ULP.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
