#!/usr/bin/env python3
"""
07c_payne_hanek_correct_accum.py

Run from the folder that CONTAINS the v12A1 working folder.

Rewrites ml_rem_pio2_large with correct integer/fractional separation.

The v2 algorithm accumulated the full product q = x*(2/pi) in double
precision. For |x| > 1e10 the integer part has hundreds of digits and
the fractional part (which we actually need) is lost.

Fix: for each term, extract the integer contribution (mod 4) and
accumulate only the fractional part via Kahan summation.

Key insight:
  - shift >= 2:  term is an integer divisible by 4 -> no contribution
  - shift == 1:  term mod 4 = 2*(prod mod 2)
  - shift == 0:  term mod 4 = prod mod 4
  - shift <  0:  term has a fractional part -> Kahan-accumulate it

Targets:
    v12A1/src/internal/payne_hanek.h

Usage:
    python3 07c_payne_hanek_correct_accum.py
    python3 07c_payne_hanek_correct_accum.py --force
"""

from __future__ import annotations
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_PAYNE_HANEK_V3"

NEW_PAYNE_HANEK_H = r"""#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H

/* MATHLIB_V12A1_PAYNE_HANEK_V3 */
/*
* RANGE REDUCTION FOR TRIGONOMETRIC FUNCTIONS
*
* Two paths:
*
* 1. |x| <= 1e6: Fast 2-term Cody-Waite with error-free transforms.
*
* 2. |x| > 1e6: Payne-Hanek reduction with correct integer/fractional
*    separation. For each term in the 2/pi table multiplication:
*      - shift >= 2: integer divisible by 4, skip
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

/* ========================================================================
* PATH 1: Cody-Waite constants (|x| <= 1e6)
* ====================================================================== */

static const double
ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

/* ========================================================================
* PATH 2: Payne-Hanek (|x| > 1e6)
*
* 2/pi stored as 24-bit chunks (standard Cephes/musl table).
*   2/pi = sum_{k=0}^{65} two_over_pi[k] * 2^(-24*(k+1))
* ====================================================================== */

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

static const double
ML_PH_PI2_HI = 0x1.921fb54442d18p+0,
ML_PH_PI2_LO = 0x1.1a62633145c07p-54;

static const double ML_PH_TWO_OVER_PI = 0.63661977236758134308;

/*
* Process one term: prod * 2^shift.
*
* prod is an exact integer in double (< 2^53).
* Extracts integer contribution mod 4 into *n.
* Accumulates fractional part into *acc via Kahan summation.
*/
static inline void ml_ph_process_term(
    double prod, int shift, int *n, double *acc, double *comp
) {
    if (shift >= 2) {
        /* prod * 2^shift is an integer divisible by 4. No contribution. */
        return;
    }
    if (shift == 1) {
        /* prod * 2 is an integer. (prod*2) mod 4 = 2*(prod mod 2). */
        double half = prod * 0.5;
        double prod_mod2 = prod - 2.0 * (half - (half - (double)(long long)half));
        /* Simpler: prod is an exact integer, prod mod 2: */
        double fl = (double)((long long)prod & 1LL);
        *n = (*n + (int)(2.0 * fl)) & 3;
        return;
    }
    if (shift == 0) {
        /* prod is an integer. prod mod 4. */
        double fl = (double)((long long)prod & 3LL);
        *n = (*n + (int)fl) & 3;
        return;
    }
    /* shift < 0: fractional term */
    double t = ml_ldexp_pure(prod, shift);
    double int_part = 0.0;
    double frac_part;
    if (shift >= -52) {
        /* t might have an integer part */
        int_part = (double)(long long)t;
        frac_part = t - int_part;
        *n = (*n + ((int)(long long)int_part & 3)) & 3;
    } else {
        /* t < 1, purely fractional */
        frac_part = t;
    }
    /* Kahan summation */
    double y_k = frac_part - *comp;
    double t_k = *acc + y_k;
    *comp = (t_k - *acc) - y_k;
    *acc = t_k;
}

/*
* ml_rem_pio2_large: Payne-Hanek reduction for |x| > 1e6.
*/
static inline int ml_rem_pio2_large(double x, double *y) {
    uint64_t bits;
    double ax = ml_fabs(x);
    int sign = (x < 0.0);
    memcpy(&bits, &ax, sizeof(uint64_t));

    int biased_e = (int)((bits >> 52) & 0x7FF);
    int E = biased_e - 1075;
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);

    /* Split m into high 28 bits and low 25 bits */
    double m_hi = (double)(m >> 25);
    double m_lo = (double)(m & 0x1FFFFFFULL);

    /* Determine relevant table range */
    int k_start = (E - 77) / 24;
    if (k_start < 0) k_start = 0;
    int k_end = (E + 53) / 24;
    if (k_end > 65) k_end = 65;

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
    double total = frac_acc;
    int extra = (int)(long long)total;
    n = (n + (extra & 3)) & 3;
    double frac = total - (double)extra;

    /* Center in [-0.5, 0.5] */
    if (frac > 0.5) {
        frac -= 1.0;
        n = (n + 1) & 3;
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

    return ml_rem_pio2_large(x, y);
}

#endif /* LIBMATHC_PAYNE_HANEK_H */
"""


def fail(msg):
    print("ERROR: " + msg)
    sys.exit(1)

def write_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")

def locate_v12a1():
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "internal" / "payne_hanek.h").is_file():
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/")

def main():
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()

    print("=========================================================")
    print("  MATHLIB v12A1: PAYNE-HANEK V3 (correct accumulation)")
    print("=========================================================")

    path = v12 / "src" / "internal" / "payne_hanek.h"
    if not path.is_file():
        fail(f"Missing: {path}")
    text = path.read_text(encoding="utf-8")
    if MARKER in text and not force:
        print(f"  [skip] {path}: already at v3")
        return 0
    write_text(path, NEW_PAYNE_HANEK_H)

    try:
        src = Path(__file__).resolve()
        dst = v12 / "scripts" / "v12a1" / src.name
        if src != dst:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists() or force:
                shutil.copy2(src, dst)
                print(f"  [archive] {dst}")
    except NameError:
        pass

    print("---------------------------------------------------------")
    print("  Payne-Hanek v3 applied.")
    print("")
    print("  What changed:")
    print("    - Each term's integer part extracted mod 4 separately")
    print("    - Only fractional parts accumulated (Kahan summation)")
    print("    - No precision loss regardless of |x| magnitude")
    print("")
    print("  Verify:")
    print("    cd v12A1 && cmake -B build && cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("    # Then re-run oracle")
    print("=========================================================")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
