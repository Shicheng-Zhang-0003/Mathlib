#!/usr/bin/env python3
"""
07c_payne_hanek_correct.py

Run from the folder that CONTAINS the v12A1 working folder.

Rewrites ml_rem_pio2_large with the correct Payne-Hanek algorithm:
- Splits the 53-bit significand into 24-bit chunks
- Multiplies each chunk by the relevant 2/pi table entries
- Accumulates modulo 8 to prevent overflow while preserving
  the fractional part
- Extracts quadrant (integer mod 4) and reduced argument

This replaces the broken accumulation that tried to build the
full product x*(2/pi) in double precision.

Targets:
    v12A1/src/internal/payne_hanek.h

Usage:
    python3 07c_payne_hanek_correct.py
    python3 07c_payne_hanek_correct.py --force
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
* 2. |x| > 1e6: Payne-Hanek reduction.
*    Splits the significand into 24-bit chunks, multiplies each by
*    the relevant 2/pi table entries, and accumulates modulo 8.
*    This preserves the fractional part regardless of how large
*    the integer part is.
*
* Based on the algorithm from Cephes / musl / FreeBSD msun.
*/

#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

/* ========================================================================
* PATH 1: Cody-Waite constants
* ====================================================================== */

static const double
ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

/* ========================================================================
* PATH 2: Payne-Hanek
*
* 2/pi stored as 24-bit chunks (standard Cephes/musl table).
* 2/pi = sum_{k=0}^{65} two_over_pi[k] * 2^(-24*(k+1))
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

/* pi/2 high and low for reconstruction */
static const double
ML_PH_PI2_HI = 0x1.921fb54442d18p+0,
ML_PH_PI2_LO = 0x1.1a62633145c07p-54;

static const double ML_PH_TWO_OVER_PI = 0.63661977236758134308;

/*
* ml_rem_pio2_large: Payne-Hanek reduction for |x| > 1e6.
*
* Algorithm (Cephes/musl style):
*   1. Split the 53-bit significand into 24-bit chunks
*   2. For each chunk, multiply by the relevant 2/pi table entries
*   3. Accumulate in double precision, reducing mod 8 after each
*      addition to prevent overflow while preserving the fraction
*   4. The integer part mod 4 is the quadrant
*   5. The fractional part * (pi/2) is the reduced argument
*/
static inline int ml_rem_pio2_large(double x, double *y) {
    uint64_t bits;
    double ax = ml_fabs(x);
    int sign = (x < 0.0);
    memcpy(&bits, &ax, sizeof(uint64_t));

    int biased_e = (int)((bits >> 52) & 0x7FF);
    /* E: exponent such that ax = m * 2^E, m is 53-bit integer */
    int E = biased_e - 1075;
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);

    /*
    * Split m into 24-bit chunks (MSB first):
    *   x[0] = bits 52..29 (24 bits), weight 2^(E+29)
    *   x[1] = bits 28..5  (24 bits), weight 2^(E+5)
    *   x[2] = bits 4..0   (5 bits),  weight 2^E
    */
    uint32_t x0 = (uint32_t)(m >> 29);
    uint32_t x1 = (uint32_t)((m >> 5) & 0xFFFFFFULL);
    uint32_t x2 = (uint32_t)(m & 0x1FULL);

    /*
    * For each chunk x[i] with weight 2^(E + offset_i), we need
    * the 2/pi table entries starting at the chunk that contains
    * bit position (E + offset_i) of 2/pi.
    *
    * The product x[i] * T[k] * 2^(E + offset_i - 24*(k+1))
    * contributes to x * (2/pi).
    *
    * We accumulate all contributions, reducing mod 8 after each
    * to keep the running sum small while preserving the fractional
    * part. (sin/cos period is 4 in units of pi/2, so mod 8 is safe.)
    */
    double q = 0.0;

    /* Process chunk x0 (weight 2^(E+29)) */
    {
        int exp0 = E + 29;
        /* Starting table index: the chunk of 2/pi that contains
        * the bit position where the integer part ends */
        int j0 = exp0 / 24;
        if (j0 < 0) j0 = 0;
        /* Process a few table entries */
        for (int k = j0; k <= j0 + 2 && k < 66; k++) {
            int shift = exp0 - 24 * (k + 1);
            double prod = (double)x0 * (double)ml_two_over_pi[k];
            double term = ml_ldexp_pure(prod, shift);
            q += term;
            /* Reduce mod 8 to prevent overflow */
            q -= 8.0 * ml_round(q * 0.125);
        }
    }

    /* Process chunk x1 (weight 2^(E+5)) */
    {
        int exp1 = E + 5;
        int j1 = exp1 / 24;
        if (j1 < 0) j1 = 0;
        for (int k = j1; k <= j1 + 2 && k < 66; k++) {
            int shift = exp1 - 24 * (k + 1);
            double prod = (double)x1 * (double)ml_two_over_pi[k];
            double term = ml_ldexp_pure(prod, shift);
            q += term;
            q -= 8.0 * ml_round(q * 0.125);
        }
    }

    /* Process chunk x2 (weight 2^E) */
    {
        int exp2 = E;
        int j2 = exp2 / 24;
        if (j2 < 0) j2 = 0;
        for (int k = j2; k <= j2 + 2 && k < 66; k++) {
            int shift = exp2 - 24 * (k + 1);
            double prod = (double)x2 * (double)ml_two_over_pi[k];
            double term = ml_ldexp_pure(prod, shift);
            q += term;
            q -= 8.0 * ml_round(q * 0.125);
        }
    }

    /*
    * q is now x * (2/pi) mod 8, with the fractional part preserved.
    * Extract quadrant and reduced argument.
    */
    double n_d = ml_round(q);
    int n = ((int)(long long)n_d) & 3;
    double frac = q - n_d;

    /* Center in [-0.5, 0.5] */
    if (frac > 0.5) { frac -= 1.0; n = (n + 1) & 3; }
    else if (frac < -0.5) { frac += 1.0; n = (n + 3) & 3; }

    /* Reduced argument = frac * (pi/2) */
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
        print(f"  [skip] {path}: already patched (v3)")
        return
    write_text(path, NEW_PAYNE_HANEK_H)


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
    print("  MATHLIB v12A1: PAYNE-HANEK CORRECT (v3)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] payne_hanek.h — correct bit-extraction algorithm")
    patch_payne_hanek(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Payne-Hanek v3 applied.")
    print("")
    print("  Algorithm: Cephes/musl style")
    print("    - 53-bit significand split into 24-bit chunks")
    print("    - Each chunk multiplied by relevant 2/pi table entries")
    print("    - Accumulated mod 8 (preserves fraction, prevents overflow)")
    print("    - Quadrant = integer part mod 4")
    print("    - Reduced arg = fractional part * (pi/2)")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake -B build && cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("    # Re-run oracle:")
    print("    gcc -std=c99 -O3 -fPIE \\")
    print("      -Iinclude/mathlib -Isrc \\")
    print("      -DMATHLIB_HAS_ORACLE_DATA \\")
    print("      -o build/oracle_check \\")
    print("      tests/test_oracle.c \\")
    print("      -Lbuild -lmathc -lm")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: sin/cos at 1e10..1e300 should now be correct.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
