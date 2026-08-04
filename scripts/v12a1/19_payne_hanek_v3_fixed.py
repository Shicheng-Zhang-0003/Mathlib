#!/usr/bin/env python3
"""
19_payne_hanek_v3_fixed.py
Run from the folder that CONTAINS the v12A1 working folder.

GIANT ERROR #1: large-argument sin/cos are garbage (up to 2.1e15 ULP).

Root cause in v2:
  ml_rem_pio2_large accumulates the FULL product q = |x| * (2/pi) in double
  precision (q_hi + q_lo), then does n = round(q), frac = q - n. For
  |x| > 1e10 the integer part of q has more digits than double can hold,
  so the fractional part (which is all we need) is destroyed.

  Additionally, the table window k_end = (E + 53) / 24 is too short: it
  drops terms whose fractional contribution is still >> 2^-53.

Fix:
  - Do NOT compute the full q. Instead, for each term, extract the integer
    contribution mod 4 into the quadrant counter, and Kahan-accumulate ONLY
    the fractional parts. The fractional accumulator stays small and precise
    regardless of how huge |x| is.
  - Widen the window to k_end = (E + 108) / 24 so all terms with fractional
    contribution >= ~2^-107 are included.

Target:
  v12A1/src/internal/payne_hanek.h

Usage:
  python3 19_payne_hanek_v3_fixed.py
  python3 19_payne_hanek_v3_fixed.py --force
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_PAYNE_HANEK_V3_FIXED"

NEW_PAYNE_HANEK_H = r"""#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H

/* MATHLIB_V12A1_PAYNE_HANEK_V3_FIXED */

#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

/* ---- Cody-Waite constants for the small-argument path ---- */
static const double
ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

/* ---- 2/pi table (66 x 24-bit chunks, standard Cephes/musl) ---- */
/* 2/pi = sum_{k=0}^{65} ml_two_over_pi[k] * 2^(-24*(k+1)) */
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
* Process one term: prod * 2^shift, where prod is an exact integer < 2^52.
*
* Adds the integer contribution (mod 4) to *n, and Kahan-accumulates the
* fractional part into *acc / *comp. This is the core of the fix: we never
* build the full product in double, so the fractional part survives no
* matter how large |x| is.
*
*   shift >= 2 : integer divisible by 4 -> contributes 0 mod 4, no fraction
*   shift == 1 : (prod*2) mod 4 = 2*(prod mod 2)
*   shift == 0 : prod mod 4
*   shift <  0 : has a fractional part -> split int / frac, accumulate frac
*/
static inline void ml_ph_process_term(
    double prod, int shift, int *n, double *acc, double *comp
) {
    if (shift >= 2) {
        return; /* multiple of 4, no quadrant or fractional contribution */
    }
    if (shift == 1) {
        double fl = (double)((long long)prod & 1LL);
        *n = (*n + (int)(2.0 * fl)) & 3;
        return;
    }
    if (shift == 0) {
        double fl = (double)((long long)prod & 3LL);
        *n = (*n + (int)fl) & 3;
        return;
    }

    /* shift < 0: fractional term */
    double t = ml_ldexp_pure(prod, shift);
    double frac_part;
    if (shift >= -52) {
        double int_part = (double)(long long)t;
        frac_part = t - int_part;
        *n = (*n + ((int)(long long)int_part & 3)) & 3;
    } else {
        frac_part = t;
    }

    /* Kahan summation of the fractional part */
    double y_k = frac_part - *comp;
    double t_k = *acc + y_k;
    *comp = (t_k - *acc) - y_k;
    *acc = t_k;
}

/*
* Payne-Hanek reduction for |x| > 1e6.
*
* Returns quadrant n in {0,1,2,3} and sets *y to the reduced argument
* in [-pi/4, pi/4].
*/
static inline int ml_rem_pio2_large(double x, double *y) {
    uint64_t bits;
    double ax = ml_fabs(x);
    int sign = (x < 0.0);
    memcpy(&bits, &ax, sizeof(uint64_t));

    int biased_e = (int)((bits >> 52) & 0x7FF);
    int E = biased_e - 1075; /* ax = m * 2^E, m is 53-bit integer */
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);

    /* Split m into high 28 bits and low 25 bits. */
    double m_hi = (double)(m >> 25);
    double m_lo = (double)(m & 0x1FFFFFFULL);

    /*
    * Table window.
    *
    * v2 used k_end = (E + 53) / 24, which is too short: it drops terms
    * whose fractional contribution is still far above 2^-53. We need all
    * terms with shift >= ~-107, i.e. k <= (E + 108) / 24.
    */
    int k_start = (E - 77) / 24;
    if (k_start < 0) k_start = 0;
    int k_end = (E + 108) / 24;
    if (k_end > 65) k_end = 65;

    int n = 0;
    double frac_acc = 0.0;
    double frac_comp = 0.0;

    for (int k = k_start; k <= k_end; k++) {
        double tk = (double)ml_two_over_pi[k];
        int base_shift = E - 24 * k - 24;

        /* m_hi term: m_hi * tk * 2^(base_shift + 25), exact (28+24 bits) */
        ml_ph_process_term(m_hi * tk, base_shift + 25,
                           &n, &frac_acc, &frac_comp);

        /* m_lo term: m_lo * tk * 2^base_shift, exact (25+24 bits) */
        ml_ph_process_term(m_lo * tk, base_shift,
                           &n, &frac_acc, &frac_comp);
    }

    /* Fold the integer part of the fractional accumulator into the quadrant. */
    double total = frac_acc;
    int extra = (int)(long long)total;
    n = (n + (extra & 3)) & 3;
    double frac = total - (double)extra;

    /* Center the fraction in [-0.5, 0.5]. */
    if (frac > 0.5) { frac -= 1.0; n = (n + 1) & 3; }
    else if (frac < -0.5) { frac += 1.0; n = (n + 3) & 3; }

    /* Reduced argument = frac * (pi/2). */
    double result = frac * ML_PH_PI2_HI + frac * ML_PH_PI2_LO;
    if (sign) { result = -result; n = (4 - n) & 3; }

    *y = result;
    return n;
}

/* ---- Unified entry point ---- */
static inline int ml_rem_pio2(double x, double *y) {
    if (ml_isnan(x) || ml_isinf(x)) { *y = ml_make_nan(); return 0; }
    double ax = ml_fabs(x);

    if (ax <= 1.0e6) {
        /* Fast 2-term Cody-Waite (unchanged, already accurate here). */
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
        print(f"  [skip] {path}: already at v3_fixed")
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
    print("  MATHLIB v12A1: PAYNE-HANEK v3_fixed (Section 19)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] payne_hanek.h — per-term int/frac separation + wide window")
    patch_payne_hanek(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Payne-Hanek v3_fixed applied.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: sin/cos at 1e10..1e300 collapse from ~1e15 ULP to <= 5.")
    print("  gamma/lgamma will STILL FAIL — that is Section 21.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
