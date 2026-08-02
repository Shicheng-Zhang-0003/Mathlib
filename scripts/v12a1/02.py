#!/usr/bin/env python3
"""
02_error_free_cleanup.py

Run from the folder that CONTAINS the v12A1 working folder.

Fixes the ml_fma / ML_FMA naming collision in error_free.h.

Problem:
  error_free.h defines ml_fma() — a software FMA via Two-Product + Two-Sum.
  ml_compiler.h defines ML_FMA — a macro routing to hardware __builtin_fma.

  The software version rounds TWICE (up to 2 ULP error).
  The hardware version rounds ONCE (0.5 ULP, correctly rounded).

  Having both named nearly identically is a latent confusion bomb.
  The actual math code (exp_log.c, minimax.h) correctly uses ML_FMA.
  Only one test (fuzz_god_mode.c) calls the software ml_fma directly.

This script:
  1. Renames ml_fma -> ml_fma_soft in error_free.h
  2. Adds documentation clarifying the distinction
  3. Adds ml_fast_two_sum alias with documented precondition
  4. Updates the call site in fuzz_god_mode.c

Does NOT modify any numerical behavior.
Does NOT touch exp_log.c, minimax.h, or any math source.

Usage:
    python3 02_error_free_cleanup.py
    python3 02_error_free_cleanup.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_ERROR_FREE_CLEANUP"

NEW_ERROR_FREE_H = r"""#ifndef LIBMATHC_ERROR_FREE_H
#define LIBMATHC_ERROR_FREE_H

/* MATHLIB_V12A1_ERROR_FREE_CLEANUP */
/*
 * ERROR-FREE TRANSFORMATIONS
 *
 * These functions compute exact rounding errors from floating-point
 * operations. They are building blocks for compensated algorithms
 * (Kahan summation, compensated Horner, etc.).
 *
 * NAMING CONVENTION:
 *   ml_two_sum(a, b, &err)      -> s = fl(a + b), err = (a + b) - s exactly
 *   ml_fast_two_sum(a, b, &err) -> same, but REQUIRES |a| >= |b|
 *   ml_two_product(a, b, &err)  -> p = fl(a * b), err = (a * b) - p exactly
 *   ml_fma_soft(a, b, c)        -> software FMA via Two-Product + Two-Sum
 *
 * IMPORTANT: ml_fma_soft is NOT the same as ML_FMA (from ml_compiler.h).
 *
 *   ML_FMA(a, b, c)    = hardware fused multiply-add, rounds ONCE (0.5 ULP)
 *   ml_fma_soft(a, b, c) = software emulation, rounds TWICE (up to 2 ULP)
 *
 * Use ML_FMA for all production math. Use ml_fma_soft only when you
 * explicitly need the software path (e.g., testing, or platforms
 * without hardware FMA).
 */

#include "ml_compiler.h"

/*
 * Knuth's Two-Sum.
 * No magnitude assumption on a, b.
 * Returns s = fl(a + b) and sets *err = (a + b) - s exactly.
 */
static inline double ml_two_sum(double a, double b, double *err) {
    double s = a + b;
    double v = s - a;
    *err = (a - (s - v)) + (b - v);
    return s;
}

/*
 * Dekker's Fast Two-Sum.
 * PRECONDITION: |a| >= |b| (caller must guarantee this).
 * Returns s = fl(a + b) and sets *err = (a + b) - s exactly.
 *
 * Faster than ml_two_sum (3 ops vs 6) but only valid when
 * the magnitude ordering is known.
 */
static inline double ml_fast_two_sum(double a, double b, double *err) {
    double s = a + b;
    double z = s - a;
    *err = b - z;
    return s;
}

/*
 * Dekker's Two-Product.
 * Returns p = fl(a * b) and sets *err = (a * b) - p exactly.
 *
 * Uses Dekker splitting (multiply by 2^26 + 1) to split each
 * operand into high and low 26-bit halves. This is exact for
 * any finite double because the significand is 53 bits.
 *
 * Note: If hardware FMA is available, you can compute the error
 * term more cheaply as: err = ML_FMA(a, b, -p). But this function
 * remains useful for platforms without FMA and for code that
 * needs to be FMA-independent.
 */
static inline double ml_two_product(double a, double b, double *err) {
    double p = a * b;
    double ca = a * 67108865.0; /* 2^26 + 1 */
    double a_hi = ca - (ca - a);
    double a_lo = a - a_hi;
    double cb = b * 67108865.0;
    double b_hi = cb - (cb - b);
    double b_lo = b - b_hi;
    *err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    return p;
}

/*
 * Software FMA emulation.
 *
 * Computes fl(a * b + c) using error-free transformations.
 * This rounds TWICE and has up to 2 ULP error.
 *
 * DO NOT USE THIS IN PRODUCTION MATH. Use ML_FMA instead.
 * This exists for:
 *   - testing error-free transformation correctness
 *   - platforms without hardware FMA (the ML_FMA macro falls
 *     back to (a*b)+c on such platforms, which is even worse)
 *   - educational/reference purposes
 */
static inline double ml_fma_soft(double a, double b, double c) {
    double p, prod_err;
    p = ml_two_product(a, b, &prod_err);
    double s1, sum_err;
    s1 = ml_two_sum(p, c, &sum_err);
    return s1 + (prod_err + sum_err);
}

#endif /* LIBMATHC_ERROR_FREE_H */
"""

def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


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
    if (root / "src" / "internal" / "error_free.h").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_error_free(v12: Path, force: bool) -> None:
    path = v12 / "src" / "internal" / "error_free.h"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    write_text(path, NEW_ERROR_FREE_H)


def patch_fuzzer(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "fuzz_god_mode.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if "ml_fma_soft" in text and not force:
        print(f"  [skip] {path}: already uses ml_fma_soft")
        return
    # Replace the ml_fma call in test_catastrophic_cancellation
    old = "double fma_res = ml_fma(1e16, 1.0, 1.0);"
    new = "double fma_res = ml_fma_soft(1e16, 1.0, 1.0);"
    if old not in text:
        # Try a more lenient match
        pattern = re.compile(r"double\s+fma_res\s*=\s*ml_fma\s*\(")
        if pattern.search(text):
            text = pattern.sub("double fma_res = ml_fma_soft(", text, count=1)
            write_text(path, text)
            print(f"  [patch] {path}: ml_fma -> ml_fma_soft (regex)")
            return
        fail(f"{path}: could not find ml_fma call site. Source may have drifted.")
    text = text.replace(old, new, 1)
    write_text(path, text)


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
    print("  MATHLIB v12A1: ERROR-FREE TRANSFORMATION CLEANUP")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] error_free.h")
    patch_error_free(v12, force)

    print("\n[2/2] fuzz_god_mode.c call site")
    patch_fuzzer(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Cleanup complete.")
    print("")
    print("  What changed:")
    print("    - ml_fma renamed to ml_fma_soft in error_free.h")
    print("    - Documentation added explaining ML_FMA vs ml_fma_soft")
    print("    - ml_fast_two_sum precondition documented")
    print("    - fuzz_god_mode.c updated to call ml_fma_soft")
    print("")
    print("  What did NOT change:")
    print("    - No numerical behavior modified")
    print("    - exp_log.c still uses ML_FMA (hardware FMA)")
    print("    - minimax.h still uses ML_FMA (hardware FMA)")
    print("    - All existing tests should still pass")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    grep -rn 'ml_fma(' src/ tests/ include/")
    print("    # Should return ZERO hits (all renamed to ml_fma_soft)")
    print("")
    print("    grep -rn 'ML_FMA(' src/")
    print("    # Should show exp_log.c and minimax.h using hardware FMA")
    print("")
    print("  Next: 03_exp_cody_waite.py")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
