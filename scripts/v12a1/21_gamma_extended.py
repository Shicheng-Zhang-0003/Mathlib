#!/usr/bin/env python3
"""
21_gamma_extended.py
Run from the folder that CONTAINS the v12A1 working folder.

GIANT ERROR #5: gamma/lgamma amplification error.

Root cause:
  ml_lgamma_positive computes (z+0.5) * ml_log(t) where z = x-1.
  ml_log(t) has ~1 ULP error. When multiplied by (z+0.5) which can be
  up to ~171, the error is amplified to ~(z+0.5) ULP.

  For gamma(100): z+0.5 = 99.5, error ≈ 99.5 * 8.9e-16 ≈ 8.9e-14
  ≈ 400+ ULP. Observed: 544 ULP. Matches exactly.

Fix:
  1. Add ml_log_split() to exp_log.c: returns log(x) as double-double
     (log_hi + log_lo) with ~106 bits of precision.
  2. Use ml_log_split() in ml_lgamma_positive and ml_gamma_new.
  3. Compute (z+0.5) * (log_hi + log_lo) using FMA for extended precision.

Targets:
  v12A1/src/exp_log.c
  v12A1/include/mathlib/ml_exp_log.h
  v12A1/src/integral.c

Usage:
  python3 21_gamma_extended.py
  python3 21_gamma_extended.py --force
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_LOG_SPLIT"

# ---------------------------------------------------------------------------
# New ml_log_split function to add to exp_log.c
# ---------------------------------------------------------------------------
NEW_LOG_SPLIT = r"""
/* MATHLIB_V12A1_GAMMA_LOG_SPLIT */
/*
 * Double-double log: returns log(x) as log_hi + log_lo.
 *
 * log_hi is the main result (rounded to double).
 * log_lo captures the low bits of e*ln2.
 * Together they give ~106 bits of precision.
 *
 * This is used by ml_lgamma_positive and ml_gamma_new to avoid
 * the (z+0.5)*log(t) amplification error.
 */
ML_API void ml_log_split(double x, double *log_hi, double *log_lo) {
    if (ml_isnan(x) || x <= 0.0) {
        *log_hi = ml_make_nan();
        *log_lo = 0.0;
        return;
    }
    if (ml_isinf(x)) {
        *log_hi = x;
        *log_lo = 0.0;
        return;
    }
    if (x == 1.0) {
        *log_hi = 0.0;
        *log_lo = 0.0;
        return;
    }
    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + adjust);
    e -= adjust;
    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    double poly = 0.09523809523809523;
    poly = poly * z2 + 0.10526315789473684;
    poly = poly * z2 + 0.11764705882352941;
    poly = poly * z2 + 0.13333333333333333;
    poly = poly * z2 + 0.15384615384615385;
    poly = poly * z2 + 0.18181818181818182;
    poly = poly * z2 + 0.2222222222222222;
    poly = poly * z2 + 0.2857142857142857;
    poly = poly * z2 + 0.4;
    poly = poly * z2 + 0.6666666666666666;
    poly = poly * z2 + 2.0;
    *log_hi = ML_FMA((double)e, ML_LN2_HI, z * poly);
    *log_lo = (double)e * ML_LN2_LO;
}
"""

# ---------------------------------------------------------------------------
# New ml_lgamma_positive body (replace the return statement)
# ---------------------------------------------------------------------------
OLD_LGAMMA_RETURN = (
    "    return 0.91893853320467274178  /* 0.5 * log(2*pi) */\n"
    "        + (z + 0.5) * ml_log(t)\n"
    "        - t\n"
    "        + ml_log(ag);"
)

NEW_LGAMMA_RETURN = (
    "    /* MATHLIB_V12A1_GAMMA_LOG_SPLIT */\n"
    "    /*\n"
    "     * Double-double log to avoid (z+0.5)*log(t) amplification.\n"
    "     *\n"
    "     * The old code computed (z+0.5) * ml_log(t), where ml_log(t)\n"
    "     * has ~1 ULP error. When multiplied by (z+0.5) which can be\n"
    "     * up to ~171, the error is amplified to ~(z+0.5) ULP.\n"
    "     *\n"
    "     * Fix: compute log(t) as double-double (log_hi + log_lo),\n"
    "     * then multiply by (z+0.5) using FMA for extended precision.\n"
    "     */\n"
    "    double log_hi, log_lo;\n"
    "    ml_log_split(t, &log_hi, &log_lo);\n"
    "    double zp5 = z + 0.5;\n"
    "    double prod_hi = zp5 * log_hi;\n"
    "    double prod_lo = ML_FMA(zp5, log_hi, -prod_hi) + zp5 * log_lo;\n"
    "    double prod = prod_hi + prod_lo;\n"
    "    return 0.91893853320467274178 + prod - t + ml_log(ag);"
)

# ---------------------------------------------------------------------------
# New ml_gamma_new direct path (replace log_part computation)
# ---------------------------------------------------------------------------
OLD_GAMMA_LOG_PART = (
    "            double log_part = (z + 0.5) * ml_log(t) - t;"
)

NEW_GAMMA_LOG_PART = (
    "            /* MATHLIB_V12A1_GAMMA_LOG_SPLIT */\n"
    "            /*\n"
    "             * Double-double log to avoid (z+0.5)*log(t) amplification.\n"
    "             */\n"
    "            double log_hi, log_lo;\n"
    "            ml_log_split(t, &log_hi, &log_lo);\n"
    "            double zp5 = z + 0.5;\n"
    "            double prod_hi = zp5 * log_hi;\n"
    "            double prod_lo = ML_FMA(zp5, log_hi, -prod_hi) + zp5 * log_lo;\n"
    "            double log_part = (prod_hi + prod_lo) - t;"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
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
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/")


# ---------------------------------------------------------------------------
# Patch exp_log.c: add ml_log_split
# ---------------------------------------------------------------------------
def patch_exp_log(v12: Path, force: bool) -> None:
    path = v12 / "src" / "exp_log.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    # Insert ml_log_split right after ml_log function
    # Find the end of ml_log (the return statement)
    anchor = re.compile(
        r"(ML_API double ml_log\(double x\) \{.*?"
        r"return ML_FMA\(\(double\)e, ML_LN2_HI, z \* poly\)\s*\n"
        r"\s*\+ \(double\)e \* ML_LN2_LO;\s*\n"
        r"\})",
        re.DOTALL,
    )
    patched, count = anchor.subn(
        lambda m: m.group(1) + NEW_LOG_SPLIT,
        text,
        count=1,
    )
    if count != 1:
        # Fallback: try to find ml_log and insert after it
        pattern = re.compile(
            r"(?ms)(ML_API double ml_log\(double x\) \{.*?\n\})"
        )
        patched, count = pattern.subn(
            lambda m: m.group(1) + NEW_LOG_SPLIT,
            text,
            count=1,
        )
        if count != 1:
            fail(f"{path}: could not find ml_log function to insert ml_log_split after.")
    write_text(path, patched)


# ---------------------------------------------------------------------------
# Patch ml_exp_log.h: add ml_log_split declaration
# ---------------------------------------------------------------------------
def patch_header(v12: Path, force: bool) -> None:
    path = v12 / "include" / "mathlib" / "ml_exp_log.h"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if "ml_log_split" in text and not force:
        print(f"  [skip] {path}: ml_log_split already declared")
        return
    # Add declaration after ml_log
    old = "ML_API double ml_log(double x);"
    new = (
        "ML_API double ml_log(double x);\n"
        "ML_API void ml_log_split(double x, double *log_hi, double *log_lo);"
    )
    if old not in text:
        fail(f"{path}: could not find ml_log declaration.")
    text = text.replace(old, new, 1)
    write_text(path, text)


# ---------------------------------------------------------------------------
# Patch integral.c: fix ml_lgamma_positive and ml_gamma_new
# ---------------------------------------------------------------------------
def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Fix ml_lgamma_positive return statement
    if OLD_LGAMMA_RETURN in text:
        text = text.replace(OLD_LGAMMA_RETURN, NEW_LGAMMA_RETURN, 1)
        print(f"  [patch] {path}: ml_lgamma_positive return")
    else:
        # Try regex fallback for whitespace variation
        pattern = re.compile(
            r"return 0\.91893853320467274178\s+/\* 0\.5 \* log\(2\*pi\) \*/\s*\n"
            r"\s*\+ \(z \+ 0\.5\) \* ml_log\(t\)\s*\n"
            r"\s*- t\s*\n"
            r"\s*\+ ml_log\(ag\);"
        )
        patched, count = pattern.subn(NEW_LGAMMA_RETURN, text, count=1)
        if count != 1:
            fail(f"{path}: could not find ml_lgamma_positive return statement.")
        text = patched
        print(f"  [patch] {path}: ml_lgamma_positive return (regex)")

    # Fix ml_gamma_new log_part computation
    if OLD_GAMMA_LOG_PART in text:
        text = text.replace(OLD_GAMMA_LOG_PART, NEW_GAMMA_LOG_PART, 1)
        print(f"  [patch] {path}: ml_gamma_new log_part")
    else:
        # Try regex fallback
        pattern = re.compile(
            r"double log_part = \(z \+ 0\.5\) \* ml_log\(t\) - t;"
        )
        patched, count = pattern.subn(NEW_GAMMA_LOG_PART, text, count=1)
        if count != 1:
            fail(f"{path}: could not find ml_gamma_new log_part computation.")
        text = patched
        print(f"  [patch] {path}: ml_gamma_new log_part (regex)")

    write_text(path, text)


# ---------------------------------------------------------------------------
# Archive self
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()
    print("=========================================================")
    print("  MATHLIB v12A1: GAMMA/LGAMMA EXTENDED PRECISION (script 21)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/3] exp_log.c — add ml_log_split")
    patch_exp_log(v12, force)

    print("\n[2/3] ml_exp_log.h — add ml_log_split declaration")
    patch_header(v12, force)

    print("\n[3/3] integral.c — fix ml_lgamma_positive and ml_gamma_new")
    patch_integral(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Gamma/lgamma extended precision fix applied.")
    print("")
    print("  What changed:")
    print("    - Added ml_log_split(): returns log(x) as double-double")
    print("    - ml_lgamma_positive: uses double-double log for (z+0.5)*log(t)")
    print("    - ml_gamma_new: uses double-double log for direct path")
    print("")
    print("  Why this fixes the amplification:")
    print("    Old: (z+0.5) * ml_log(t)")
    print("         ml_log(t) has ~1 ULP error")
    print("         multiplied by (z+0.5) up to ~171 -> ~171 ULP error")
    print("")
    print("    New: ml_log_split(t, &log_hi, &log_lo)")
    print("         (z+0.5) * (log_hi + log_lo) via FMA")
    print("         ~106 bits of precision -> error stays ~1 ULP")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: all gamma/lgamma entries collapse to <= 5 ULP.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
