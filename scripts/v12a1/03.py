#!/usr/bin/env python3
"""
03_exp_cody_waite.py

Run from the folder that CONTAINS the v12A1 working folder.

Fixes the Cody-Waite reduction in ml_exp().

Problem:
  The current code computes the reduced argument as:

    double r = x - n * 0.69314718036912381649
                 - n * 1.90821490974462528503e-10;

  This is TWO separate rounded operations:
    1. x - n * hi   (multiply rounds, subtract rounds)
    2. result - n * lo   (multiply rounds, subtract rounds)

  Each step loses bits. The residual r can be off by up to 2 ULP
  from the true value x - n * ln(2).

Fix:
  Use ML_FMA (hardware fused multiply-add) for each step:

    double r = ML_FMA(-n, 0.69314718036912381649, x);
    r = ML_FMA(-n, 1.90821490974462528503e-10, r);

  Each FMA rounds ONCE. The residual is exact to 0.5 ULP.

  The 2-term split of ln(2) gives ~106 bits of precision.
  For |n| <= 1024 (the exp overflow limit), this is more than
  sufficient. No 3-term split needed.

Targets:
  v12A1/src/exp_log.c

Does NOT change:
  - The polynomial evaluation (still Taylor, swapped in script 06)
  - The overflow/underflow thresholds
  - The ml_ldexp_pure reconstruction
  - Any other function

Usage:
    python3 03_exp_cody_waite.py
    python3 03_exp_cody_waite.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_EXP_FMA_REDUCTION"

OLD_REDUCTION = (
    "    double r = x - n * 0.69314718036912381649"
    " - n * 1.90821490974462528503e-10;"
)

NEW_REDUCTION = (
    "    /* " + MARKER + " */\n"
    "    /*\n"
    "     * Error-free Cody-Waite reduction.\n"
    "     *\n"
    "     * The old code used two separate rounded subtractions:\n"
    "     *   r = x - n * hi - n * lo\n"
    "     *\n"
    "     * Each multiply-subtract pair rounds twice, losing up to\n"
    "     * 2 ULP of precision in the residual.\n"
    "     *\n"
    "     * FMA rounds once per step. The 2-term split of ln(2)\n"
    "     * provides ~106 bits, which is sufficient for |n| <= 1024.\n"
    "     */\n"
    "    double r = ML_FMA(-n, 0.69314718036912381649, x);\n"
    "    r = ML_FMA(-n, 1.90821490974462528503e-10, r);"
)


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
    if (root / "src" / "exp_log.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_exp_reduction(v12: Path, force: bool) -> None:
    path = v12 / "src" / "exp_log.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Try exact match first
    if OLD_REDUCTION in text:
        text = text.replace(OLD_REDUCTION, NEW_REDUCTION, 1)
        write_text(path, text)
        return

    # Try regex for whitespace variation
    pattern = re.compile(
        r"double\s+r\s*=\s*x\s*-\s*n\s*\*\s*0\.69314718036912381649\s*"
        r"-\s*n\s*\*\s*1\.90821490974462528503e-10\s*;"
    )
    patched, count = pattern.subn(NEW_REDUCTION, text, count=1)
    if count != 1:
        fail(
            f"{path}: could not find the Cody-Waite reduction line. "
            "Source may have drifted."
        )
    write_text(path, patched)


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
    print("  MATHLIB v12A1: EXP CODY-WAITE FMA FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] exp_log.c reduction")
    patch_exp_reduction(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  What changed:")
    print("    ml_exp() Cody-Waite reduction now uses ML_FMA")
    print("    Old: r = x - n*hi - n*lo  (4 rounding ops)")
    print("    New: r = FMA(-n,hi,x); r = FMA(-n,lo,r)  (2 rounding ops)")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    grep -n 'ML_FMA(-n' src/exp_log.c")
    print("    # Should show 2 hits inside ml_exp")
    print("")
    print("    cmake -B build && cmake --build build")
    print("    ./build/test_core")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("  Next: 04_log_reconstruction.py")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
