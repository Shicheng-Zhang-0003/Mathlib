#!/usr/bin/env python3
"""
04_log_reconstruction.py

Run from the folder that CONTAINS the v12A1 working folder.

Fixes the ml_log() final reconstruction.

Problem:
  The current return statement is:

    return ML_FMA((double)e, ML_LN2, z * poly);

  This single FMA rounds e * ln(2) to 53 bits before adding z * poly.
  For |e| up to ~1024, the rounding error in e * ln(2) is up to
  1024 * 2^-53 ~ 1.1e-13, which is several ULP of the final result.

Fix:
  Split ln(2) into high and low parts (same technique already used
  in ml_exp after script 03). Use two FMAs to capture the low bits:

    return ML_FMA((double)e, ML_LN2_HI, z * poly)
         + (double)e * ML_LN2_LO;

  The first FMA computes e * ln2_hi + z * poly with one rounding.
  The second term adds back the bits lost by truncating ln(2).

  Together, ln2_hi + ln2_lo represent ln(2) to ~106 bits.

Targets:
    v12A1/src/exp_log.c

Does NOT change:
  - The polynomial evaluation
  - The frexp / mantissa adjustment
  - The z = (m-1)/(m+1) transform
  - Any other function

Usage:
    python3 04_log_reconstruction.py
    python3 04_log_reconstruction.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_LOG_COMPENSATED_RECONSTRUCT"

LN2_SPLIT_BLOCK = r"""/* MATHLIB_V12A1_LOG_COMPENSATED_RECONSTRUCT */
/*
 * Split ln(2) into high and low parts for compensated reconstruction.
 *
 * ML_LN2_HI has the low 26 bits of its significand zeroed.
 * ML_LN2_LO captures the remaining bits.
 * Together they represent ln(2) to ~106 bits of precision.
 *
 * These are the same values used by musl libc and match the
 * 2-term Cody-Waite split already used in ml_exp (script 03).
 */
#ifndef ML_LN2_HI
#define ML_LN2_HI 6.93147180369123816490e-01
#endif
#ifndef ML_LN2_LO
#define ML_LN2_LO 1.90821492927058500170e-10
#endif
"""

OLD_RETURN = (
    "    /* MATHLIB_CLOSURE_P2_LOG_FMA_RECONSTRUCT */\n"
    "    return ML_FMA((double)e, ML_LN2, z * poly);"
)

NEW_RETURN = (
    "    /* MATHLIB_V12A1_LOG_COMPENSATED_RECONSTRUCT */\n"
    "    /*\n"
    "     * Compensated reconstruction.\n"
    "     *\n"
    "     * Old: ML_FMA((double)e, ML_LN2, z * poly)\n"
    "     *   -> rounds e * ln(2) to 53 bits, losing up to 1e-13\n"
    "     *      for large |e|.\n"
    "     *\n"
    "     * New: split ln(2) and use two terms.\n"
    "     *   -> FMA(e, ln2_hi, z*poly) captures the high product\n"
    "     *      with one rounding.\n"
    "     *   -> e * ln2_lo adds back the truncated low bits.\n"
    "     */\n"
    "    return ML_FMA((double)e, ML_LN2_HI, z * poly)\n"
    "         + (double)e * ML_LN2_LO;"
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


def patch_log_reconstruction(v12: Path, force: bool) -> None:
    path = v12 / "src" / "exp_log.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Step 1: Insert the ln(2) split defines after the existing macro block.
    # Anchor on the last existing #endif in the macro section.
    if "ML_LN2_HI" not in text:
        # Find the ML_LOG_UNDERFLOW block and insert after its #endif
        anchor = re.compile(
            r"(#ifndef ML_LOG_UNDERFLOW\s*\n"
            r"#define ML_LOG_UNDERFLOW[^\n]*\n"
            r"#endif\s*\n)"
        )
        patched, count = anchor.subn(
            lambda m: m.group(1) + "\n" + LN2_SPLIT_BLOCK,
            text,
            count=1,
        )
        if count != 1:
            # Fallback: insert after the first #include block
            anchor2 = re.compile(r'(#include "internal/pow_util\.h"[^\n]*\n)')
            patched, count = anchor2.subn(
                lambda m: m.group(1) + "\n" + LN2_SPLIT_BLOCK,
                text,
                count=1,
            )
            if count != 1:
                fail(f"{path}: could not find anchor for ln(2) split defines.")
        text = patched
        print(f"  [patch] {path}: inserted ML_LN2_HI / ML_LN2_LO defines")

    # Step 2: Replace the return statement.
    if OLD_RETURN in text:
        text = text.replace(OLD_RETURN, NEW_RETURN, 1)
        write_text(path, text)
        return

    # Regex fallback for whitespace variation
    pattern = re.compile(
        r"/\* MATHLIB_CLOSURE_P2_LOG_FMA_RECONSTRUCT \*/\s*\n"
        r"\s*return ML_FMA\(\(double\)e,\s*ML_LN2,\s*z \* poly\);"
    )
    patched, count = pattern.subn(NEW_RETURN, text, count=1)
    if count != 1:
        fail(
            f"{path}: could not find ml_log return statement. "
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
    print("  MATHLIB v12A1: LOG COMPENSATED RECONSTRUCTION")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] exp_log.c — ml_log return")
    patch_log_reconstruction(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  What changed:")
    print("    Old: return ML_FMA(e, ML_LN2, z * poly)")
    print("    New: return ML_FMA(e, ML_LN2_HI, z * poly)")
    print("              + e * ML_LN2_LO")
    print("")
    print("  Why:")
    print("    ML_LN2 is 53 bits. e * ML_LN2 rounds, losing up to")
    print("    1024 * 2^-53 ~ 1.1e-13 for large exponents.")
    print("    ML_LN2_HI + ML_LN2_LO gives ~106 bits.")
    print("    The second term recovers the truncated low bits.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake -B build && cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("  Next: 09_pow_extended.py (now unblocked)")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
