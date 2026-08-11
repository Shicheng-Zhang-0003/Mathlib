#!/usr/bin/env python3
"""
45_a1_gamma_1step_recurrence.py

Fixes the remaining 3 oracle failures:
  gamma(0.1)   12 ULP
  gamma(0.01)  12 ULP
  gamma(0.001) 10 ULP

ROOT CAUSE:
  Direct Lanczos for x < 0.5 has cancellation factor ~3.9x,
  pushing lgamma to ~5 ULP. exp() amplifies to 10-12 ULP on gamma.

FIX:
  Use 1-step recurrence in log-space:
    lgamma(x) = lgamma(x+1) - log(x)

  For x = 0.1: lgamma(1.1) ~ -0.05, log(0.1) ~ -2.30.
  This is ADDITION (both terms become positive), not subtraction.
  Cancellation factor ~1.04 (almost zero).

  For x = 0.001: cancellation factor ~1.0001.

  The old script 31 tried this but used g/x for gamma (division rounds).
  This script uses ml_exp_dd(L) instead, staying in DD until final exp.

Targets:
  v12A1/src/integral.c
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE"


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
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/.")


def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # ----------------------------------------------------------------
    # Fix 1: Replace direct Lanczos in ml_lgamma with 1-step recurrence
    # ----------------------------------------------------------------
    lgamma_block = re.compile(
        r"/\* " + MARKER.replace("1STEP", "DIRECT_LANCZOS").replace("MATHLIB_V12A1_GAMMA_", "MATHLIB_V12A1_GAMMA_") + r" \*/\s*\n"
        r"/\* x < 0\.5: use Lanczos directly.*?\*/\s*\n"
        r"if \(x < 0\.5\) \{\s*\n"
        r"\s*ml_dd_t L = ml_lgamma_positive_dd\(x\);\s*\n"
        r"\s*return L\.hi \+ L\.lo;\s*\n"
        r"\s*\}",
        re.MULTILINE | re.DOTALL,
    )

    # Simpler: match by the code signature
    lgamma_simple = re.compile(
        r"if \(x < 0\.5\) \{\s*\n"
        r"\s*ml_dd_t L = ml_lgamma_positive_dd\(x\);\s*\n"
        r"\s*return L\.hi \+ L\.lo;\s*\n"
        r"\s*\}",
        re.MULTILINE,
    )

    lgamma_replacement = (
        "/* " + MARKER + " */\n"
        "/* x < 0.5: 1-step recurrence in log-space.\n"
        " * lgamma(x) = lgamma(x+1) - log(x)\n"
        " * For x=0.1: lgamma(1.1)~-0.05, log(0.1)~-2.30.\n"
        " * This is addition (both positive), not subtraction.\n"
        " * Cancellation factor ~1.04 (almost zero).\n"
        " */\n"
        "if (x < 0.5) {\n"
        "    ml_dd_t L = ml_lgamma_positive_dd(x + 1.0);\n"
        "    L = ml_dd_sub(L, ml_log_dd(x));\n"
        "    return L.hi + L.lo;\n"
        "}"
    )

    matches = lgamma_simple.findall(text)
    if len(matches) == 0:
        fail(
            f"{path}: could not find the lgamma x<0.5 direct Lanczos block.\n"
            f"  Expected: if (x < 0.5) {{ ml_lgamma_positive_dd(x); return L.hi + L.lo; }}"
        )
    if len(matches) > 1:
        fail(f"{path}: found {len(matches)} lgamma blocks, expected 1.")

    text = lgamma_simple.sub(lgamma_replacement, text, count=1)
    print("  [ok] Replaced lgamma direct Lanczos with 1-step recurrence")

    # ----------------------------------------------------------------
    # Fix 2: Replace direct Lanczos in ml_gamma_new with 1-step recurrence
    # ----------------------------------------------------------------
    gamma_simple = re.compile(
        r"if \(x < 0\.5\) \{\s*\n"
        r"\s*return ml_gamma_positive\(x\);\s*\n"
        r"\s*\}",
        re.MULTILINE,
    )

    gamma_replacement = (
        "/* " + MARKER + " */\n"
        "/* x < 0.5: 1-step recurrence in log-space, then exp.\n"
        " * Same cancellation fix as ml_lgamma above.\n"
        " * Uses ml_exp_dd(L) instead of g/x to avoid division rounding.\n"
        " */\n"
        "if (x < 0.5) {\n"
        "    ml_dd_t L = ml_lgamma_positive_dd(x + 1.0);\n"
        "    L = ml_dd_sub(L, ml_log_dd(x));\n"
        "    return ml_exp_dd(L);\n"
        "}"
    )

    matches = gamma_simple.findall(text)
    if len(matches) == 0:
        fail(
            f"{path}: could not find the gamma_new x<0.5 direct Lanczos block.\n"
            f"  Expected: if (x < 0.5) {{ return ml_gamma_positive(x); }}"
        )
    if len(matches) > 1:
        fail(f"{path}: found {len(matches)} gamma blocks, expected 1.")

    text = gamma_simple.sub(gamma_replacement, text, count=1)
    print("  [ok] Replaced gamma_new direct Lanczos with 1-step recurrence")

    write_text(path, text)


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
    print("  MATHLIB v12A1: 1-STEP RECURRENCE (kill the last ULPs)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — replace direct Lanczos with 1-step recurrence")
    patch_integral(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  ROOT CAUSE:")
    print("    Direct Lanczos for x<0.5 has cancellation factor ~3.9x.")
    print("    This pushes lgamma to ~5 ULP, which exp() amplifies")
    print("    to 10-12 ULP on gamma.")
    print("")
    print("  FIX:")
    print("    lgamma(x) = lgamma(x+1) - log(x)")
    print("    For x=0.1: lgamma(1.1)~-0.05, log(0.1)~-2.30.")
    print("    This is ADDITION, not subtraction.")
    print("    Cancellation factor drops to ~1.04.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: 211 passed, 0 failed.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
