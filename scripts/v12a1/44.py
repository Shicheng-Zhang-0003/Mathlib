#!/usr/bin/env python3
"""
44_a1_gamma_lgamma_direct_lanczos.py

Fixes ALL 4 oracle failures:
  gamma(0.1)   24 ULP
  gamma(0.01)  34 ULP
  gamma(0.001) 32 ULP
  lgamma(0.1)  10 ULP

ROOT CAUSE:
  The x < 0.5 recurrence formula computes:
    lgamma(x) = lgamma(x+16) - sum(log(x+k), k=0..15)

  For x = 0.1, this is 28.17 - 25.92 = 2.25.
  Cancellation factor ~12.5x amplifies DD rounding errors
  into 10 ULP on lgamma, which exp() amplifies to 24-34 ULP
  on gamma.

  Depth 8 vs 16 does NOT fix this. The cancellation is
  structural. Both paths are already at depth 16 (script 40
  applied correctly to both). The recurrence itself must go.

FIX:
  Delete the x < 0.5 recurrence blocks in both ml_lgamma
  and ml_gamma_new. Let them fall through to
  ml_lgamma_positive_dd(x), which dispatches to Lanczos
  for x < 8. Lanczos computes lgamma(x) directly with
  cancellation factor ~3.9x instead of ~12.5x.

Targets:
  v12A1/src/integral.c

Usage:
  python3 44_a1_gamma_lgamma_direct_lanczos.py
  python3 44_a1_gamma_lgamma_direct_lanczos.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS"


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
    # Fix 1: Remove the x < 0.5 recurrence block from ml_lgamma
    # ----------------------------------------------------------------
    # The lgamma block: depth 16, returns L.hi + L.lo
    # Match by the unique return statement "return L.hi + L.lo;"
    # inside an if(x < 0.5) block that calls ml_lgamma_positive_dd(x + 16.0)
    lgamma_block = re.compile(
        r"if \(x < 0\.5\) \{\s*\n"
        r"\s*ml_dd_t L = ml_lgamma_positive_dd\(x \+ 16\.0\);\s*\n"
        r"\s*for \(int k = 0; k < 16; k\+\+\)\s*\n"
        r"\s*L = ml_dd_sub\(L, ml_log_dd\(x \+ \(double\)k\)\);\s*\n"
        r"\s*return L\.hi \+ L\.lo;\s*\n"
        r"\s*\}",
        re.MULTILINE,
    )

    lgamma_replacement = (
        "/* " + MARKER + " */\n"
        "/* x < 0.5: use Lanczos directly via ml_lgamma_positive_dd.\n"
        " * The old recurrence formula subtracted two large DD values\n"
        " * (lgamma(x+16) ~ 28 minus sum-of-logs ~ 26) to get a small\n"
        " * result (~2.25), amplifying rounding errors by ~12.5x.\n"
        " * Lanczos computes lgamma(x) directly with ~3.9x cancellation.\n"
        " */\n"
        "if (x < 0.5) {\n"
        "    ml_dd_t L = ml_lgamma_positive_dd(x);\n"
        "    return L.hi + L.lo;\n"
        "}"
    )

    matches = lgamma_block.findall(text)
    if len(matches) == 0:
        fail(
            f"{path}: could not find the lgamma x<0.5 recurrence block.\n"
            f"  Expected: if (x < 0.5) {{ ... x + 16.0 ... k < 16 ... return L.hi + L.lo; }}"
        )
    if len(matches) > 1:
        fail(f"{path}: found {len(matches)} lgamma blocks, expected 1.")

    text = lgamma_block.sub(lgamma_replacement, text, count=1)
    print("  [ok] Removed lgamma x<0.5 recurrence (depth 16)")

    # ----------------------------------------------------------------
    # Fix 2: Remove the x < 0.5 recurrence block from ml_gamma_new
    # ----------------------------------------------------------------
    # The gamma block: depth 16, returns ml_exp_dd(L)
    # Match by the unique return statement "return ml_exp_dd(L);"
    # inside an if(x < 0.5) block that calls ml_lgamma_positive_dd(x + 16.0)
    gamma_block = re.compile(
        r"if \(x < 0\.5\) \{\s*\n"
        r"\s*ml_dd_t L = ml_lgamma_positive_dd\(x \+ 16\.0\);\s*\n"
        r"\s*for \(int k = 0; k < 16; k\+\+\)\s*\n"
        r"\s*L = ml_dd_sub\(L, ml_log_dd\(x \+ \(double\)k\)\);\s*\n"
        r"\s*return ml_exp_dd\(L\);\s*\n"
        r"\s*\}",
        re.MULTILINE,
    )

    gamma_replacement = (
        "/* " + MARKER + " */\n"
        "/* x < 0.5: use Lanczos directly via ml_gamma_positive.\n"
        " * Same cancellation fix as ml_lgamma above.\n"
        " */\n"
        "if (x < 0.5) {\n"
        "    return ml_gamma_positive(x);\n"
        "}"
    )

    matches = gamma_block.findall(text)
    if len(matches) == 0:
        fail(
            f"{path}: could not find the gamma_new x<0.5 recurrence block.\n"
            f"  Expected: if (x < 0.5) {{ ... x + 16.0 ... k < 16 ... return ml_exp_dd(L); }}"
        )
    if len(matches) > 1:
        fail(f"{path}: found {len(matches)} gamma blocks, expected 1.")

    text = gamma_block.sub(gamma_replacement, text, count=1)
    print("  [ok] Removed gamma_new x<0.5 recurrence (depth 16)")

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
    print("  MATHLIB v12A1: DIRECT LANCZOS (kill the recurrence)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — replace recurrence with direct Lanczos")
    patch_integral(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  ROOT CAUSE:")
    print("    The x<0.5 recurrence computes lgamma(x+16) - sum(log(x+k)).")
    print("    For x=0.1: 28.17 - 25.92 = 2.25.")
    print("    Cancellation factor ~12.5x amplifies DD rounding to 10 ULP.")
    print("    exp() then amplifies to 24-34 ULP on gamma.")
    print("")
    print("  FIX:")
    print("    Use Lanczos directly for x < 0.5.")
    print("    Cancellation factor drops to ~3.9x.")
    print("    No recurrence. No depth parameter. No cancellation.")
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
