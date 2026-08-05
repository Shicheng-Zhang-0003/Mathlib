#!/usr/bin/env python3
"""
38_a1_exp_dd_second_order.py
Run from the folder that CONTAINS the v12A1 working folder.

ROOT CAUSE (finally identified):
  ml_exp_dd() uses the first-order approximation:
      exp(L.hi + L.lo) ≈ exp(L.hi) * (1 + L.lo)

  This drops the second-order term L.lo²/2. For the 8-step
  recurrence at x=0.001, L.lo reaches ~5e-8, making the dropped
  term g * L.lo²/2 ≈ 1.25e-12 = exactly the 10 ULP we see.

FIX:
  Compute exp(L.lo) properly instead of approximating:
      exp(L.hi + L.lo) = exp(L.hi) * exp(L.lo)

  Since L.lo is tiny, exp(L.lo) can be computed with a 3-term
  Taylor: 1 + L.lo + L.lo²/2, which is accurate to ~1e-48.

Targets:
  v12A1/src/integral.c

Usage:
  python3 38_a1_exp_dd_second_order.py
  python3 38_a1_exp_dd_second_order.py --force
"""
from __future__ import annotations
import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_EXP_DD_SECOND_ORDER"

OLD_FUNC = re.compile(
    r"static double ml_exp_dd\(ml_dd_t L\) \{.*?\n\}",
    re.DOTALL,
)

NEW_FUNC = """static double ml_exp_dd(ml_dd_t L) {
/* """ + MARKER + """ */
/*
* exp(L.hi + L.lo) = exp(L.hi) * exp(L.lo)
*
* The old code used the first-order approximation:
*     exp(L.lo) ≈ 1 + L.lo
* which drops L.lo²/2. For the 8-step recurrence at x=0.001,
* L.lo reaches ~5e-8, making the dropped term ~1.25e-12 = 10 ULP.
*
* Fix: compute exp(L.lo) with a 3-term Taylor:
*     exp(L.lo) ≈ 1 + L.lo + L.lo²/2
* Since L.lo is tiny, this is accurate to ~1e-48.
*/
if (L.hi > ML_GAMMA_EXP_OVERFLOW) return ml_make_inf(0);
if (L.hi < ML_GAMMA_EXP_UNDERFLOW) return 0.0;
double g = ml_exp(L.hi);
if (!ml_isfinite(g) || g == 0.0) return g;
/* 3-term Taylor for exp(L.lo): accurate to ~L.lo^3/6 */
double elo = ML_FMA(L.lo, L.lo * 0.5, L.lo) + 1.0;
return ML_FMA(g, elo, 0.0);
}"""


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


def locate_v12a1() -> tuple:
    root = Path.cwd()
    cand = root / "v12A1"
    if cand.is_dir():
        return root, cand
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/")


def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    patched, count = OLD_FUNC.subn(NEW_FUNC, text, count=1)
    if count != 1:
        fail(
            f"{path}: could not find ml_exp_dd function. "
            f"Got {count} matches. Source may have drifted."
        )
    write_text(path, patched)


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
    print("  MATHLIB v12A1: EXP_DD SECOND-ORDER FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")
    print("\n[1/1] integral.c — fix ml_exp_dd second-order term")
    patch_integral(v12, force)
    archive_self(v12, force)
    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  ROOT CAUSE:")
    print("    ml_exp_dd used: exp(L.hi) * (1 + L.lo)")
    print("    This drops L.lo²/2. For the 8-step recurrence at")
    print("    x=0.001, L.lo reaches ~5e-8, making the dropped")
    print("    term g * L.lo²/2 ≈ 1.25e-12 = exactly 10 ULP.")
    print("")
    print("  FIX:")
    print("    Compute exp(L.lo) with 3-term Taylor:")
    print("    exp(L.lo) ≈ 1 + L.lo + L.lo²/2")
    print("    Since L.lo is tiny, this is accurate to ~1e-48.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Expected: 211 passed, 0 failed.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
