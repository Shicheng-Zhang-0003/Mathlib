#!/usr/bin/env python3
"""
33_a1_gamma_divide_fix.py
Run from the folder that CONTAINS the v12A1 working folder.

Fixes the LAST oracle failure: gamma(0.001) at 11 ULP.

Root cause:
  ml_gamma_positive() for x < 8 uses product-then-divide:
    gy = ml_exp_dd(Ly)    <- rounds to double
    p  = P.hi + P.lo      <- rounds to double
    return gy / p          <- division rounds AGAIN

  Two extra roundings after the DD chain. The lgamma path
  (which passes) uses log-subtract-then-exp, keeping everything
  in DD until the final ml_exp_dd.

Fix:
  Replace the product-then-divide path with:
    L = ml_lgamma_positive_dd(x)   <- DD log subtraction
    return ml_exp_dd(L)            <- single DD->double at the end

Targets:
  v12A1/src/integral.c

Usage:
  python3 33_a1_gamma_divide_fix.py
  python3 33_a1_gamma_divide_fix.py --force
"""
from __future__ import annotations
import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_DIVIDE_FIX"

OLD_FUNC = re.compile(
    r"static double ml_gamma_positive\(double x\) \{.*?\n\}",
    re.DOTALL,
)

NEW_FUNC = """static double ml_gamma_positive(double x) {
/* """ + MARKER + """ */
/*
* For x >= 8: Stirling DD -> exp.
* For half-integers: exact product formula.
* For 0 < x < 8: use the DD lgamma path (log subtraction),
*   then exponentiate. This avoids the product-then-divide
*   approach which introduced two extra double roundings.
*/
if (x >= 8.0) {
ml_dd_t L = ml_stirling_lgamma_dd(x);
return ml_exp_dd(L);
}
if (ml_is_half_integer(x)) {
return ml_gamma_half_positive(x);
}
/* DD log-subtract-then-exp: same path as lgamma, proven accurate */
ml_dd_t L = ml_lgamma_positive_dd(x);
return ml_exp_dd(L);
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
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    patched, count = OLD_FUNC.subn(NEW_FUNC, text, count=1)
    if count != 1:
        fail(
            f"{path}: could not find ml_gamma_positive function. "
            f"Got {count} matches. Source may have drifted."
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
    print("  MATHLIB v12A1: GAMMA DIVIDE FIX (final ULP)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")
    print("\n[1/1] integral.c — replace product-divide with DD lgamma path")
    patch_integral(v12, force)
    archive_self(v12, force)
    print("\n---------------------------------------------------------")
    print("  Final gamma fix applied.")
    print("")
    print("  What changed:")
    print("    OLD: gy = ml_exp_dd(Ly); p = P.hi+P.lo; return gy/p;")
    print("         (two extra double roundings)")
    print("    NEW: L = ml_lgamma_positive_dd(x); return ml_exp_dd(L);")
    print("         (everything stays in DD until final exp)")
    print("")
    print("  Why this fixes gamma(0.001):")
    print("    The lgamma path for x=0.001 PASSES the oracle.")
    print("    It uses DD log subtraction: lgamma(8.001) - sum(log(x+k)).")
    print("    The gamma path now uses the SAME DD chain, then exponentiates.")
    print("    No intermediate double rounding. No division rounding.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Expected: 211 passed, 0 failed. ALL EDGE TESTS PASSED.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
