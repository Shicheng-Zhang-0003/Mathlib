#!/usr/bin/env python3
"""
37_a1_gamma_new_logspace.py
Run from the folder that CONTAINS the v12A1 working folder.

Fixes the LAST oracle failure: gamma(0.001) at 11 ULP.

Root cause:
  ml_gamma_new() for x < 0.5 still uses product-then-divide:
    gy = ml_exp_dd(Ly)     <- rounds to double
    p  = P.hi + P.lo       <- rounds to double
    return gy / p           <- division rounds AGAIN

  Script 33 fixed ml_gamma_positive but NOT ml_gamma_new.
  The fix: same log-space approach used in ml_lgamma (which passes).

Targets:
  v12A1/src/integral.c

Usage:
  python3 37_a1_gamma_new_logspace.py
  python3 37_a1_gamma_new_logspace.py --force
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_NEW_LOGSPACE"

OLD_BLOCK = """\
/* x < 0.5: 8-step recurrence via product */
if (x < 0.5) {
ml_dd_t Ly = ml_lgamma_positive_dd(x + 8.0);
double gy = ml_exp_dd(Ly);
ml_dd_t P = ml_dd_from_d(1.0);
for (int k = 0; k < 8; k++)
P = ml_dd_mul_d(P, x + (double)k);
double p = P.hi + P.lo;
if (p == 0.0) return ml_make_inf(0);
return gy / p;
}"""

NEW_BLOCK = """\
/* """ + MARKER + """ */
/* x < 0.5: 8-step recurrence in log-space (no product-divide)
*
* Old: gy = exp(lgamma(x+8)); p = prod(x+k); return gy/p;
*   -> two extra double roundings + division rounding = 11 ULP
*
* New: lgamma(x+8) - sum(log(x+k)) then exp
*   -> everything stays in DD until final exp, same as ml_lgamma
*/
if (x < 0.5) {
ml_dd_t L = ml_lgamma_positive_dd(x + 8.0);
for (int k = 0; k < 8; k++)
L = ml_dd_sub(L, ml_log_dd(x + (double)k));
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
    if OLD_BLOCK not in text:
        # Try with different whitespace
        alt = text.replace("    ", "")
        alt_old = OLD_BLOCK.replace("    ", "")
        if alt_old not in alt:
            fail(
                f"{path}: could not find the x<0.5 product-divide block.\n"
                f"  Searched for:\n"
                f"    /* x < 0.5: 8-step recurrence via product */\n"
                f"  Current file may have drifted."
            )
        # Whitespace-insensitive replacement
        import re
        pattern = re.compile(
            r"/\* x < 0\.5: 8-step recurrence via product \*/\s*"
            r"if \(x < 0\.5\) \{\s*"
            r"ml_dd_t Ly = ml_lgamma_positive_dd\(x \+ 8\.0\);\s*"
            r"double gy = ml_exp_dd\(Ly\);\s*"
            r"ml_dd_t P = ml_dd_from_d\(1\.0\);\s*"
            r"for \(int k = 0; k < 8; k\+\+\)\s*"
            r"P = ml_dd_mul_d\(P, x \+ \(double\)k\);\s*"
            r"double p = P\.hi \+ P\.lo;\s*"
            r"if \(p == 0\.0\) return ml_make_inf\(0\);\s*"
            r"return gy / p;\s*"
            r"\}",
            re.MULTILINE,
        )
        patched, count = pattern.subn(NEW_BLOCK, text, count=1)
        if count != 1:
            fail(f"{path}: regex replacement failed ({count} matches).")
        write_text(path, patched)
        return
    text = text.replace(OLD_BLOCK, NEW_BLOCK, 1)
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
    print("  MATHLIB v12A1: GAMMA_NEW LOG-SPACE FIX (final)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")
    print("\n[1/1] integral.c — replace product-divide in ml_gamma_new")
    patch_integral(v12, force)
    archive_self(v12, force)
    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  What changed:")
    print("    OLD: Ly=lgamma(x+8); gy=exp(Ly); P=prod(x+k); return gy/P")
    print("         (two double roundings + division rounding)")
    print("    NEW: L=lgamma(x+8)-sum(log(x+k)); return exp(L)")
    print("         (everything stays in DD until final exp)")
    print("")
    print("  This matches the ml_lgamma path for x<0.5 which PASSES.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Expected: 211 passed, 0 failed. ALL TESTS PASSED.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
