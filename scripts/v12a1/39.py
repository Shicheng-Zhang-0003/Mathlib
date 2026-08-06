#!/usr/bin/env python3
"""
39_a1_gamma_recurrence_depth16.py

Fixes the LAST oracle failure: gamma(0.001) at 10 ULP.

Root cause:
The 8-step recurrence for x < 0.5 evaluates Stirling at x+8 = 8.001.
At x = 8.001, the Stirling correction (~0.0104) is still large enough
that the DD chain accumulates ~1 ULP of error in the log domain,
which exp() amplifies to ~10 ULP in the gamma domain.

Fix:
Increase recurrence depth from 8 to 16 for x < 0.5.
Stirling is now evaluated at x+16 = 16.001, where the correction
is ~0.0052 (half as large), reducing the DD chain error.

Also fixes the diagnostic's wrong oracle_lgamma value.

Targets:
  v12A1/src/integral.c
  v12A1/tests/diag_gamma_upstream.c
"""

from __future__ import annotations
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH16"

# --- Fix 1: integral.c recurrence depth 8 -> 16 ---

OLD_RECURRENCE = """\
        if (x < 0.5) {
            ml_dd_t L = ml_lgamma_positive_dd(x + 8.0);
            for (int k = 0; k < 8; k++)
                L = ml_dd_sub(L, ml_log_dd(x + (double)k));
            return ml_exp_dd(L);
        }"""

NEW_RECURRENCE = """\
        /* """ + MARKER + """ */
        /* x < 0.5: 16-step recurrence in log-space.
         * Depth 16 pushes Stirling to x+16 where the correction
         * is smaller, reducing DD chain error from ~1 ULP to <0.5 ULP.
         */
        if (x < 0.5) {
            ml_dd_t L = ml_lgamma_positive_dd(x + 16.0);
            for (int k = 0; k < 16; k++)
                L = ml_dd_sub(L, ml_log_dd(x + (double)k));
            return ml_exp_dd(L);
        }"""

# Same fix for ml_lgamma's x < 0.5 path
OLD_LGAMMA_REC = """\
        if (x < 0.5) {
            ml_dd_t L = ml_lgamma_positive_dd(x + 8.0);
            for (int k = 0; k < 8; k++)
                L = ml_dd_sub(L, ml_log_dd(x + (double)k));
            return L.hi + L.lo;
        }"""

NEW_LGAMMA_REC = """\
        /* """ + MARKER + """ */
        if (x < 0.5) {
            ml_dd_t L = ml_lgamma_positive_dd(x + 16.0);
            for (int k = 0; k < 16; k++)
                L = ml_dd_sub(L, ml_log_dd(x + (double)k));
            return L.hi + L.lo;
        }"""

# --- Fix 2: diagnostic oracle_lgamma ---
OLD_ORACLE = 'double oracle_lgamma = 6.90775527898213682e+00;'
NEW_ORACLE = 'double oracle_lgamma = 6.90717888538385250e+00; /* ln(gamma(0.001)), NOT ln(1000) */'


def fail(msg):
    print("ERROR: " + msg)
    sys.exit(1)

def normalize(t):
    return t.replace("\r\n", "\n").replace("\r", "\n")

def write_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")

def locate_v12a1():
    root = Path.cwd()
    cand = root / "v12A1"
    if cand.is_dir():
        return root, cand
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/")

def patch_integral(v12, force):
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    count = 0
    if OLD_RECURRENCE in text:
        text = text.replace(OLD_RECURRENCE, NEW_RECURRENCE, 1)
        count += 1
    if OLD_LGAMMA_REC in text:
        text = text.replace(OLD_LGAMMA_REC, NEW_LGAMMA_REC, 1)
        count += 1
    if count == 0:
        fail(f"{path}: could not find recurrence blocks. Source may have drifted.")
    write_text(path, text)
    print(f"  [info] Replaced {count} recurrence block(s)")

def patch_diagnostic(v12, force):
    path = v12 / "tests" / "diag_gamma_upstream.c"
    if not path.is_file():
        print(f"  [skip] {path}: not found")
        return
    text = normalize(path.read_text(encoding="utf-8"))
    if OLD_ORACLE in text:
        text = text.replace(OLD_ORACLE, NEW_ORACLE, 1)
        write_text(path, text)
    else:
        print(f"  [skip] {path}: oracle already fixed or drifted")

def archive_self(v12, force):
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

def main():
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()
    print("=========================================================")
    print("  MATHLIB v12A1: GAMMA RECURRENCE DEPTH 16")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")
    print("\n[1/2] integral.c — recurrence depth 8 -> 16")
    patch_integral(v12, force)
    print("\n[2/2] diag_gamma_upstream.c — fix wrong oracle_lgamma")
    patch_diagnostic(v12, force)
    archive_self(v12, force)
    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  What changed:")
    print("    - x < 0.5 recurrence depth: 8 -> 16")
    print("    - Stirling evaluated at x+16 instead of x+8")
    print("    - Diagnostic oracle_lgamma corrected")
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
