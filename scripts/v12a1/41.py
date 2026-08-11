#!/usr/bin/env python3
"""
41_a1_gamma_new_depth16.py

Fixes the LAST oracle failure: gamma(0.001) at ~10 ULP.

Root cause:
  Script 40 was supposed to patch BOTH ml_gamma_new and ml_lgamma
  to depth 16. But script 37 had already replaced the ml_gamma_new
  x<0.5 block and added a MATHLIB_V12A1_GAMMA_NEW_LOGSPACE comment
  header above it. Script 40's string match expected the bare block
  without that comment, so it silently skipped ml_gamma_new.

  Result:
    ml_lgamma  x < 0.5 -> depth 16  (script 39/40 applied)
    ml_gamma_new x < 0.5 -> depth 8  (script 40 MISSED)

  At x = 0.001, depth 8 evaluates Stirling at 8.001 where the
  correction is ~0.0104. The DD chain accumulates ~1 ULP in the
  log domain, which exp() amplifies to ~10 ULP.

Fix:
  Change ml_gamma_new's x < 0.5 recurrence from depth 8 to depth 16,
  matching ml_lgamma. Stirling is now evaluated at x+16 = 16.001
  where the correction is ~0.0052, halving the injected error.

Targets:
  v12A1/src/integral.c

Usage:
  python3 41_a1_gamma_new_depth16.py
  python3 41_a1_gamma_new_depth16.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_NEW_DEPTH16"

# The current ml_gamma_new x<0.5 block, as it exists after scripts 37+38.
# We match the if-block itself (not the comment above it) to avoid
# the same string-drift failure that killed script 40.
OLD_GAMMA_NEW = re.compile(
    r"(if \(x < 0\.5\) \{\s*\n"
    r"\s*ml_dd_t L = ml_lgamma_positive_dd\(x \+ 8\.0\);\s*\n"
    r"\s*for \(int k = 0; k < 8; k\+\+\)\s*\n"
    r"\s*L = ml_dd_sub\(L, ml_log_dd\(x \+ \(double\)k\)\);\s*\n"
    r"\s*return ml_exp_dd\(L\);\s*\n"
    r"\s*\})",
    re.MULTILINE,
)

NEW_GAMMA_NEW = """\
/* """ + MARKER + """ */
/* x < 0.5: 16-step recurrence in log-space.
* Depth 16 pushes Stirling to x+16 where the correction
* is ~0.0052 instead of ~0.0104, halving the DD chain error.
* This matches ml_lgamma's depth-16 path (script 39/40).
*/
if (x < 0.5) {
    ml_dd_t L = ml_lgamma_positive_dd(x + 16.0);
    for (int k = 0; k < 16; k++)
        L = ml_dd_sub(L, ml_log_dd(x + (double)k));
    return ml_exp_dd(L);
}"""


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

    # Verify the lgamma path is already at depth 16 (sanity check).
    if "MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH16" not in text:
        fail(
            f"{path}: ml_lgamma depth-16 marker not found. "
            f"Run scripts 39/40 first."
        )

    # Count how many depth-8 gamma blocks exist.
    # There should be exactly ONE: the ml_gamma_new x<0.5 block.
    # (ml_lgamma's x<0.5 was already patched to depth 16.)
    matches = OLD_GAMMA_NEW.findall(text)
    if len(matches) == 0:
        fail(
            f"{path}: could not find the depth-8 gamma_new block.\n"
            f"  The source may have drifted. Expected:\n"
            f"    if (x < 0.5) {{\n"
            f"        ml_dd_t L = ml_lgamma_positive_dd(x + 8.0);\n"
            f"        for (int k = 0; k < 8; k++)\n"
            f"            ...\n"
            f"        return ml_exp_dd(L);\n"
            f"    }}"
        )
    if len(matches) > 1:
        fail(
            f"{path}: found {len(matches)} depth-8 blocks, expected 1. "
            f"Refusing to patch ambiguous source."
        )

    patched, count = OLD_GAMMA_NEW.subn(NEW_GAMMA_NEW, text, count=1)
    if count != 1:
        fail(f"{path}: regex replacement failed ({count} matches).")

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
    print("  MATHLIB v12A1: GAMMA_NEW DEPTH 16 (script 40 follow-up)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — gamma_new recurrence depth 8 -> 16")
    patch_integral(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Fix applied.")
    print("")
    print("  What changed:")
    print("    ml_gamma_new x < 0.5 recurrence depth: 8 -> 16")
    print("    Stirling evaluated at x+16 instead of x+8")
    print("    Now matches ml_lgamma's depth-16 path")
    print("")
    print("  Why script 40 missed this:")
    print("    Script 37 added a comment header above the block.")
    print("    Script 40's string match expected the bare block.")
    print("    The match failed silently. Only ml_lgamma was patched.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: gamma(0.001) drops from ~10 ULP to <= 5 ULP.")
    print("  Expected: 211 passed, 0 failed.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
