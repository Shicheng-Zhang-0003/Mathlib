#!/usr/bin/env python3
"""
49_a1_integral_comment_cleanup.py

Removes stale comments from integral.c left behind by scripts 44/45.

After the 1-step recurrence fix (script 45), the following comment
blocks are now stale and misleading:

  1. ml_lgamma x<0.5:
     MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS comments
     (replaced by MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE)

  2. ml_gamma_new x<0.5:
     MATHLIB_V12A1_GAMMA_NEW_LOGSPACE comments
     MATHLIB_V12A1_GAMMA_NEW_DEPTH16 comments
     (replaced by MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE)

  3. Top-level gamma section comment:
     "x < 0.5 uses 8-step recurrence (not 1-step)"
     should now say "x < 0.5 uses 1-step recurrence"

Targets:
  v12A1/src/integral.c

Usage:
  python3 49_a1_integral_comment_cleanup.py
  python3 49_a1_integral_comment_cleanup.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_INTEGRAL_COMMENT_CLEANUP"


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

    changes = 0

    # ----------------------------------------------------------------
    # Fix 1: Remove stale DIRECT_LANCZOS comments from ml_lgamma
    # ----------------------------------------------------------------
    # These comments were added by script 44 and are now superseded
    # by the 1STEP_RECURRENCE comments from script 45.
    direct_lanczos_lgamma = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS \*/\s*\n"
        r"/\* x < 0\.5: use Lanczos directly via ml_lgamma_positive_dd\.\s*\n"
        r" \* The old recurrence formula subtracted two large DD values\s*\n"
        r" \* \(lgamma\(x\+16\) ~ 28 minus sum-of-logs ~ 26\) to get a small\s*\n"
        r" \* result \(~2\.25\), amplifying rounding errors by ~12\.5x\.\s*\n"
        r" \* Lanczos computes lgamma\(x\) directly with ~3\.9x cancellation\.\s*\n"
        r" \*/\s*\n",
        re.MULTILINE,
    )

    matches = direct_lanczos_lgamma.findall(text)
    if matches:
        text = direct_lanczos_lgamma.sub("", text)
        changes += 1
        print("  [ok] Removed stale DIRECT_LANCZOS comments from ml_lgamma")
    else:
        print("  [skip] DIRECT_LANCZOS comments not found in ml_lgamma (already clean)")

    # ----------------------------------------------------------------
    # Fix 2: Remove stale LOGSPACE + DEPTH16 comments from ml_gamma_new
    # ----------------------------------------------------------------
    # These comments were added by scripts 37 and 40 and are now
    # superseded by the 1STEP_RECURRENCE comments from script 45.

    # LOGSPACE comment block (multi-line)
    logspace_block = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_NEW_LOGSPACE \*/\s*\n"
        r"/\* x < 0\.5: 8-step recurrence in log-space.*?\*/\s*\n",
        re.MULTILINE | re.DOTALL,
    )

    matches = logspace_block.findall(text)
    if matches:
        text = logspace_block.sub("", text)
        changes += 1
        print("  [ok] Removed stale LOGSPACE comments from ml_gamma_new")
    else:
        print("  [skip] LOGSPACE comments not found in ml_gamma_new (already clean)")

    # DEPTH16 comment block
    depth16_block = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_NEW_DEPTH16 \*/\s*\n"
        r"/\* x < 0\.5: 16-step recurrence in log-space\.\s*\n"
        r"\* Depth 16 pushes Stirling to x\+16 where the correction\s*\n"
        r"\* is ~0\.0052 instead of ~0\.0104, halving the DD chain error\.\s*\n"
        r"\* This matches ml_lgamma's depth-16 path \(script 39/40\)\.\s*\n"
        r"\*/\s*\n",
        re.MULTILINE,
    )

    matches = depth16_block.findall(text)
    if matches:
        text = depth16_block.sub("", text)
        changes += 1
        print("  [ok] Removed stale DEPTH16 comments from ml_gamma_new")
    else:
        print("  [skip] DEPTH16 comments not found in ml_gamma_new (already clean)")

    # ----------------------------------------------------------------
    # Fix 3: Update top-level gamma section comment
    # ----------------------------------------------------------------
    old_line = "*   - x < 0.5 uses 8-step recurrence (not 1-step)."
    new_line = "*   - x < 0.5 uses 1-step recurrence: lgamma(x) = lgamma(x+1) - log(x)."

    if old_line in text:
        text = text.replace(old_line, new_line, 1)
        changes += 1
        print("  [ok] Updated top-level gamma section comment")
    else:
        print("  [skip] Top-level gamma comment already updated or not found")

    # ----------------------------------------------------------------
    # Add cleanup marker
    # ----------------------------------------------------------------
    # Add the marker to the top-level gamma section comment
    marker_line = "/* " + MARKER + " */\n"
    halfint_marker = "/* MATHLIB_V12A1_GAMMA_HALFINT_V7 */"
    if halfint_marker in text and marker_line not in text:
        text = text.replace(halfint_marker, marker_line + halfint_marker, 1)
        changes += 1

    if changes == 0:
        print(f"  [skip] {path}: no stale comments found")
        return

    write_text(path, text)
    print(f"  [info] Applied {changes} cleanup change(s)")


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
    print("  MATHLIB v12A1: INTEGRAL COMMENT CLEANUP")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — remove stale comments from scripts 44/45")
    patch_integral(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Comment cleanup applied.")
    print("")
    print("  What was removed:")
    print("    - MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS comments (ml_lgamma)")
    print("    - MATHLIB_V12A1_GAMMA_NEW_LOGSPACE comments (ml_gamma_new)")
    print("    - MATHLIB_V12A1_GAMMA_NEW_DEPTH16 comments (ml_gamma_new)")
    print("")
    print("  What was updated:")
    print("    - Top-level gamma section: '8-step recurrence' -> '1-step recurrence'")
    print("")
    print("  No build or test changes required.")
    print("  This is a documentation-only cleanup.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
