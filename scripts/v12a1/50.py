#!/usr/bin/env python3
"""
49_a1_integral_comment_cleanup.py

Removes stale comments from integral.c left behind by scripts 39-45.

After the 1-step recurrence fix (script 45), the following comment
blocks are now stale and misleading:

  1. ml_lgamma x<0.5:
     "/* x < 0.5: 8-step recurrence */"
     "/* MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH16 */"
     (replaced by MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE)

  2. ml_gamma_new x<0.5:
     "/* MATHLIB_V12A1_GAMMA_NEW_LOGSPACE */" + multi-line comment
     "/* MATHLIB_V12A1_GAMMA_NEW_DEPTH16 */" + multi-line comment
     (replaced by MATHLIB_V12A1_GAMMA_1STEP_RECURRENCE)

  3. Any MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS comments from script 44

  4. Top-level gamma section comment:
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
    # Fix 1: Remove stale "8-step recurrence" comment from ml_lgamma
    # ----------------------------------------------------------------
    old_lgamma_comment = "/* x < 0.5: 8-step recurrence */\n"
    if old_lgamma_comment in text:
        text = text.replace(old_lgamma_comment, "", 1)
        changes += 1
        print("  [ok] Removed stale '8-step recurrence' comment from ml_lgamma")

    # ----------------------------------------------------------------
    # Fix 2: Remove stale RECURRENCE_DEPTH16 comment from ml_lgamma
    # ----------------------------------------------------------------
    depth16_lgamma = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH16 \*/\s*\n",
        re.MULTILINE,
    )
    matches = depth16_lgamma.findall(text)
    if matches:
        text = depth16_lgamma.sub("", text)
        changes += 1
        print("  [ok] Removed stale RECURRENCE_DEPTH16 comment from ml_lgamma")

    # ----------------------------------------------------------------
    # Fix 3: Remove stale LOGSPACE comment block from ml_gamma_new
    # ----------------------------------------------------------------
    logspace_block = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_NEW_LOGSPACE \*/\s*\n"
        r"/\* x < 0\.5: 8-step recurrence in log-space.*?\*/\s*\n",
        re.MULTILINE | re.DOTALL,
    )
    matches = logspace_block.findall(text)
    if matches:
        text = logspace_block.sub("", text)
        changes += 1
        print("  [ok] Removed stale LOGSPACE comment block from ml_gamma_new")

    # ----------------------------------------------------------------
    # Fix 4: Remove stale DEPTH16 comment block from ml_gamma_new
    # ----------------------------------------------------------------
    depth16_gamma = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_NEW_DEPTH16 \*/\s*\n"
        r"/\* x < 0\.5: 16-step recurrence in log-space\..*?\*/\s*\n",
        re.MULTILINE | re.DOTALL,
    )
    matches = depth16_gamma.findall(text)
    if matches:
        text = depth16_gamma.sub("", text)
        changes += 1
        print("  [ok] Removed stale DEPTH16 comment block from ml_gamma_new")

    # ----------------------------------------------------------------
    # Fix 5: Remove stale DIRECT_LANCZOS comments from script 44
    # ----------------------------------------------------------------
    direct_lanczos = re.compile(
        r"/\* MATHLIB_V12A1_GAMMA_DIRECT_LANCZOS \*/\s*\n"
        r"/\* x < 0\.5: use Lanczos directly.*?\*/\s*\n",
        re.MULTILINE | re.DOTALL,
    )
    matches = direct_lanczos.findall(text)
    if matches:
        text = direct_lanczos.sub("", text)
        changes += 1
        print("  [ok] Removed stale DIRECT_LANCZOS comments from script 44")

    # ----------------------------------------------------------------
    # Fix 6: Update top-level gamma section comment
    # ----------------------------------------------------------------
    old_line = "*   - x < 0.5 uses 8-step recurrence (not 1-step)."
    new_line = "*   - x < 0.5 uses 1-step recurrence: lgamma(x) = lgamma(x+1) - log(x)."
    if old_line in text:
        text = text.replace(old_line, new_line, 1)
        changes += 1
        print("  [ok] Updated top-level gamma section comment")

    # ----------------------------------------------------------------
    # Add cleanup marker
    # ----------------------------------------------------------------
    halfint_marker = "/* MATHLIB_V12A1_GAMMA_HALFINT_V7 */"
    if halfint_marker in text and MARKER not in text:
        text = text.replace(
            halfint_marker,
            "/* " + MARKER + " */\n" + halfint_marker,
            1,
        )
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

    print("\n[1/1] integral.c — remove stale comments from scripts 39-45")
    patch_integral(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Comment cleanup applied.")
    print("")
    print("  What was removed:")
    print("    - '8-step recurrence' comment (ml_lgamma)")
    print("    - RECURRENCE_DEPTH16 marker (ml_lgamma)")
    print("    - LOGSPACE comment block (ml_gamma_new)")
    print("    - DEPTH16 comment block (ml_gamma_new)")
    print("    - DIRECT_LANCZOS comments (script 44 leftovers)")
    print("")
    print("  What was updated:")
    print("    - Top-level gamma section: '8-step recurrence'")
    print("      -> '1-step recurrence: lgamma(x) = lgamma(x+1) - log(x)'")
    print("")
    print("  No build or test changes required.")
    print("  This is a documentation-only cleanup.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
