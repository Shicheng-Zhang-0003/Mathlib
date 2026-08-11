#!/usr/bin/env python3
"""
47_a1_remove_diag_gamma.py

Removes the gamma diagnostic tool that is no longer needed.

Background:
  diag_gamma_upstream.c was a forensic diagnostic built during
  the gamma(0.001) nightmare to trace where the ULP error was
  being injected. It compared ml_lgamma, ml_gamma_new, libc exp,
  and ml_exp at specific inputs to isolate the error source.

  The gamma issue is now resolved (211 passed, 0 failed).
  The diagnostic is dead weight.

  It is NOT compiled by CMake.
  It is NOT run by run_all_tests.py.
  It is NOT part of any edge suite.
  It is scaffolding from the gamma war.

Action:
  Delete v12A1/tests/diag_gamma_upstream.c

Usage:
  python3 47_a1_remove_diag_gamma.py
  python3 47_a1_remove_diag_gamma.py --force
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_DIAG_GAMMA_REMOVED"


def fail(msg: str) -> None:
    print("ERROR: " + msg)
    sys.exit(1)


def locate_v12a1() -> tuple:
    root = Path.cwd()
    cand = root / "v12A1"
    if cand.is_dir():
        return root, cand
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/.")


def remove_diagnostic(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "diag_gamma_upstream.c"

    if not path.is_file():
        print(f"  [skip] {path}: already removed")
        return

    path.unlink()
    print(f"  [delete] {path}")


def write_removal_record(v12: Path, force: bool) -> None:
    """Write a small record so the removal is documented."""
    record_path = v12 / "tests" / "regression" / "README.md"
    if not record_path.is_file():
        return

    text = record_path.read_text(encoding="utf-8")
    if MARKER in text and not force:
        print(f"  [skip] {record_path}: removal already recorded")
        return

    text += (
        "\n"
        "<!-- " + MARKER + " -->\n"
        "\n"
        "## Removed: diag_gamma_upstream.c\n"
        "\n"
        "The gamma diagnostic tool was removed during A1 closure.\n"
        "It was a forensic instrument used to trace the gamma(0.001)\n"
        "ULP error during scripts 30-45. The issue is resolved.\n"
        "The diagnostic is no longer needed.\n"
    )

    record_path.write_text(text, encoding="utf-8", newline="\n")
    print(f"  [record] {record_path}")


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
    print("  MATHLIB v12A1: REMOVE GAMMA DIAGNOSTIC")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] Removing diag_gamma_upstream.c")
    remove_diagnostic(v12, force)

    print("\n[2/2] Recording removal")
    write_removal_record(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Diagnostic removed.")
    print("")
    print("  What was removed:")
    print("    tests/diag_gamma_upstream.c")
    print("")
    print("  Why:")
    print("    It was a forensic tool for the gamma(0.001) ULP hunt.")
    print("    The issue is resolved. The tool is dead weight.")
    print("    It was never compiled by CMake or run by any test suite.")
    print("")
    print("  No build or test changes required.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
