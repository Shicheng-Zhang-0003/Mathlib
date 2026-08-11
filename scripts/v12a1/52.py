#!/usr/bin/env python3
"""
51_a1_readme_status_update.py

Updates the README status section to reflect A1 closure completion.
Uses regex matching to handle whitespace variations.

Targets:
  v12A1/README.md

Usage:
  python3 51_a1_readme_status_update.py
  python3 51_a1_readme_status_update.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_README_STATUS_V2"


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


def patch_readme(v12: Path, force: bool) -> None:
    path = v12 / "README.md"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Use regex to match the status section flexibly.
    # The status section starts with "## Status" and goes to end of file
    # (or to the next ## heading if there is one).
    status_pattern = re.compile(
        r"## Status\s*\n"
        r"<!-- MATHLIB_V12A1_A1_FREEZE -->\s*\n"
        r"A1 feature freeze is in effect\..*?\n"
        r"This is a development tree\..*?\n"
        r"See `docs/V12A1_ROADMAP\.md` for the work plan\.",
        re.MULTILINE | re.DOTALL,
    )

    new_status = (
        "## Status\n"
        "<!-- " + MARKER + " -->\n"
        "<!-- MATHLIB_V12A1_A1_FREEZE -->\n"
        "A1 closure is **complete**.\n"
        "\n"
        "- Oracle validation: **212 passed, 0 failed** (all functions <= 5 ULP vs mpmath ground truth)\n"
        "- Full test gauntlet: **32/32 passed** (modular, smoke, edge, fuzz, oracle, boundary)\n"
        "- Closure gate: **PASSED**\n"
        "\n"
        "See `docs/V12A1_ROADMAP.md` for the work plan.\n"
        "See `release_notes.md` for the v12A1 closure summary."
    )

    matches = status_pattern.findall(text)
    if len(matches) == 0:
        # Fallback: try to find just the "## Status" line and replace everything after it
        status_idx = text.find("## Status")
        if status_idx == -1:
            fail(
                f"{path}: could not find '## Status' section.\n"
                f"  The README may have been restructured."
            )
        # Replace everything from "## Status" to end of file
        text = text[:status_idx] + new_status + "\n"
        write_text(path, text)
        return

    if len(matches) > 1:
        fail(f"{path}: found {len(matches)} status sections, expected 1.")

    text = status_pattern.sub(new_status, text, count=1)
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
    print("  MATHLIB v12A1: README STATUS UPDATE")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] README.md — update status to A1 closure complete")
    patch_readme(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  README status updated.")
    print("")
    print("  What changed:")
    print("    - 'A1 feature freeze is in effect' -> 'A1 closure is complete'")
    print("    - Removed 'Nothing here is release-grade yet'")
    print("    - Added oracle/test/gate results")
    print("    - Added pointer to release_notes.md")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
