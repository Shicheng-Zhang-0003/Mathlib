#!/usr/bin/env python3
"""
24_a1_1_1_doctrine_freeze_fix.py

Run from the folder that CONTAINS the v12A1 working folder.

Corrective follow-up to 23_a1_1_1_feature_freeze.py.

The previous freeze script marked TEAM_DOCTRINE but did not replace the old
"new algorithms are in scope" current-state block.

This script replaces the TEAM_DOCTRINE current-state block with the true
A1 closure-freeze state.

Targets:
    v12A1/docs/TEAM_DOCTRINE.md

Usage:
    python3 24_a1_1_1_doctrine_freeze_fix.py
    python3 24_a1_1_1_doctrine_freeze_fix.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_A1_FREEZE_DOCTRINE_FIX"

NEW_BLOCK = f"""Current state:
> v12A1 A1 closure freeze.
> v11S shipped. A1 feature development is frozen.

Therefore:
- no new modules, math families, or public APIs,
- only A1 closure fixes, tests, validation, docs, and process hygiene,
- correctness contracts still apply,
- script-only change policy still applies,
- zero allocation and thread safety still apply.

<!-- {MARKER} -->"""


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


def locate_v12a1() -> tuple[Path, Path]:
    root = Path.cwd()
    candidate = root / "v12A1"

    if candidate.is_dir():
        return root, candidate

    if (root / "docs" / "TEAM_DOCTRINE.md").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root

    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_team_doctrine(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "TEAM_DOCTRINE.md"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: doctrine freeze fix already present")
        return

    pattern = re.compile(
        r"(?ms)^[ \t]*Current state:.*?"
        r"^[ \t]*- zero allocation and thread safety still apply\.[ \t]*$"
    )

    patched, count = pattern.subn(lambda m: NEW_BLOCK, text, count=1)

    if count != 1:
        heading_pattern = re.compile(r"(?m)^#[ \t]+MathLib Team Doctrine[ \t]*$")
        patched, count = heading_pattern.subn(
            lambda m: m.group(0) + "\n\n" + NEW_BLOCK,
            text,
            count=1,
        )
        if count != 1:
            fail(f"{path}: could not find Current state block or doctrine heading.")

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
    print("  MATHLIB v12A1: A1 DOCTRINE FREEZE FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    patch_team_doctrine(v12, force)
    archive_self(v12, force)

    print("---------------------------------------------------------")
    print("  Doctrine freeze fix applied.")
    print("")
    print("  Verify:")
    print("")
    print('    grep -n "v12A1 A1 closure freeze" v12A1/docs/TEAM_DOCTRINE.md')
    print('    grep -n "new algorithms are in scope" v12A1/docs/TEAM_DOCTRINE.md')
    print("")
    print("  Expected:")
    print("    first grep: at least one hit")
    print("    second grep: no hits")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
