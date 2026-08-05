#!/usr/bin/env python3
"""
23_a1_1_1_feature_freeze.py

Run from the folder that CONTAINS the v12A1 working folder.

A1 closure table subsection 1.1:
    Feature freeze / closure lock.

Purpose:
    v12A1 must stop accepting new feature work before A2.
    This script makes the A1 freeze explicit and durable.

Targets:
    v12A1/docs/V12A1_A1_FREEZE.md
    v12A1/README.md
    v12A1/docs/V12A1_ROADMAP.md
    v12A1/docs/TEAM_DOCTRINE.md

Usage:
    python3 23_a1_1_1_feature_freeze.py
    python3 23_a1_1_1_feature_freeze.py --force
"""

from __future__ import annotations

import shutil
import sys
from datetime import date
from pathlib import Path

MARKER = "MATHLIB_V12A1_A1_FREEZE"
TODAY = date.today().isoformat()

NEW_FREEZE_DOC = f"""# v12A1 A1 Feature Freeze

<!-- {MARKER} -->

Subsection: 1.1
Effective: {TODAY}

v12A1 is now in **A1 closure freeze**.

## Allowed during A1 closure

Only the following are permitted:

1. Correctness fixes required by the A1 closure table.
2. Test and oracle expansion needed to prove those fixes.
3. Validation and gate work.
4. Documentation alignment.
5. Script/process hygiene directly tied to A1 closure.

## Not allowed during A1 closure

The following are frozen:

1. New modules.
2. New public APIs.
3. New math families.
4. New performance experiments.
5. Speculative refactors.
6. Feature creep of any kind.
7. Anything whose primary purpose is to make A2 easier.

## Script rule

Every A1 closure change must be applied through a numbered script in:

```text
v12A1/scripts/v12a1/
```

and that script must correspond to a specific A1 closure subsection.

## Closure mindset

The goal is no longer expansion.

The goal is:

> make the current tree true, tested, documented, and stable enough to close A1.
"""

README_FREEZE_BLOCK = f"""
<!-- {MARKER} -->
A1 feature freeze is in effect. Only A1 closure fixes, tests, validation, and documentation alignment are allowed.
"""

ROADMAP_FREEZE_SECTION = f"""
<!-- {MARKER} -->
## A1 Closure Freeze (Subsection 1.1)

- Effective: {TODAY}
- No new modules.
- No new public APIs.
- No new math families.
- No speculative features.
- Only A1 closure table fixes, tests, oracle expansion, validation, docs alignment, and script/process hygiene are allowed.
- Each change must be applied by a numbered script corresponding to an A1 subsection.

"""

TEAM_OLD_STATE = """Current state:
> v12A1 early A1 development wave.
> v11S shipped. Closure documents archived.
Therefore:
- new algorithms are in scope,
- approximation redesigns are the primary mission,
- correctness contracts still apply,
- script-only change policy still applies,
- zero allocation and thread safety still apply."""

TEAM_NEW_STATE = f"""Current state:
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

    # Convenience: allow running from inside v12A1 itself.
    if (
        (root / "src" / "core.c").is_file()
        and (root / "include" / "mathlib").is_dir()
    ):
        print("  [note] Running from inside v12A1.")
        return root.parent, root

    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def write_freeze_doc(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "V12A1_A1_FREEZE.md"

    if path.exists() and not force:
        try:
            old = normalize(path.read_text(encoding="utf-8"))
            if MARKER in old:
                print(f"  [skip] {path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(path, NEW_FREEZE_DOC)


def patch_readme(v12: Path, force: bool) -> None:
    path = v12 / "README.md"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: freeze notice already present")
        return

    anchor = "## Status"
    if anchor in text:
        text = text.replace(anchor, anchor + README_FREEZE_BLOCK, 1)
    else:
        text = text.rstrip() + "\n\n## Status\n" + README_FREEZE_BLOCK

    write_text(path, text)


def patch_roadmap(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "V12A1_ROADMAP.md"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: freeze section already present")
        return

    title = "# v12A1 Development Roadmap"
    if title in text:
        text = text.replace(title, title + ROADMAP_FREEZE_SECTION, 1)
    else:
        text = "# v12A1 Development Roadmap\n" + ROADMAP_FREEZE_SECTION + text

    write_text(path, text)


def patch_team_doctrine(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "TEAM_DOCTRINE.md"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: freeze state already present")
        return

    if TEAM_OLD_STATE in text:
        text = text.replace(TEAM_OLD_STATE, TEAM_NEW_STATE, 1)
    else:
        text = text.rstrip() + f"""

<!-- {MARKER} -->
## A1 Closure Freeze

A1 feature freeze is in effect.

Only A1 closure fixes, tests, validation, documentation alignment, and script/process hygiene are allowed.
"""

    write_text(path, text)


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
    print("  MATHLIB v12A1: A1 SUBSECTION 1.1")
    print("  FEATURE FREEZE / CLOSURE LOCK")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/4] docs/V12A1_A1_FREEZE.md")
    write_freeze_doc(v12, force)

    print("\n[2/4] README.md")
    patch_readme(v12, force)

    print("\n[3/4] docs/V12A1_ROADMAP.md")
    patch_roadmap(v12, force)

    print("\n[4/4] docs/TEAM_DOCTRINE.md")
    patch_team_doctrine(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Subsection 1.1 applied.")
    print("")
    print("  What changed:")
    print("    - A1 freeze policy written to docs/V12A1_A1_FREEZE.md")
    print("    - README status now states the freeze")
    print("    - Roadmap now contains the A1 freeze section")
    print("    - Team Doctrine switched from A1 feature wave to A1 closure freeze")
    print("")
    print("  Verify:")
    print("    grep -R \"MATHLIB_V12A1_A1_FREEZE\" v12A1")
    print("")
    print("  Next: 1.2")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
