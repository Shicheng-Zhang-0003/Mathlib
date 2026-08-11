#!/usr/bin/env python3
"""
48_a1_oracle_lgamma_0001.py

Adds lgamma(0.001) to the oracle validation set.

Background:
  gamma(0.001) was THE nightmare input that haunted v12A1 through
  scripts 30-45. The oracle has gamma(0.001) but is missing the
  corresponding lgamma(0.001) entry. Given that the 1-step
  recurrence fix (script 45) changed the lgamma path for x < 0.5,
  this input must be validated.

Action:
  1. Patch scripts/oracles/generate_oracles.py to add 0.001
     to lgamma_inputs.
  2. Run the generator to regenerate tests/oracle_data.h.
  3. Print verification instructions.

Targets:
  v12A1/scripts/oracles/generate_oracles.py
  v12A1/tests/oracle_data.h (regenerated)

Usage:
  python3 48_a1_oracle_lgamma_0001.py
  python3 48_a1_oracle_lgamma_0001.py --force
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_ORACLE_LGAMMA_0001"


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


def patch_generator(v12: Path, force: bool) -> None:
    path = v12 / "scripts" / "oracles" / "generate_oracles.py"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # The current lgamma_inputs block
    old_block = (
        "lgamma_inputs = [\n"
        "    0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 100.0, 171.0,\n"
        "    0.1, 0.01,\n"
        "    -0.5, -1.5, -2.5, -0.1,\n"
        "]"
    )

    new_block = (
        "# " + MARKER + "\n"
        "lgamma_inputs = [\n"
        "    0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 100.0, 171.0,\n"
        "    0.1, 0.01, 0.001,\n"
        "    -0.5, -1.5, -2.5, -0.1,\n"
        "]"
    )

    if old_block not in text:
        fail(
            f"{path}: could not find the lgamma_inputs block.\n"
            f"  Expected: 'lgamma_inputs = [\\n    0.5, 1.0, ...\\n    0.1, 0.01,\\n    ...]'\n"
            f"  Source may have drifted."
        )

    text = text.replace(old_block, new_block, 1)
    write_text(path, text)


def run_generator(v12: Path) -> None:
    gen_script = v12 / "scripts" / "oracles" / "generate_oracles.py"
    oracle_header = v12 / "tests" / "oracle_data.h"

    if not gen_script.is_file():
        fail(f"Missing generator: {gen_script}")

    print(f"  [run] python3 {gen_script}")
    result = subprocess.run(
        [sys.executable, str(gen_script)],
        cwd=str(v12),
        capture_output=True,
        text=True,
        timeout=120,
    )

    if result.returncode != 0:
        print(f"  [stdout] {result.stdout}")
        print(f"  [stderr] {result.stderr}")
        fail("Oracle generator failed. Is mpmath installed? (pip install mpmath)")

    print(f"  [stdout] {result.stdout.strip()}")

    if not oracle_header.is_file():
        fail(f"Generator did not produce: {oracle_header}")

    # Verify the new entry is present
    content = oracle_header.read_text(encoding="utf-8")
    if "oracle_lgamma_count = 17" not in content:
        fail(
            f"{oracle_header}: expected oracle_lgamma_count = 17 "
            f"(16 original + 1 new lgamma(0.001) entry). "
            f"Generator may not have run correctly."
        )

    print(f"  [ok] oracle_data.h regenerated with 17 lgamma entries")


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
    print("  MATHLIB v12A1: ADD lgamma(0.001) TO ORACLE")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] generate_oracles.py — add 0.001 to lgamma_inputs")
    patch_generator(v12, force)

    print("\n[2/2] Regenerating oracle_data.h")
    run_generator(v12)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Oracle expanded.")
    print("")
    print("  What changed:")
    print("    - lgamma_inputs now includes 0.001")
    print("    - oracle_lgamma_count: 16 -> 17")
    print("    - oracle_data.h regenerated with mpmath ground truth")
    print("")
    print("  Why:")
    print("    gamma(0.001) was the nightmare input. The lgamma path")
    print("    for x < 0.5 was rewritten (1-step recurrence, script 45).")
    print("    The oracle must validate lgamma at the same point.")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("")
    print("  Expected: 212 passed, 0 failed.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
