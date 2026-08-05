#!/usr/bin/env python3
"""
23_a1_payne_hanek_exp_constants.py

Run from the folder that CONTAINS the v12A1 working folder.

A1 oracle triage fixes:

1. Payne-Hanek large-argument reduction:
   The current k_start/k_end heuristic drops table terms that still
   contribute fractional bits. Replace it with a full-table sweep.
   The existing process_term() logic skips irrelevant terms.

2. ml_exp() Cody-Waite constants:
   Use ML_LN2_HI / ML_LN2_LO instead of hardcoded literals.
   The old exp path used a different low constant than the log path.

Targets:
    v12A1/src/internal/payne_hanek.h
    v12A1/src/exp_log.c

Usage:
    python3 23_a1_payne_hanek_exp_constants.py
    python3 23_a1_payne_hanek_exp_constants.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER_PH = "MATHLIB_V12A1_PAYNE_HANEK_V4_FULL_TABLE"
MARKER_EXP = "MATHLIB_V12A1_EXP_LN2_SPLIT_MACROS"


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

    if (
        (root / "src" / "internal" / "payne_hanek.h").is_file()
        and (root / "src" / "exp_log.c").is_file()
    ):
        print("  [note] Running from inside v12A1.")
        return root.parent, root

    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_payne_hanek(v12: Path, force: bool) -> None:
    path = v12 / "src" / "internal" / "payne_hanek.h"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER_PH in text and not force:
        print(f"  [skip] {path}: full-table marker already present")
        return

    # 1. Replace the restrictive k bounds with a full-table sweep.
    range_pattern = re.compile(
        r"int k_start = \(E - 77\) / 24;\s*\n"
        r"\s*if \(k_start < 0\) k_start = 0;\s*\n"
        r"\s*int k_end = \(E \+ 53\) / 24;\s*\n"
        r"\s*if \(k_end > 65\) k_end = 65;"
    )

    new_range = f"""/* {MARKER_PH}
     *
     * Full-table sweep.
     *
     * The old heuristic k bounds dropped table terms that could still
     * contribute fractional bits for some large arguments. The skip
     * logic inside ml_ph_process_term() already ignores irrelevant
     * integer-multiple-of-4 terms, so a full sweep is safe.
     */
    int k_start = 0;
    int k_end = 65;"""

    patched, count = range_pattern.subn(lambda m: new_range, text, count=1)
    if count != 1:
        fail(f"{path}: could not find Payne-Hanek k_start/k_end block.")

    # 2. Improve reduced-argument reconstruction with FMA.
    old_reconstruct = "double result = frac * ML_PH_PI2_HI + frac * ML_PH_PI2_LO;"
    new_reconstruct = (
        f"double result = ML_FMA(frac, ML_PH_PI2_HI, frac * ML_PH_PI2_LO); "
        f"/* {MARKER_PH} */"
    )

    if old_reconstruct in patched:
        patched = patched.replace(old_reconstruct, new_reconstruct, 1)
    else:
        print(f"  [warn] {path}: could not find result reconstruction line.")

    write_text(path, patched)


def patch_exp_constants(v12: Path, force: bool) -> None:
    path = v12 / "src" / "exp_log.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER_EXP in text and not force:
        print(f"  [skip] {path}: exp ln2 split macro marker already present")
        return

    pattern = re.compile(
        r"double r = ML_FMA\(-n, 0\.69314718036912381649, x\);\s*\n"
        r"\s*r = ML_FMA\(-n, 1\.90821490974462528503e-10, r\);"
    )

    replacement = f"""/* {MARKER_EXP}
     *
     * Use the canonical ln(2) split shared with ml_log().
     * The previous hardcoded low constant did not match ML_LN2_LO.
     */
    double r = ML_FMA(-n, ML_LN2_HI, x);
    r = ML_FMA(-n, ML_LN2_LO, r);"""

    patched, count = pattern.subn(lambda m: replacement, text, count=1)
    if count != 1:
        # If the old literal block is gone but the marker is present, okay.
        if "ML_FMA(-n, ML_LN2_HI, x)" in text:
            print(f"  [skip] {path}: exp already appears to use ML_LN2_HI.")
            return
        fail(f"{path}: could not find ml_exp Cody-Waite constant block.")

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
    print("  MATHLIB v12A1: PAYNE-HANEK + EXP LN2 SPLIT FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] payne_hanek.h — full-table sweep")
    patch_payne_hanek(v12, force)

    print("\n[2/2] exp_log.c — use ML_LN2_HI / ML_LN2_LO")
    patch_exp_constants(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Fixes applied.")
    print("")
    print("  Rebuild and rerun oracle:")
    print("")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check > /tmp/oracle_out2.txt || true")
    print("    grep -n FAIL /tmp/oracle_out2.txt")
    print("    tail -n 30 /tmp/oracle_out2.txt")
    print("")
    print("  Expected:")
    print("    - large sin/cos failures should collapse dramatically")
    print("    - gamma may still have residue; that becomes the next script")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
