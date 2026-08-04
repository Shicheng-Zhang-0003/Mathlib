#!/usr/bin/env python3
"""
17_fix_oracle_build.py
Run from the folder that CONTAINS the v12A1 working folder.

GIANT ERROR #1: the oracle verification path cannot be built reliably.

Root cause:
  build/libmathc.a is compiled with MATHLIB_SANITIZERS=ON, so every
  object contains ASan/UBSan instrumentation. Three different places
  then tell you to link oracle_check MANUALLY, and the instructions
  printed by the 15_*.py scripts omit the sanitizer flags. Result:
  hundreds of undefined references to __asan_* / __ubsan_* (error.txt).

Fix:
  1. CMakeLists.txt: oracle_check becomes a first-class, always-built
     target. It inherits add_link_options() automatically, so the
     flags can never mismatch again.
  2. closure_gate.sh step 6: use ./build/oracle_check (no gcc line).
  3. run_all_tests.py phase_oracle: use ./build/oracle_check.

Targets:
  v12A1/CMakeLists.txt
  v12A1/closure_gate.sh
  v12A1/run_all_tests.py

Usage:
  python3 17_fix_oracle_build.py
  python3 17_fix_oracle_build.py --force
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_ORACLE_BUILD_FIX"

NEW_CMAKE_BLOCK = """# MATHLIB_V12A1_ORACLE_BUILD_FIX
# oracle_check is a first-class target, always built.
# It inherits sanitizer link flags from add_link_options() when
# MATHLIB_SANITIZERS=ON, eliminating the manual link lines that
# produced undefined __asan_* / __ubsan_* references.
add_executable(oracle_check tests/test_oracle.c)
target_compile_definitions(oracle_check PRIVATE MATHLIB_HAS_ORACLE_DATA)
target_include_directories(oracle_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(oracle_check mathc m)
"""

NEW_GATE_STEP6 = """echo "[6/7] Running mpmath oracle validation..."
# MATHLIB_V12A1_ORACLE_BUILD_FIX: oracle_check is produced by the CMake
# build above and inherits sanitizer link flags. No manual gcc line.
./build/oracle_check"""

NEW_PHASE_ORACLE = r'''def phase_oracle(report: TestReport, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 6: MPMATH ORACLE VALIDATION")
    print("=" * 72)

    # MATHLIB_V12A1_ORACLE_BUILD_FIX
    # oracle_check is a first-class CMake target. It inherits the
    # sanitizer link flags from the global configuration, so the
    # manual gcc line (which caused undefined __asan_* references
    # when flags mismatched) is gone.
    oracle_binary = Path(BUILD_DIR) / "oracle_check"
    if not oracle_binary.exists():
        report.add(TestResult(
            "oracle_check", "Oracle", False, 0.0,
            error="Binary not found. Run the build phase first."
        ))
        return

    run_test(report, "oracle_check (mpmath ground truth)", "Oracle",
             [str(oracle_binary)], timeout=120,
             verbose=verbose, fail_fast=fail_fast)


'''


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
    if (root / "CMakeLists.txt").is_file() and (root / "src").is_dir():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def patch_cmake(v12: Path, force: bool) -> None:
    path = v12 / "CMakeLists.txt"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    pattern = re.compile(
        r"(?ms)^option\(MATHLIB_ORACLE_TESTS[^\n]*\)\s*\n"
        r"if\(MATHLIB_ORACLE_TESTS\).*?^endif\(\)\s*$"
    )
    patched, count = pattern.subn(
        lambda m: NEW_CMAKE_BLOCK.rstrip("\n"), text, count=1
    )
    if count != 1:
        fail(f"{path}: could not find MATHLIB_ORACLE_TESTS block. Source may have drifted.")
    write_text(path, patched)


def patch_gate(v12: Path, force: bool) -> None:
    path = v12 / "closure_gate.sh"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    pattern = re.compile(
        r'(?ms)^echo "\[6/7\] Running mpmath oracle validation\.\.\."\s*\n'
        r".*?\n\./build/oracle_check\s*$"
    )
    patched, count = pattern.subn(
        lambda m: NEW_GATE_STEP6, text, count=1
    )
    if count != 1:
        fail(f"{path}: could not find step [6/7] oracle block. Source may have drifted.")
    write_text(path, patched)


def patch_run_all_tests(v12: Path, force: bool) -> None:
    path = v12 / "run_all_tests.py"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    pattern = re.compile(
        r"(?ms)^def phase_oracle\(.*?(?=^# -{20,}\n# Phase 7)"
    )
    patched, count = pattern.subn(
        lambda m: NEW_PHASE_ORACLE, text, count=1
    )
    if count != 1:
        fail(f"{path}: could not find phase_oracle(). Source may have drifted.")
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
    print("  MATHLIB v12A1: FIX ORACLE BUILD (script 17)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/3] CMakeLists.txt — oracle_check as first-class target")
    patch_cmake(v12, force)

    print("\n[2/3] closure_gate.sh — step 6 uses the CMake binary")
    patch_gate(v12, force)

    print("\n[3/3] run_all_tests.py — phase_oracle uses the CMake binary")
    patch_run_all_tests(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Oracle build path fixed.")
    print("")
    print("  What changed:")
    print("    - oracle_check is always built by CMake")
    print("    - sanitizer link flags are inherited, never hand-typed")
    print("    - all manual gcc oracle link lines removed")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC \\")
    print("      -DCMAKE_BUILD_TYPE=Debug -DMATHLIB_SANITIZERS=ON")
    print("    cmake --build build")
    print("    ./build/oracle_check")
    print("")
    print("  Next: 18_minimax_truth.py")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
