
from __future__ import annotations

import re
import shutil
import sys
from datetime import date
from pathlib import Path

MARKER = "MATHLIB_V12A1_BOOTSTRAP"
TODAY = date.today().isoformat()

# ---------------------------------------------------------------------------
# New file contents
# ---------------------------------------------------------------------------

NEW_VERSION_H = r"""#ifndef MATHLIB_VERSION_H
#define MATHLIB_VERSION_H

/* MATHLIB_V12A1_BOOTSTRAP */
#define MATHLIB_VERSION_MAJOR 12
#define MATHLIB_VERSION_MINOR 1
#define MATHLIB_VERSION_PATCH 0
#define MATHLIB_VERSION_STRING "12.1.0-a1 (v12A1 development)"
#define MATHLIB_VERSION_TAG "v12A1"

#endif /* MATHLIB_VERSION_H */
"""

NEW_README = r"""# MathLib v12A1 Development Tree

v12A1 is the architectural evolution cycle following the v11S stable release.

v11S proved the foundations:
- deterministic C99 implementation
- strict compiler hygiene
- explicit numerical contracts
- zero-allocation workspace design
- reproducible validation workflow

v12A1 replaces approximations with the real thing:
- true minimax polynomials (replacing Taylor series)
- true Payne-Hanek range reduction (removing the 1e15 wall)
- Lanczos gamma function (replacing the degree-8 sketch)
- extended-precision pow
- error-free Cody-Waite reductions

## Build

```bash
cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC
cmake --build build
```

## Test

```bash
python3 run_all_tests.py
```

## Status

This is a development tree. Nothing here is release-grade yet.
See `docs/V12A1_ROADMAP.md` for the work plan.
"""

NEW_ROADMAP = r"""# v12A1 Development Roadmap

## Theme

v12A1 is the architectural evolution cycle.
v11S proved the foundations. v12A1 replaces approximations with the real thing.

## Bootstrap

- [x] Identity transition (this document created by 00_v12a1_bootstrap.py)
- [x] v11S closure documents archived
- [x] Version bumped to 12.1.0-a1
- [x] Banner strings updated

## Work Items

### 1. True Minimax Polynomials (P0)
- Run compute_minimax.py (it exists, it was never used)
- Replace Taylor coefficients in src/internal/minimax.h
- Target: true Remez or Chebyshev economized polynomials
- Validate: oracle ULP distance must not regress

### 2. Extended Range Reduction (P0)
- The 1e15 wall in payne_hanek.h is a domain clamp, not Payne-Hanek
- Implement true Payne-Hanek or extend Cody-Waite to full double range
- Remove the NaN return for sin(1e50)
- This is the single biggest limitation in v11S

### 3. Gamma Function Redesign (P0)
- Replace the rough degree-8 polynomial on [1,2]
- Implement Lanczos approximation (g=7, n=9)
- Add reflection formula for negative arguments
- Add ml_lgamma as a new API
- Target: <= 5 ULP like the rest of the transcendentals

### 4. Error-Free Cody-Waite in ml_exp (P0)
- Current: two separate rounded subtractions
- Fix: use ML_FMA for exact residual computation
- Or: 3-term split of ln(2)

### 5. Extended-Precision pow (P1)
- Split ml_log into high/low parts
- Compute y * log(x) with FMA
- Add integer-exponent fast path
- Add near-integer result detection

### 6. Word-at-a-Time fmod (P1)
- Current: O(quotient) loop, up to 2046 iterations
- Fix: process in 64-bit chunks

### 7. Iterative Refinement in Linear Algebra (P1)
- After LU solve: compute residual, solve correction, update
- Cost: one extra matvec + one extra triangular solve

### 8. Fixed-Point CORDIC Upgrade (P1)
- Extend from 16 to 24 iterations
- Extend atan table
- Tighten test tolerances

### 9. Better Fast-Math Polynomials (P2)
- ml_fast_log2: degree 3 -> degree 5
- ml_fast_exp2: degree 5 -> degree 7

### 10. SIMD Dispatch Evaluation (P2)
- Decision document: is the Quake rsqrt hack worth keeping?
- No code change unless decision is to replace or remove

## Not In Scope

- New math families (unless justified by existing module gaps)
- Performance experiments before correctness is established
- Feature creep during A1
- Mixed-radix FFT (deferred to v12A2 or later)
- Adaptive ODE solvers (deferred)

## Script Sequence

| #  | Script | Section |
|----|--------|---------|
| 00 | 00_v12a1_bootstrap.py | Identity (this script) |
| 01 | 01_minimax_pipeline.py | Minimax generation |
| 02 | 02_error_free_cleanup.py | FMA / error-free layer |
| 03 | 03_exp_cody_waite.py | Exp reduction fix |
| 04 | 04_log_reconstruction.py | Log reconstruction fix |
| 05 | 05_trig_minimax.py | Trig coefficient swap |
| 06 | 06_explog_minimax.py | Exp/log coefficient swap |
| 07 | 07_payne_hanek.py | True range reduction |
| 08 | 08_gamma_lanczos.py | Gamma redesign |
| 09 | 09_pow_extended.py | Extended-precision pow |
| 10 | 10_fmod_fast.py | Word-at-a-time fmod |
| 11 | 11_linalg_refinement.py | Iterative refinement |
| 12 | 12_cordic_24iter.py | CORDIC upgrade |
| 13 | 13_fastmath_polys.py | Fast-math polynomials |
| 14 | 14_simd_evaluation.py | SIMD decision doc |
| 15 | 15_oracle_expansion.py | Oracle test expansion |
| 16 | 16_closure_gate.py | v12A1 closure gate |

## Closure Rule

v12A1 is not stable until:
1. all P0 items are implemented,
2. oracle validation passes with <= 5 ULP,
3. edge tests pass,
4. sanitizers pass,
5. documentation matches code,
6. strict closure gate passes.
"""

NEW_SCRIPTS_README = r"""# MathLib v12A1 Development Scripts

This directory contains scripts that apply controlled change sets to the v12A1 tree.

Scripts are run from the folder that CONTAINS `v12A1`.

Example:

    python3 00_v12a1_bootstrap.py

---

## Script Sequence

See `docs/V12A1_ROADMAP.md` for the full table.

Scripts are numbered in dependency order.
Scripts 07, 08, 10, 11, 12, 14 can run in parallel after bootstrap.

---

## Rule

Do not manually drift the tree when a script can express the change atomically.
"""

NEW_CLOSURE_GATE = r"""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "========================================================="
echo "  MATHLIB v12A1: DEVELOPMENT GATE"
echo "========================================================="

CC="${CC:-gcc}"
SEED="${MATHLIB_ULTIMATE_SEED:-123456789}"

echo "[1/7] Configuring with ASan + UBSan..."
rm -rf build
cmake -B build \
    -DMATHLIB_PROFILE=SCIENTIFIC \
    -DCMAKE_BUILD_TYPE=Debug \
    -DMATHLIB_SANITIZERS=ON

echo "[2/7] Building..."
cmake --build build

echo "[3/7] Running modular tests..."
./build/test_core
./build/test_trig
./build/test_linalg
./build/test_dsp

echo "[4/7] Running edge tests with sanitizers..."
MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh

echo "[5/7] Running boundary gauntlet..."
./build/fuzz_boundary

echo "[6/7] Running mpmath oracle validation..."
"$CC" -std=c99 -O3 -fPIE \
    -fsanitize=address,undefined -fno-omit-frame-pointer \
    -Iinclude/mathlib -Isrc \
    -DMATHLIB_HAS_ORACLE_DATA \
    -o build/oracle_check \
    tests/test_oracle.c \
    -Lbuild -lmathc -lm
./build/oracle_check

echo "[7/7] Running ultimate fuzzer..."
"$CC" -std=c99 -O3 \
    -fsanitize=address,undefined -fno-omit-frame-pointer \
    -Iinclude/mathlib -Isrc \
    -o build/ultimate_fuzzer \
    tests/ultimate_fuzzer.c \
    -Lbuild -lmathc -lm
./build/ultimate_fuzzer "$SEED"

echo "========================================================="
echo "  v12A1 DEVELOPMENT GATE PASSED"
echo "========================================================="
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def read_text(path: Path) -> str:
    return normalize(path.read_text(encoding="utf-8"))


def locate_v12a1() -> tuple[Path, Path]:
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    # Convenience: allow running from inside v12A1 itself.
    if (root / "src" / "core.c").is_file() and (root / "include" / "mathlib").is_dir():
        print("  [note] Running from inside v12A1; treating current directory as v12A1.")
        return root.parent, root
    fail(
        "Run this script from the folder that CONTAINS the v12A1 directory, "
        "or from inside v12A1 itself."
    )


# ---------------------------------------------------------------------------
# Patch operations
# ---------------------------------------------------------------------------

def patch_version(v12: Path, force: bool) -> None:
    path = v12 / "include" / "mathlib" / "version.h"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = read_text(path)
    if MARKER in text and not force:
        print(f"  [skip] {path}: already bootstrapped")
        return
    write_text(path, NEW_VERSION_H)


def patch_cmake(v12: Path, force: bool) -> None:
    path = v12 / "CMakeLists.txt"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = read_text(path)
    if "mathlib_v12A1" in text and not force:
        print(f"  [skip] {path}: already bootstrapped")
        return
    patched = text.replace(
        "project(mathlib_v11S VERSION 11.0.0 LANGUAGES C)",
        "project(mathlib_v12A1 VERSION 12.1.0 LANGUAGES C)",
    )
    if patched == text:
        fail(f"{path}: could not find project() line to patch. Source may have drifted.")
    write_text(path, patched)


def archive_v11s_docs(v12: Path, force: bool) -> None:
    archive_dir = v12 / "docs" / "archive" / "v11S"
    docs_to_archive = [
        "CLOSURE_PUNCHLIST.md",
        "CLOSURE_P0_APPLIED.md",
        "CLOSURE_P1_APPLIED.md",
        "V11S_CLOSURE_SUMMARY.md",
        "V11S_VS_V12A1.md",
    ]
    archive_dir.mkdir(parents=True, exist_ok=True)
    for name in docs_to_archive:
        src = v12 / "docs" / name
        dst = archive_dir / name
        if not src.exists():
            print(f"  [skip] {src}: not found (already archived or never existed)")
            continue
        if dst.exists() and not force:
            print(f"  [skip] {dst}: already archived")
            continue
        shutil.move(str(src), str(dst))
        print(f"  [archive] {src.name} -> docs/archive/v11S/")


def patch_release_notes(v12: Path, force: bool) -> None:
    path = v12 / "RELEASE_NOTES.md"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = read_text(path)
    if "Shipped" in text and not force:
        print(f"  [skip] {path}: already finalized")
        return
    new_text = text.rstrip() + f"""

## Shipped

v11S closure gate passed. Promoted to stable.
Ship date: {TODAY}

---

*This file is the v11S ship record. v12A1 development notes live in
`docs/V12A1_ROADMAP.md`.*
"""
    write_text(path, new_text)


def write_readme(v12: Path, force: bool) -> None:
    path = v12 / "README.md"
    if path.exists() and not force:
        text = read_text(path)
        if "v12A1 Development Tree" in text:
            print(f"  [skip] {path}: already bootstrapped")
            return
    write_text(path, NEW_README)


def write_roadmap(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "V12A1_ROADMAP.md"
    if path.exists() and not force:
        print(f"  [skip] {path}: already exists")
        return
    write_text(path, NEW_ROADMAP)


def patch_team_doctrine(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "TEAM_DOCTRINE.md"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = read_text(path)
    if "v12A1 early A1" in text and not force:
        print(f"  [skip] {path}: already updated")
        return

    old_state = """Current state:
> v11S is operationally between late A3 and true S.

Therefore:
- no new features,
- no new math families,
- no new public APIs,
- only correctness, tests, documentation, and closure hygiene."""

    new_state = """Current state:
> v12A1 early A1 development wave.
> v11S shipped. Closure documents archived.

Therefore:
- new algorithms are in scope,
- approximation redesigns are the primary mission,
- correctness contracts still apply,
- script-only change policy still applies,
- zero allocation and thread safety still apply."""

    if old_state not in text:
        # Try a more lenient match
        pattern = re.compile(
            r"Current state:\s*\n"
            r"> v11S is operationally between late A3 and true S\.\s*\n\s*\n"
            r"Therefore:\s*\n"
            r"- no new features,\s*\n"
            r"- no new math families,\s*\n"
            r"- no new public APIs,\s*\n"
            r"- only correctness, tests, documentation, and closure hygiene\."
        )
        patched, count = pattern.subn(new_state, text, count=1)
        if count != 1:
            fail(f"{path}: could not find 'Current state' block. Source may have drifted.")
        write_text(path, patched)
    else:
        write_text(path, text.replace(old_state, new_state))


def patch_banners(v12: Path, force: bool) -> None:
    """Replace v11S banner strings with v12A1 in test/fuzz/soak infrastructure."""

    replacements = [
        # (relative_path, old_string, new_string)
        ("soak_test.sh",
         "MATHLIB v11S: DETERMINISTIC SOAK TEST",
         "MATHLIB v12A1: DETERMINISTIC SOAK TEST"),
        ("closure_gate.sh",
         None,  # full rewrite
         None),
        ("run_all_tests.py",
         "MATHLIB v11S: FULL TEST SUITE (THE NUCLEAR OPTION)",
         "MATHLIB v12A1: FULL TEST SUITE (THE NUCLEAR OPTION)"),
        ("run_all_tests.py",
         "MATHLIB v11S: FULL TEST SUITE RESULTS",
         "MATHLIB v12A1: FULL TEST SUITE RESULTS"),
        ("run_all_tests.py",
         "ALL TESTS PASSED. v11S is clean.",
         "ALL TESTS PASSED. v12A1 is clean."),
        ("run_all_tests.py",
         "MathLib v11S: Run every test in the tree.",
         "MathLib v12A1: Run every test in the tree."),
        ("run_all_tests.py",
         "Run from INSIDE the v11S folder:",
         "Run from INSIDE the v12A1 folder:"),
        ("run_all_tests.py",
         "cd v11S",
         "cd v12A1"),
        ("run_all_tests.py",
         'ERROR: Run this script from INSIDE the v11S folder.',
         'ERROR: Run this script from INSIDE the v12A1 folder.'),
        ("tests/run_tests.py",
         "MATHLIB v11S MODULAR CI/CD TEST RUNNER",
         "MATHLIB v12A1 MODULAR CI/CD TEST RUNNER"),
        ("tests/run_soak.sh",
         "MATHLIB v11S: DETERMINISTIC SOAK TEST",
         "MATHLIB v12A1: DETERMINISTIC SOAK TEST"),
        ("tests/ultimate_fuzzer.c",
         "MATHLIB v11S: THE ULTIMATE FUZZER (ASan + UBSan)",
         "MATHLIB v12A1: THE ULTIMATE FUZZER (ASan + UBSan)"),
        ("tests/fuzz_boundary_gauntlet.c",
         "MATHLIB v11S: BOUNDARY & INVARIANT GAUNTLET",
         "MATHLIB v12A1: BOUNDARY & INVARIANT GAUNTLET"),
        ("tests/test_oracle.c",
         "MATHLIB v11S: ORACLE VALIDATION (mpmath ground truth)",
         "MATHLIB v12A1: ORACLE VALIDATION (mpmath ground truth)"),
        ("tests/test_oracle.c",
         "MATHLIB v11S: ORACLE VALIDATION (SKIPPED)",
         "MATHLIB v12A1: ORACLE VALIDATION (SKIPPED)"),
        ("tests/test.c",
         "MathLib v11S: Monolithic Smoke Test",
         "MathLib v12A1: Monolithic Smoke Test"),
    ]

    # Group by file to minimize reads/writes
    file_patches: dict[str, list[tuple[str, str]]] = {}
    for rel_path, old, new in replacements:
        if old is None:
            continue  # handled separately (full rewrite)
        file_patches.setdefault(rel_path, []).append((old, new))

    for rel_path, patches in file_patches.items():
        path = v12 / rel_path
        if not path.is_file():
            print(f"  [skip] {path}: not found")
            continue
        text = read_text(path)
        changed = False
        for old, new in patches:
            if old in text:
                text = text.replace(old, new)
                changed = True
        if changed:
            write_text(path, text)
        else:
            print(f"  [skip] {path}: banners already updated")

    # Full rewrite of closure_gate.sh
    gate_path = v12 / "closure_gate.sh"
    if gate_path.is_file():
        text = read_text(gate_path)
        if "v12A1: DEVELOPMENT GATE" in text and not force:
            print(f"  [skip] {gate_path}: already updated")
        else:
            write_text(gate_path, NEW_CLOSURE_GATE)
            gate_path.chmod(0o755)


def create_scripts_dir(v12: Path, force: bool) -> None:
    scripts_dir = v12 / "scripts" / "v12a1"
    readme = scripts_dir / "README.md"
    if readme.exists() and not force:
        print(f"  [skip] {readme}: already exists")
        return
    scripts_dir.mkdir(parents=True, exist_ok=True)
    write_text(readme, NEW_SCRIPTS_README)


def archive_self(v12: Path, force: bool) -> None:
    try:
        source_script = Path(__file__).resolve()
        archived_script = v12 / "scripts" / "v12a1" / source_script.name
        if source_script == archived_script:
            return
        if archived_script.exists() and not force:
            print(f"  [skip] {archived_script}: already archived")
            return
        archived_script.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_script, archived_script)
        print(f"  [archive] {archived_script}")
    except NameError:
        print("  [note] Could not archive script because __file__ is unavailable.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()

    print("=========================================================")
    print("  MATHLIB v12A1: IDENTITY BOOTSTRAP")
    print("=========================================================")
    print(f"  Root:   {root}")
    print(f"  v12A1:  {v12}")
    print(f"  Force:  {force}")
    print(f"  Date:   {TODAY}")
    print("---------------------------------------------------------")

    print("\n[1/9] Version identity")
    patch_version(v12, force)

    print("\n[2/9] CMakeLists.txt")
    patch_cmake(v12, force)

    print("\n[3/9] Archive v11S closure documents")
    archive_v11s_docs(v12, force)

    print("\n[4/9] Finalize v11S release notes")
    patch_release_notes(v12, force)

    print("\n[5/9] README.md")
    write_readme(v12, force)

    print("\n[6/9] v12A1 roadmap")
    write_roadmap(v12, force)

    print("\n[7/9] Team Doctrine")
    patch_team_doctrine(v12, force)

    print("\n[8/9] Banner strings")
    patch_banners(v12, force)

    print("\n[9/9] scripts/v12a1/ directory")
    create_scripts_dir(v12, force)
    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Bootstrap complete.")
    print("")
    print("  Verification:")
    print("")
    print("    cd v12A1")
    print("    grep MATHLIB_VERSION_STRING include/mathlib/version.h")
    print("    grep 'project(' CMakeLists.txt")
    print("    ls docs/archive/v11S/")
    print("    ls docs/V12A1_ROADMAP.md")
    print("    ls scripts/v12a1/")
    print("")
    print("  Build check:")
    print("")
    print("    cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC")
    print("    cmake --build build")
    print("    ./build/test_core")
    print("")
    print("  Next script:")
    print("")
    print("    01_minimax_pipeline.py")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
