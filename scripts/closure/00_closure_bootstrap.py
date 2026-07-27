#!/usr/bin/env python3
"""
00_closure_bootstrap.py

Run this from the folder that CONTAINS the v11S working folder.

Expected layout:

    current_directory/
        00_closure_bootstrap.py
        v11S/
            src/
            include/
            tests/
            docs/
            ...

This script creates the v11S closure control system:

    v11S/docs/CLOSURE_PUNCHLIST.md
    v11S/docs/TEAM_DOCTRINE.md
    v11S/docs/SCRIPT_CHANGE_POLICY.md
    v11S/scripts/closure/README.md
    v11S/closure_gate.sh

It deliberately does NOT modify math source code yet.

Usage:

    python3 00_closure_bootstrap.py

Force overwrite:

    python3 00_closure_bootstrap.py --force
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


PUNCHLIST = r"""# MathLib v11S Closure Punchlist

This file is the authoritative short-term closure queue for v11S.

Current operational state:

> Late v11A3 / v11S release-candidate territory.
> Not yet a true stable v11S release.

---

## P0 Blockers

These must be fixed before v11S can be treated as stable.

### 1. `ml_pow` is not closure-grade

Current implementation is effectively:

    exp(y * log(x))

That is insufficient for:

- negative base with integer exponent,
- signed zero behavior,
- `0^0`,
- `1^NaN`,
- `NaN^0`,
- infinite exponent behavior,
- infinite base behavior.

Required action:

- add integer-exponent classification helpers,
- implement explicit special-case decision tree,
- preserve signed zero where appropriate,
- keep overflow/underflow behavior consistent with `ml_exp`.

---

### 2. `ml_cplx_power` lacks required special cases

Current implementation uses only:

    exp(b * log(a))

That cannot satisfy required special cases such as:

- `0^0 = 1`
- `1^NaN = 1`
- `NaN^0 = 1`
- `0^i = NaN`
- `0^-1 = inf`
- `inf^i = NaN`
- negative real principal-branch powers

Required action:

- add explicit complex special-case decision tree,
- preserve principal-branch behavior for finite nonzero bases,
- document conservative infinite-exponent contract.

---

### 3. `ml_cplx_arg` has a quadrant bug

The third-quadrant path currently divides by the wrong component.

Required action:

- replace hand-rolled logic with `ml_atan2(imag, real)`.

This also aligns implementation with the documented claim that complex arg is atan2-based.

---

### 4. `ml_exp` needs explicit NaN / infinity guards

`ml_exp` must classify NaN and infinity before range reduction.

Required action:

- return NaN for NaN input,
- return positive infinity for positive infinity,
- return positive zero for negative infinity,
- avoid casting NaN to integer.

---

### 5. `ml_log` must handle positive infinity

Current behavior can fall through into the finite path and produce NaN.

Required action:

- return positive infinity for positive infinity,
- preserve negative infinity as domain NaN,
- preserve zero as negative infinity.

---

### 6. Official verification must run edge tests

The edge suite exists but is not part of the main closure gate.

Required action:

- use `closure_gate.sh` as the strict gate,
- later fold the edge suite into the official verification path,
- do not claim closure while edge tests are skipped.

---

## P1 High-Priority Issues

These should be fixed soon after P0.

### 1. `ml_sinh` loses tiny inputs

Add a small-input branch consistent with other hyperbolics.

---

### 2. AVX2 batch `rsqrt` needs semantic guarding

The scalar fallback handles NaN / non-positive inputs.

The AVX2 path should either:

- fall back to scalar semantics for exceptional inputs,
- or explicitly document the restricted fast-path contract.

Preferred:

- guard the AVX2 path and fall back for non-finite / non-positive / subnormal inputs.

---

### 3. Matrix and tensor indexing should use wide size arithmetic internally

Public API may remain `int`, but internal offset arithmetic should avoid signed overflow.

Required action:

- convert matrix loops to size-based indexing internally,
- convert tensor offset arithmetic to wide size arithmetic.

---

### 4. CI should run deterministic fuzzing and edge tests

Required action:

- add edge tests to CI or create a dedicated closure gate job,
- use fixed fuzzer seeds in CI.

---

### 5. Documentation must match implementation

After source fixes:

- update `API_STATUS.md`,
- update `KNOWN_LIMITATIONS.md`,
- update `V11S_CLOSURE_SUMMARY.md`,
- ensure README verification flow matches actual verification scripts.

---

## Script Sequence

### Script 00

`00_closure_bootstrap.py`

Creates:

- closure punchlist,
- team doctrine,
- script-only change policy,
- strict closure gate.

This script does not modify math source.

---

### Script 01

`01_p0_exp_log_pow_complex.py`

Intended to fix:

- `src/internal/pow_util.h`
- `src/exp_log.c`
- `src/complex.c`

This is the first real source-change script.

---

### Script 02

`02_p1_simd_index_ci_docs.py`

Intended to fix:

- SIMD batch guarding,
- internal size arithmetic,
- CI determinism,
- documentation alignment.

---

## Closure Rule

v11S is not stable until:

1. all P0 items are fixed,
2. edge tests pass,
3. sanitizers pass,
4. documentation matches code,
5. strict closure gate passes.
"""


TEAM_DOCTRINE = r"""# MathLib Team Doctrine

This document defines how the multi-model development workflow operates.

The goal is not merely to generate code.

The goal is to generate, verify, understand, and own the system.

---

## Release Philosophy

All releases follow:

    S -> A1 -> A2 -> A3 -> S

Meaning:

- `S` = stable release
- `A1` = first major development wave
- `A2` = second major development wave
- `A3` = freeze, audit, bug discovery, hardening
- next `S` = promoted stable release

Current state:

> v11S is operationally between late A3 and true S.

Therefore:

- no new features,
- no new math families,
- no new public APIs,
- only correctness, tests, documentation, and closure hygiene.

---

## Blue Team

Owner: Qwen

Purpose:

- large-scale implementation,
- refactoring,
- test scaffolding,
- documentation generation,
- script-based codebase changes,
- turning contracts into code.

Rules:

- Blue implements against contracts.
- Blue does not silently invent new contracts.
- Blue preserves zero allocation, thread safety, and C99 compatibility.
- Blue changes are applied through single scripts whenever possible.

---

## Red Team

Members:

- ChatGPT
- Gemini
- GLM
- Kimi, when available

Purpose:

- mathematical correctness review,
- IEEE-754 special-case review,
- standards review,
- test-gap analysis,
- documentation/code mismatch detection.

Rules:

- Red does not generate large implementation patches unless explicitly asked.
- Red produces defect lists, severity rankings, and missing tests.
- Red disagreements become test cases.
- Red is invoked after pseudo-blue-red and before final release promotion.

---

## Yellow Team

Owner: Microsoft Copilot

Purpose:

- Windows portability,
- MSVC / clang-cl compatibility,
- CMake presets,
- install profiles,
- packaging,
- Windows CI,
- DLL / static library validation.

Rules:

- Yellow does not change numerical behavior.
- Yellow preserves the deterministic scientific contract.
- Yellow must validate behavior, not merely compilation.

---

## Violet Team

Owner: Grok

Purpose:

- extreme adversarial review,
- hostile input discovery,
- UB traps,
- security and robustness edge cases,
- fuzzer blind-spot discovery.

Rules:

- Violet does not propose feature creep during A3.
- Violet produces attack vectors, not roadmaps.
- Violet findings become tests or documented limitations.

---

## GitHub Copilot

Optional role:

- holistic reviewer,
- documentation drift checker,
- small regression-test writer,
- learning companion.

It is not the final mathematical authority.

---

## Decision Rules

1. Tests decide correctness.
2. If Red models disagree, convert the disagreement into a test.
3. If behavior cannot be fixed now, document it as a known limitation.
4. Documentation must never overclaim implementation.
5. During A3, correctness beats novelty.
6. A module is mathematically owned only when the orchestrator can explain its contract, failure modes, and tests.

---

## Handoff Packet

Before real Red Team review, prepare:

1. changed files,
2. contract decisions,
3. test evidence,
4. sanitizer output,
5. known limitations,
6. deferred items.

No tree should be handed to Red Team without these.
"""


SCRIPT_POLICY = r"""# MathLib Script-Only Change Policy

All future codebase changes conducted by Blue Team should be applied through a single script whenever practical.

The script is run from the folder that CONTAINS the `v11S` working folder.

Expected layout:

    parent_folder/
        some_change_script.py
        v11S/
            src/
            include/
            tests/
            docs/
            ...

---

## Goals

1. Atomicity
   A change set is applied as one logical operation.

2. Reproducibility
   The same script can be rerun on a clean tree.

3. Auditability
   The script itself documents what changed.

4. Safety
   Scripts should avoid silent manual edits.

5. Orchestrator control
   The human orchestrator remains in control of every change batch.

---

## Script Naming Convention

Scripts should be numbered:

    00_closure_bootstrap.py
    01_p0_exp_log_pow_complex.py
    02_p1_simd_index_ci_docs.py
    03_...

The number defines intended order.

---

## Script Requirements

Every script should:

- verify that it is run from the correct parent directory,
- verify that `v11S/` exists,
- print what it writes or modifies,
- avoid interactive prompts,
- be idempotent where possible,
- use `--force` if overwrite is desired,
- preserve LF line endings,
- avoid introducing fast-math, heap allocation, or global state.

---

## Source Modification Rules

Scripts that modify source code should:

- prefer small, explicit replacements,
- fail loudly if expected text is missing,
- not apply partial changes if a required anchor is absent,
- update tests and docs in the same change set when behavior changes,
- preserve public API signatures unless explicitly authorized.

---

## Forbidden During A3

Unless explicitly authorized:

- no new modules,
- no new public APIs,
- no new math families,
- no performance experiments,
- no speculative refactors,
- no feature creep.

A3 is for closure.
"""


CLOSURE_README = r"""# MathLib Closure Scripts

This directory contains scripts that apply controlled change sets to the v11S tree.

Scripts are run from the folder that CONTAINS `v11S`.

Example:

    python3 00_closure_bootstrap.py

---

## Current Script Sequence

### 00_closure_bootstrap.py

Creates closure control documents and strict closure gate.

Does not modify math source.

Creates:

    v11S/docs/CLOSURE_PUNCHLIST.md
    v11S/docs/TEAM_DOCTRINE.md
    v11S/docs/SCRIPT_CHANGE_POLICY.md
    v11S/scripts/closure/README.md
    v11S/closure_gate.sh

---

### 01_p0_exp_log_pow_complex.py

Intended next script.

Purpose:

- fix P0 exp/log/pow/complex blockers.

Likely targets:

    v11S/src/internal/pow_util.h
    v11S/src/exp_log.c
    v11S/src/complex.c

---

### 02_p1_simd_index_ci_docs.py

Intended follow-up script.

Purpose:

- SIMD semantic guarding,
- internal size arithmetic,
- CI determinism,
- documentation alignment.

---

## Rule

Do not manually drift the tree when a script can express the change atomically.
"""


GATE = r"""#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

echo "========================================================="
echo "  MATHLIB v11S: STRICT CLOSURE GATE"
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
echo "  STRICT CLOSURE GATE PASSED"
echo "========================================================="
"""


def write_file(path: Path, content: str, force: bool) -> bool:
    if path.exists() and not force:
        print(f"[skip] {path} already exists")
        return False

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)

    print(f"[write] {path}")
    return True


def main() -> int:
    force = "--force" in sys.argv[1:]

    root = Path.cwd()
    v11s = root / "v11S"

    if not v11s.is_dir():
        print("ERROR: Run this script from the folder that CONTAINS the v11S directory.")
        print("Expected to find: ./v11S")
        return 1

    docs = v11s / "docs"
    closure_scripts = v11s / "scripts" / "closure"

    print("=========================================================")
    print("  MathLib v11S Closure Bootstrap")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    write_file(docs / "CLOSURE_PUNCHLIST.md", PUNCHLIST, force)
    write_file(docs / "TEAM_DOCTRINE.md", TEAM_DOCTRINE, force)
    write_file(docs / "SCRIPT_CHANGE_POLICY.md", SCRIPT_POLICY, force)
    write_file(closure_scripts / "README.md", CLOSURE_README, force)

    gate_path = v11s / "closure_gate.sh"
    wrote_gate = write_file(gate_path, GATE, force)

    if wrote_gate or gate_path.exists():
        gate_path.chmod(0o755)
        if not wrote_gate:
            print(f"[chmod] {gate_path}")

    # Archive this script into the repository if possible.
    try:
        source_script = Path(__file__).resolve()
        archived_script = closure_scripts / "00_closure_bootstrap.py"

        if source_script != archived_script:
            if archived_script.exists() and not force:
                print(f"[skip] {archived_script} already exists")
            else:
                archived_script.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_script, archived_script)
                print(f"[archive] {archived_script}")
    except NameError:
        print("[note] Could not archive script because __file__ is unavailable.")

    print("---------------------------------------------------------")
    print("Bootstrap complete.")
    print("")
    print("Next recommended script:")
    print("    01_p0_exp_log_pow_complex.py")
    print("")
    print("Strict gate:")
    print("    cd v11S && ./closure_gate.sh")
    print("")
    print("NOTE: The strict gate is expected to FAIL until P0 source fixes are applied.")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
