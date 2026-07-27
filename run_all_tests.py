#!/usr/bin/env python3
"""
run_all_tests.py
================
The Nuclear Option.

Runs EVERY test, fuzzer, oracle, and gauntlet in the v11S tree.
No shortcuts. No skips. Everything passes or you know exactly what didn't.

Run from INSIDE the v11S folder:

    cd v11S
    python3 run_all_tests.py

Options:
    --no-build        Skip the build phase (use existing binaries)
    --no-soak         Skip the 10,000-iteration soak test (default: skipped)
    --soak            Include the soak test
    --soak-n N        Set soak iteration count (default: 10000)
    --seed N          Override fuzzer seed (default: 123456789)
    --jobs N          Parallel build jobs (default: auto)
    --verbose         Show full output from every test (default: summary only)
    --fail-fast       Stop at the first failure instead of running everything

Exit codes:
    0  All tests passed
    1  One or more tests failed
    2  Build failed
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_SEED = 123456789
DEFAULT_SOAK_ITERATIONS = 10000
BUILD_DIR = "build"

CMAKE_CONFIGURE_CMD = [
    "cmake",
    "-B", BUILD_DIR,
    "-DMATHLIB_PROFILE=SCIENTIFIC",
    "-DCMAKE_BUILD_TYPE=Debug",
    "-DMATHLIB_SANITIZERS=ON",
]

EDGE_TEST_CFLAGS = [
    "-std=c99", "-O3",
    "-Wall", "-Wextra", "-Wconversion", "-Wshadow", "-Wpedantic", "-Werror",
    "-fno-fast-math", "-ffp-contract=off",
    "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
    "-Iinclude/mathlib", "-Isrc",
]

LIB_SOURCES = [
    "src/core.c",
    "src/trig.c",
    "src/exp_log.c",
    "src/complex.c",
    "src/linalg.c",
    "src/fft.c",
    "src/cpu_dispatch.c",
    "src/combinatorics.c",
    "src/quadratics.c",
    "src/polynomial.c",
    "src/numerical.c",
    "src/statistics.c",
    "src/integral.c",
    "src/ode.c",
    "src/optimization.c",
    "src/quaternion.c",
    "src/fixed_point.c",
]


# ---------------------------------------------------------------------------
# Result tracking
# ---------------------------------------------------------------------------

@dataclass
class TestResult:
    name: str
    phase: str
    passed: bool
    duration: float
    output: str = ""
    error: str = ""


@dataclass
class TestReport:
    results: list[TestResult] = field(default_factory=list)
    build_ok: bool = False

    def add(self, result: TestResult) -> None:
        self.results.append(result)

    @property
    def total(self) -> int:
        return len(self.results)

    @property
    def passed_count(self) -> int:
        return sum(1 for r in self.results if r.passed)

    @property
    def failed_count(self) -> int:
        return sum(1 for r in self.results if not r.passed)

    @property
    def total_duration(self) -> float:
        return sum(r.duration for r in self.results)

    def print_summary(self) -> None:
        print()
        print("=" * 72)
        print("  MATHLIB v11S: FULL TEST SUITE RESULTS")
        print("=" * 72)

        current_phase = ""
        for r in self.results:
            if r.phase != current_phase:
                current_phase = r.phase
                print(f"\n  [{current_phase}]")
            status = "PASS" if r.passed else "FAIL"
            icon = "\u2705" if r.passed else "\u274c"
            print(f"    {icon} {status}  {r.name:<52s} {r.duration:7.2f}s")

        print()
        print("-" * 72)
        print(f"  Total: {self.total}  |  "
              f"Passed: {self.passed_count}  |  "
              f"Failed: {self.failed_count}  |  "
              f"Time: {self.total_duration:.1f}s")
        print("-" * 72)

        if self.failed_count == 0:
            print("  \U0001f389 ALL TESTS PASSED. v11S is clean.")
        else:
            print("  \u26a0\ufe0f  FAILURES DETECTED:")
            for r in self.results:
                if not r.passed:
                    print(f"    \u274c [{r.phase}] {r.name}")
                    if r.error:
                        for line in r.error.strip().splitlines()[:5]:
                            print(f"       {line}")
        print("=" * 72)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_cmd(
    cmd: list[str],
    cwd: Optional[Path] = None,
    timeout: int = 600,
    env_extra: Optional[dict] = None,
) -> tuple[int, str, str]:
    """Run a command, return (returncode, stdout, stderr)."""
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    try:
        proc = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
        return proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return -1, "", f"TIMEOUT after {timeout}s"
    except FileNotFoundError:
        return -2, "", f"Command not found: {cmd[0]}"


def run_test(
    report: TestReport,
    name: str,
    phase: str,
    cmd: list[str],
    cwd: Optional[Path] = None,
    timeout: int = 600,
    env_extra: Optional[dict] = None,
    verbose: bool = False,
    fail_fast: bool = False,
) -> bool:
    """Run a single test, record result, return True if passed."""
    t0 = time.time()
    rc, stdout, stderr = run_cmd(cmd, cwd=cwd, timeout=timeout, env_extra=env_extra)
    duration = time.time() - t0

    passed = (rc == 0)
    result = TestResult(
        name=name,
        phase=phase,
        passed=passed,
        duration=duration,
        output=stdout,
        error=stderr if not passed else "",
    )
    report.add(result)

    if verbose or not passed:
        icon = "\u2705" if passed else "\u274c"
        print(f"  {icon} {name} ({duration:.2f}s)")
        if not passed:
            if stdout.strip():
                for line in stdout.strip().splitlines()[-10:]:
                    print(f"     stdout: {line}")
            if stderr.strip():
                for line in stderr.strip().splitlines()[-10:]:
                    print(f"     stderr: {line}")

    if not passed and fail_fast:
        print(f"\n  \u274c FAIL-FAST: Stopping at {name}")
        report.print_summary()
        sys.exit(1)

    return passed


# ---------------------------------------------------------------------------
# Phase 1: Build
# ---------------------------------------------------------------------------

def phase_build(report: TestReport, verbose: bool, jobs: Optional[int]) -> bool:
    print("\n" + "=" * 72)
    print("  PHASE 1: BUILD (CMake + ASan + UBSan)")
    print("=" * 72)

    # Configure
    t0 = time.time()
    rc, stdout, stderr = run_cmd(CMAKE_CONFIGURE_CMD)
    if rc != 0:
        report.build_ok = False
        print(f"  \u274c CMake configure FAILED")
        print(stderr[-2000:] if stderr else stdout[-2000:])
        return False
    print(f"  \u2705 CMake configured ({time.time() - t0:.1f}s)")

    # Build
    build_cmd = ["cmake", "--build", BUILD_DIR]
    if jobs:
        build_cmd += ["-j", str(jobs)]
    t0 = time.time()
    rc, stdout, stderr = run_cmd(build_cmd, timeout=300)
    if rc != 0:
        report.build_ok = False
        print(f"  \u274c Build FAILED")
        print(stderr[-3000:] if stderr else stdout[-3000:])
        return False
    report.build_ok = True
    print(f"  \u2705 Build succeeded ({time.time() - t0:.1f}s)")
    return True


# ---------------------------------------------------------------------------
# Phase 2: Modular CI Tests
# ---------------------------------------------------------------------------

def phase_modular(report: TestReport, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 2: MODULAR CI TESTS")
    print("=" * 72)

    tests = ["test_core", "test_trig", "test_linalg", "test_dsp"]
    for t in tests:
        binary = Path(BUILD_DIR) / t
        if not binary.exists():
            report.add(TestResult(t, "Modular CI", False, 0.0, error="Binary not found"))
            continue
        run_test(report, t, "Modular CI", [str(binary)],
                 verbose=verbose, fail_fast=fail_fast)


# ---------------------------------------------------------------------------
# Phase 3: Monolithic Smoke Test
# ---------------------------------------------------------------------------

def phase_smoke(report: TestReport, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 3: MONOLITHIC SMOKE TEST")
    print("=" * 72)

    binary = Path(BUILD_DIR) / "test"
    if not binary.exists():
        report.add(TestResult("test (smoke)", "Smoke", False, 0.0, error="Binary not found"))
        return
    run_test(report, "test (monolithic smoke)", "Smoke", [str(binary)],
             verbose=verbose, fail_fast=fail_fast)


# ---------------------------------------------------------------------------
# Phase 4: Edge Tests (compile + run each test_edge_*.c)
# ---------------------------------------------------------------------------

def phase_edge(report: TestReport, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 4: EDGE TESTS (compile + run)")
    print("=" * 72)

    cc = os.environ.get("CC", "gcc")
    edge_sources = sorted(glob.glob("tests/test_edge_*.c"))

    if not edge_sources:
        print("  \u26a0\ufe0f  No edge test sources found.")
        return

    os.makedirs(BUILD_DIR, exist_ok=True)

    for src in edge_sources:
        name = Path(src).stem
        binary = Path(BUILD_DIR) / name

        # Compile
        compile_cmd = [cc] + EDGE_TEST_CFLAGS + ["-o", str(binary), src] + LIB_SOURCES + ["-lm"]
        t0 = time.time()
        rc, stdout, stderr = run_cmd(compile_cmd, timeout=120)
        compile_time = time.time() - t0

        if rc != 0:
            report.add(TestResult(
                f"{name} (compile)",
                "Edge Tests",
                False,
                compile_time,
                error=stderr[-1000:],
            ))
            if fail_fast:
                print(f"  \u274c FAIL-FAST: {name} failed to compile")
                report.print_summary()
                sys.exit(1)
            continue

        # Run
        run_test(report, name, "Edge Tests", [str(binary)],
                 verbose=verbose, fail_fast=fail_fast)


# ---------------------------------------------------------------------------
# Phase 5: Fuzzers
# ---------------------------------------------------------------------------

def phase_fuzzers(report: TestReport, seed: int, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 5: FUZZERS (deterministic seeds)")
    print("=" * 72)

    # God Mode Fuzzer
    god_mode = Path(BUILD_DIR) / "fuzz_god_mode"
    if god_mode.exists():
        run_test(
            report,
            f"fuzz_god_mode (seed={seed})",
            "Fuzzers",
            [str(god_mode), str(seed)],
            timeout=300,
            verbose=verbose,
            fail_fast=fail_fast,
        )
    else:
        report.add(TestResult("fuzz_god_mode", "Fuzzers", False, 0.0, error="Binary not found"))

    # Boundary Gauntlet
    boundary = Path(BUILD_DIR) / "fuzz_boundary"
    if boundary.exists():
        run_test(
            report,
            "fuzz_boundary (boundary gauntlet)",
            "Fuzzers",
            [str(boundary)],
            timeout=120,
            verbose=verbose,
            fail_fast=fail_fast,
        )
    else:
        report.add(TestResult("fuzz_boundary", "Fuzzers", False, 0.0, error="Binary not found"))

    # Ultimate Fuzzer
    ultimate = Path(BUILD_DIR) / "ultimate_fuzzer"
    if not ultimate.exists():
        # Compile it manually if CMake didn't produce it
        cc = os.environ.get("CC", "gcc")
        ultimate_src = "tests/ultimate_fuzzer.c"
        if Path(ultimate_src).exists():
            compile_cmd = (
                [cc, "-std=c99", "-O3",
                 "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
                 "-Iinclude/mathlib", "-Isrc",
                 "-o", str(ultimate), ultimate_src]
                + LIB_SOURCES + ["-lm"]
            )
            rc, _, stderr = run_cmd(compile_cmd, timeout=120)
            if rc != 0:
                report.add(TestResult(
                    "ultimate_fuzzer (compile)", "Fuzzers", False, 0.0, error=stderr[-500:]
                ))

    if ultimate.exists():
        run_test(
            report,
            f"ultimate_fuzzer (seed={seed})",
            "Fuzzers",
            [str(ultimate), str(seed)],
            timeout=300,
            verbose=verbose,
            fail_fast=fail_fast,
        )
    else:
        report.add(TestResult("ultimate_fuzzer", "Fuzzers", False, 0.0, error="Binary not found"))


# ---------------------------------------------------------------------------
# Phase 6: Oracle Validation
# ---------------------------------------------------------------------------

def phase_oracle(report: TestReport, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 6: MPMATH ORACLE VALIDATION")
    print("=" * 72)

    cc = os.environ.get("CC", "gcc")
    oracle_binary = Path(BUILD_DIR) / "oracle_check"
    oracle_src = "tests/test_oracle.c"

    if not Path(oracle_src).exists():
        report.add(TestResult("oracle_check", "Oracle", False, 0.0, error="Source not found"))
        return

    # Compile with oracle data
    compile_cmd = (
        [cc, "-std=c99", "-O3", "-fPIE",
         "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
         "-Iinclude/mathlib", "-Isrc",
         "-DMATHLIB_HAS_ORACLE_DATA",
         "-o", str(oracle_binary), oracle_src,
         f"-L{BUILD_DIR}", "-lmathc", "-lm"]
    )
    t0 = time.time()
    rc, stdout, stderr = run_cmd(compile_cmd, timeout=120)
    if rc != 0:
        report.add(TestResult(
            "oracle_check (compile)", "Oracle", False,
            time.time() - t0, error=stderr[-500:]
        ))
        return

    run_test(report, "oracle_check (mpmath ground truth)", "Oracle",
             [str(oracle_binary)], timeout=120,
             verbose=verbose, fail_fast=fail_fast)


# ---------------------------------------------------------------------------
# Phase 7: Soak Test (optional)
# ---------------------------------------------------------------------------

def phase_soak(
    report: TestReport,
    seed: int,
    iterations: int,
    verbose: bool,
    fail_fast: bool,
) -> None:
    print("\n" + "=" * 72)
    print(f"  PHASE 7: SOAK TEST ({iterations} iterations)")
    print("=" * 72)

    god_mode = Path(BUILD_DIR) / "fuzz_god_mode"
    if not god_mode.exists():
        report.add(TestResult("soak_test", "Soak", False, 0.0, error="fuzz_god_mode not found"))
        return

    t0 = time.time()
    failed_at = -1
    for i in range(1, iterations + 1):
        iter_seed = seed + i
        rc, stdout, stderr = run_cmd(
            [str(god_mode), str(iter_seed)],
            timeout=60,
        )
        if rc != 0:
            failed_at = i
            report.add(TestResult(
                f"soak_test (failed at iter {i}, seed {iter_seed})",
                "Soak",
                False,
                time.time() - t0,
                error=stderr[-500:] if stderr else stdout[-500:],
            ))
            if fail_fast:
                break
            return

        if i % 1000 == 0:
            elapsed = time.time() - t0
            print(f"    [SOAK] {i}/{iterations} passed ({elapsed:.0f}s elapsed)")

    if failed_at < 0:
        report.add(TestResult(
            f"soak_test ({iterations} iterations)",
            "Soak",
            True,
            time.time() - t0,
        ))
        print(f"  \u2705 Soak test passed: {iterations} iterations, 0 failures")


# ---------------------------------------------------------------------------
# Phase 8: Python test runner (run_tests.py)
# ---------------------------------------------------------------------------

def phase_python_runner(report: TestReport, verbose: bool, fail_fast: bool) -> None:
    print("\n" + "=" * 72)
    print("  PHASE 8: PYTHON TEST RUNNER (run_tests.py)")
    print("=" * 72)

    runner = Path("tests/run_tests.py")
    if not runner.exists():
        print("  \u26a0\ufe0f  tests/run_tests.py not found, skipping.")
        return

    # run_tests.py expects to be run from build/
    run_test(
        report,
        "run_tests.py (modular runner)",
        "Python Runner",
        [sys.executable, os.path.join("..", "tests", "run_tests.py")],
        cwd=Path(BUILD_DIR),
        timeout=120,
        verbose=verbose,
        fail_fast=fail_fast,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="MathLib v11S: Run every test in the tree.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--no-build", action="store_true",
                        help="Skip build phase (use existing binaries)")
    parser.add_argument("--soak", action="store_true",
                        help="Include the soak test (default: skipped)")
    parser.add_argument("--no-soak", action="store_true",
                        help="Explicitly skip soak test")
    parser.add_argument("--soak-n", type=int, default=DEFAULT_SOAK_ITERATIONS,
                        help=f"Soak iteration count (default: {DEFAULT_SOAK_ITERATIONS})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"Fuzzer seed (default: {DEFAULT_SEED})")
    parser.add_argument("--jobs", "-j", type=int, default=None,
                        help="Parallel build jobs")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show output from every test")
    parser.add_argument("--fail-fast", action="store_true",
                        help="Stop at first failure")
    args = parser.parse_args()

    # Verify we're in the right place
    if not Path("src/core.c").exists() or not Path("include/mathlib").is_dir():
        print("ERROR: Run this script from INSIDE the v11S folder.")
        print("Expected to find: ./src/core.c and ./include/mathlib/")
        return 2

    print("=" * 72)
    print("  MATHLIB v11S: FULL TEST SUITE (THE NUCLEAR OPTION)")
    print("=" * 72)
    print(f"  Seed:       {args.seed}")
    print(f"  Soak:       {'YES (' + str(args.soak_n) + ' iters)' if args.soak else 'NO'}")
    print(f"  Build:      {'SKIP' if args.no_build else 'YES'}")
    print(f"  Verbose:    {args.verbose}")
    print(f"  Fail-fast:  {args.fail_fast}")
    print(f"  Started:    {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 72)

    report = TestReport()
    t_start = time.time()

    # Phase 1: Build
    if not args.no_build:
        if not phase_build(report, args.verbose, args.jobs):
            print("\n  \u274c Build failed. Cannot run tests.")
            return 2
    else:
        report.build_ok = True
        print("\n  [SKIP] Build phase (--no-build)")

    # Phase 2: Modular CI
    phase_modular(report, args.verbose, args.fail_fast)

    # Phase 3: Smoke test
    phase_smoke(report, args.verbose, args.fail_fast)

    # Phase 4: Edge tests
    phase_edge(report, args.verbose, args.fail_fast)

    # Phase 5: Fuzzers
    phase_fuzzers(report, args.seed, args.verbose, args.fail_fast)

    # Phase 6: Oracle
    phase_oracle(report, args.verbose, args.fail_fast)

    # Phase 7: Soak (optional)
    if args.soak and not args.no_soak:
        phase_soak(report, args.seed, args.soak_n, args.verbose, args.fail_fast)
    else:
        print("\n  [SKIP] Soak test (use --soak to enable)")

    # Phase 8: Python runner
    phase_python_runner(report, args.verbose, args.fail_fast)

    # Summary
    total_time = time.time() - t_start
    report.print_summary()
    print(f"\n  Wall-clock time: {total_time:.1f}s")
    print(f"  Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    return 0 if report.failed_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
