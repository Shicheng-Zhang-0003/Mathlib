#!/usr/bin/env python3
"""
26_a1_gamma_recurrence_depth.py

Run from the folder that CONTAINS the v12A1 working folder.

Fixes the final gamma/lgamma residue (gamma(0.001) 11 ULP, gamma(-0.9) 6 ULP).

Root cause:
    For x near 1 (e.g. x=0.001 after one recurrence step, or x=1.9 in the
    reflection path), the Lanczos sum Ag = c0 + sum(c_i/(z+i)) has large
    alternating terms that cancel, losing precision.

Fix:
    Use 8 recurrence steps for x < 8, so the Lanczos approximation always
    sees z >= 7, where the terms are smaller and cancellation is reduced.

        lgamma(x) = lgamma(x+8) - sum( log(x+i) for i in 0..7 )

    The product/sum is accumulated in double-double, so no precision is lost.

Targets:
    v12A1/src/integral.c

Usage:
    python3 26_a1_gamma_recurrence_depth.py
    python3 26_a1_gamma_recurrence_depth.py --force
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH"

NEW_LGAMMA_POSITIVE_DD = '''static ml_dd_t ml_lgamma_positive_dd(double x) {
    /* MATHLIB_V12A1_GAMMA_RECURRENCE_DEPTH */
    /*
     * For x < 8, use 8 recurrence steps to shift the argument
     * into the range [8, 16], where the Lanczos approximation
     * is more accurate and less prone to cancellation.
     *
     * Gamma(x+8) = Gamma(x) * x * (x+1) * ... * (x+7)
     * lgamma(x)  = lgamma(x+8) - sum(log(x+i), i=0..7)
     */
    if (x < 8.0) {
        ml_dd_t L = ml_lgamma_lanczos_dd(x + 8.0);
        for (int i = 0; i < 8; i++) {
            ml_dd_t log_term = ml_log_dd(x + (double)i);
            L = ml_dd_sub(L, log_term);
        }
        return L;
    }
    return ml_lgamma_lanczos_dd(x);
}
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
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


def find_function_range(text: str, func_signature: str) -> tuple[int, int]:
    """Find start/end of a function via brace counting."""
    start = text.find(func_signature)
    if start == -1:
        return -1, -1
    brace_start = text.find("{", start)
    if brace_start == -1:
        return -1, -1
    depth = 0
    i = brace_start
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    return -1, -1


def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    func_sig = "static ml_dd_t ml_lgamma_positive_dd(double x)"
    start, end = find_function_range(text, func_sig)
    if start == -1:
        fail(f"{path}: could not find ml_lgamma_positive_dd function")

    new_text = text[:start] + NEW_LGAMMA_POSITIVE_DD + text[end:]
    write_text(path, new_text)


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
    print("  MATHLIB v12A1: GAMMA RECURRENCE DEPTH FIX")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/1] integral.c — increase recurrence depth to 8")
    patch_integral(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Recurrence depth fix applied.")
    print("")
    print("  What changed:")
    print("    - ml_lgamma_positive_dd now uses 8 recurrence steps for x < 8")
    print("    - Lanczos always sees z >= 7, reducing term cancellation")
    print("    - Fixes gamma(0.001) and gamma(-0.9) residue")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/oracle_check > /tmp/oracle_out5.txt || true")
    print("    grep -n FAIL /tmp/oracle_out5.txt")
    print("    tail -n 20 /tmp/oracle_out5.txt")
    print("")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Expected: oracle and edge suites collapse to zero failures.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
