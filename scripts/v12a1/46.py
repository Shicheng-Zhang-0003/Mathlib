#!/usr/bin/env python3
"""
46_a1_docs_alignment.py

Documentation alignment for v12A1 A1 closure.

Fixes:
  1. README.md: "true minimax polynomials" claim is wrong.
     The minimax coefficients are DORMANT. Active kernels are
     still Maclaurin. Also updates gamma description.

  2. API_STATUS.md: ml_integral.h says "EXPERIMENTAL" and
     "Positive-domain gamma only." Now has full gamma/lgamma
     with Lanczos, Stirling, half-integers, 1-step recurrence,
     and reflection formula.

  3. KNOWN_LIMITATIONS.md: Only 4 lines. Needs v12A1 context:
     A1 freeze, dormant minimax, 1e15 wall removal.

  4. CLOSURE_PUNCHLIST.md: P0 items lack DONE markers.

Targets:
  v12A1/README.md
  v12A1/docs/API_STATUS.md
  v12A1/docs/KNOWN_LIMITATIONS.md
  v12A1/docs/archive/v11S/CLOSURE_PUNCHLIST.md

Usage:
  python3 46_a1_docs_alignment.py
  python3 46_a1_docs_alignment.py --force
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_DOCS_ALIGNMENT"


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


# ----------------------------------------------------------------
# Fix 1: README.md
# ----------------------------------------------------------------
def patch_readme(v12: Path, force: bool) -> None:
    path = v12 / "README.md"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    old_block = (
        "v12A1 replaces approximations with the real thing:\n"
        "- true minimax polynomials (replacing Taylor series)\n"
        "- true Payne-Hanek range reduction (removing the 1e15 wall)\n"
        "- Lanczos gamma function (replacing the degree-8 sketch)\n"
        "- extended-precision pow\n"
        "- error-free Cody-Waite reductions"
    )

    new_block = (
        "<!-- " + MARKER + " -->\n"
        "v12A1 replaces approximations with the real thing:\n"
        "- validated Maclaurin kernels (minimax swap deferred to v12A2)\n"
        "- true Payne-Hanek range reduction (removing the 1e15 wall)\n"
        "- hybrid gamma: Lanczos DD + Stirling DD + exact half-integers\n"
        "  + 1-step recurrence for x < 0.5 (<=5 ULP vs mpmath oracle)\n"
        "- extended-precision pow\n"
        "- error-free Cody-Waite reductions"
    )

    if old_block not in text:
        fail(
            f"{path}: could not find the v12A1 feature list.\n"
            f"  Expected: 'v12A1 replaces approximations with the real thing:'\n"
            f"  Source may have drifted."
        )

    text = text.replace(old_block, new_block, 1)
    write_text(path, text)


# ----------------------------------------------------------------
# Fix 2: API_STATUS.md
# ----------------------------------------------------------------
def patch_api_status(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "API_STATUS.md"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    old_line = (
        "| `ml_integral.h` | **EXPERIMENTAL** | Positive-domain gamma only; "
        "traditional integrator is experimental |"
    )

    new_line = (
        "| `ml_integral.h` | **STABLE** | "
        "<!-- " + MARKER + " --> "
        "Full gamma/lgamma: Lanczos DD, Stirling DD, exact half-integers, "
        "1-step recurrence (x<0.5), reflection formula, <=5 ULP oracle-validated; "
        "traditional integrator remains experimental |"
    )

    if old_line not in text:
        fail(
            f"{path}: could not find the ml_integral.h row.\n"
            f"  Expected: '| `ml_integral.h` | **EXPERIMENTAL** | Positive-domain gamma only; ...'\n"
            f"  Source may have drifted."
        )

    text = text.replace(old_line, new_line, 1)
    write_text(path, text)


# ----------------------------------------------------------------
# Fix 3: KNOWN_LIMITATIONS.md
# ----------------------------------------------------------------
def patch_known_limitations(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "KNOWN_LIMITATIONS.md"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    new_content = (
        "# Known Limitations\n"
        "\n"
        "<!-- " + MARKER + " -->\n"
        "\n"
        "v12A1 intentionally does not claim:\n"
        "\n"
        "- universal correctly rounded transcendental functions\n"
        "- identical output across every architecture\n"
        "- replacement of specialized high precision libraries\n"
        "\n"
        "These limitations are design choices, not hidden defects.\n"
        "\n"
        "## v12A1 Specific\n"
        "\n"
        "- **Minimax coefficients are DORMANT.** The generated coefficients in\n"
        "  `src/internal/minimax_coeffs.h` are validated but not active.\n"
        "  The running kernels use Maclaurin (Taylor) series, which already\n"
        "  meet the <=5 ULP oracle gate. The minimax swap is deferred to v12A2.\n"
        "\n"
        "- **The 1e15 wall is removed.** Payne-Hanek V6 handles the full\n"
        "  double range up to ~1.8e308. `sin(1e50)` and `cos(1e300)` produce\n"
        "  finite, correct results.\n"
        "\n"
        "- **A1 closure freeze is in effect.** No new modules, APIs, or math\n"
        "  families are permitted. Only closure fixes, tests, validation,\n"
        "  and documentation alignment are allowed.\n"
        "\n"
        "- **Gamma uses Lanczos g=7 n=9.** This coefficient set has ~1e-15\n"
        "  intrinsic approximation error. Half-integers bypass Lanczos entirely\n"
        "  via exact product/sum formulas. The 1-step recurrence for x < 0.5\n"
        "  avoids the cancellation that plagued the old multi-step recurrence.\n"
        "\n"
        "- **`ml_integral_traditional` remains experimental.** It is a simple\n"
        "  Riemann sum integrator and is not part of the validated numerical\n"
        "  core.\n"
    )

    write_text(path, new_content)


# ----------------------------------------------------------------
# Fix 4: CLOSURE_PUNCHLIST.md
# ----------------------------------------------------------------
def patch_closure_punchlist(v12: Path, force: bool) -> None:
    path = v12 / "docs" / "archive" / "v11S" / "CLOSURE_PUNCHLIST.md"
    if not path.is_file():
        fail(f"Missing: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Add a closure banner at the top, after the title
    old_header = "# MathLib v11S Closure Punchlist\n"
    new_header = (
        "# MathLib v11S Closure Punchlist\n"
        "\n"
        "<!-- " + MARKER + " -->\n"
        "> **STATUS: ALL P0 AND P1 ITEMS RESOLVED.**\n"
        "> v11S shipped 2026-08-02. All blockers below were fixed,\n"
        "> tested, and validated before promotion to stable.\n"
        "> This file is archived for historical reference.\n"
    )

    if old_header not in text:
        fail(
            f"{path}: could not find the title line.\n"
            f"  Expected: '# MathLib v11S Closure Punchlist'\n"
            f"  Source may have drifted."
        )

    text = text.replace(old_header, new_header, 1)

    # Mark each P0 item as DONE
    p0_items = [
        ("### 1. `ml_pow` is not closure-grade",
         "### 1. `ml_pow` is not closure-grade ✅ DONE"),
        ("### 2. `ml_cplx_power` lacks required special cases",
         "### 2. `ml_cplx_power` lacks required special cases ✅ DONE"),
        ("### 3. `ml_cplx_arg` has a quadrant bug",
         "### 3. `ml_cplx_arg` has a quadrant bug ✅ DONE"),
        ("### 4. `ml_exp` needs explicit NaN / infinity guards",
         "### 4. `ml_exp` needs explicit NaN / infinity guards ✅ DONE"),
        ("### 5. `ml_log` must handle positive infinity",
         "### 5. `ml_log` must handle positive infinity ✅ DONE"),
        ("### 6. Official verification must run edge tests",
         "### 6. Official verification must run edge tests ✅ DONE"),
    ]

    for old, new in p0_items:
        if old in text:
            text = text.replace(old, new, 1)

    # Mark P1 items as DONE
    p1_items = [
        ("### 1. `ml_sinh` loses tiny inputs",
         "### 1. `ml_sinh` loses tiny inputs ✅ DONE"),
        ("### 2. AVX2 batch `rsqrt` needs semantic guarding",
         "### 2. AVX2 batch `rsqrt` needs semantic guarding ✅ DONE"),
        ("### 3. Matrix and tensor indexing should use wide size arithmetic internally",
         "### 3. Matrix and tensor indexing should use wide size arithmetic internally ✅ DONE"),
        ("### 4. CI should run deterministic fuzzing and edge tests",
         "### 4. CI should run deterministic fuzzing and edge tests ✅ DONE"),
        ("### 5. Documentation must match implementation",
         "### 5. Documentation must match implementation ✅ DONE"),
    ]

    for old, new in p1_items:
        if old in text:
            text = text.replace(old, new, 1)

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
    print("  MATHLIB v12A1: DOCUMENTATION ALIGNMENT")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/4] README.md — fix minimax claim, update gamma description")
    patch_readme(v12, force)

    print("\n[2/4] API_STATUS.md — ml_integral.h EXPERIMENTAL -> STABLE")
    patch_api_status(v12, force)

    print("\n[3/4] KNOWN_LIMITATIONS.md — expand for v12A1")
    patch_known_limitations(v12, force)

    print("\n[4/4] CLOSURE_PUNCHLIST.md — mark P0/P1 items as DONE")
    patch_closure_punchlist(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Documentation alignment applied.")
    print("")
    print("  What changed:")
    print("    README.md:")
    print("      - 'true minimax polynomials' -> 'validated Maclaurin")
    print("        kernels (minimax swap deferred to v12A2)'")
    print("      - 'Lanczos gamma function' -> 'hybrid gamma: Lanczos DD")
    print("        + Stirling DD + exact half-integers + 1-step recurrence'")
    print("")
    print("    API_STATUS.md:")
    print("      - ml_integral.h: EXPERIMENTAL -> STABLE")
    print("      - 'Positive-domain gamma only' -> full gamma/lgamma description")
    print("")
    print("    KNOWN_LIMITATIONS.md:")
    print("      - Added v12A1-specific limitations section")
    print("      - Documented dormant minimax, removed 1e15 wall,")
    print("        A1 freeze, Lanczos g=7 n=9 intrinsic error")
    print("")
    print("    CLOSURE_PUNCHLIST.md:")
    print("      - Added closure banner: ALL P0 AND P1 ITEMS RESOLVED")
    print("      - Marked all 6 P0 items as DONE")
    print("      - Marked all 5 P1 items as DONE")
    print("")
    print("  Verify:")
    print("    Review the four files for accuracy.")
    print("    No build or test changes required.")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
