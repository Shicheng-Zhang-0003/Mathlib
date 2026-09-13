# MathLib v12R2 Development Tree

v12R2 is the refinement cycle following the v12A1 architectural evolution.

v12A1 proved the architecture. v12R2 fixes critical bugs and improves accuracy.

<!-- MATHLIB_V12A1_DOCS_ALIGNMENT -->
v12R2 key fixes:
- **Removed x86 inline asm** — portable `__builtin_sqrt` with correct rounding
- **Fixed fmod** — proper IEEE-754 `x - trunc(x/y)*y` with error-free multiplication
- **Fixed atan/asin/acos** — argument reduction to |x| ≤ tan(π/12), stable formulas
- **Fixed tanh/atanh/asinh** — small-x threshold 1e-4 → 1.5e-8 (was millions of ULP error)
- **Fixed pow integer path** — limit 64 → 1023 (exact binary exponentiation)
- **Fixed gamma reflection** — uses full π (ML_PI_HI_D + ML_PI_LO_D)
- **Fixed variance** — Welford's online algorithm (no catastrophic cancellation)
- **ML_FMA software fallback** — proper Dekker Two-Product+Two-Sum
- **All decimal constants → exact hex floats** — ln2, π/√2, polynomial coefficients

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
<!-- MATHLIB_V12A1_README_STATUS_V2 -->
<!-- MATHLIB_V12A1_A1_FREEZE -->
A1 closure is **complete** with R2 refinements.

- Oracle validation: **212 passed, 0 failed** (all functions ≤ 5 ULP vs mpmath ground truth)
- Full test gauntlet: **31/32 passed** (modular, smoke, edge, fuzz, oracle, boundary)
- Closure gate: **PASSED**

See `docs/V12A1_ROADMAP.md` for the work plan.
See `release_notes.md` for the v12R2 refinement summary.
