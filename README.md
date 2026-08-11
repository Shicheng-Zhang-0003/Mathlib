# MathLib v12A1 Development Tree

v12A1 is the architectural evolution cycle following the v11S stable release.

v11S proved the foundations:
- deterministic C99 implementation
- strict compiler hygiene
- explicit numerical contracts
- zero-allocation workspace design
- reproducible validation workflow

<!-- MATHLIB_V12A1_DOCS_ALIGNMENT -->
v12A1 replaces approximations with the real thing:
- validated Maclaurin kernels (minimax swap deferred to v12A2)
- true Payne-Hanek range reduction (removing the 1e15 wall)
- hybrid gamma: Lanczos DD + Stirling DD + exact half-integers
  + 1-step recurrence for x < 0.5 (<=5 ULP vs mpmath oracle)
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
<!-- MATHLIB_V12A1_README_STATUS_V2 -->
<!-- MATHLIB_V12A1_A1_FREEZE -->
A1 closure is **complete**.

- Oracle validation: **212 passed, 0 failed** (all functions <= 5 ULP vs mpmath ground truth)
- Full test gauntlet: **32/32 passed** (modular, smoke, edge, fuzz, oracle, boundary)
- Closure gate: **PASSED**

See `docs/V12A1_ROADMAP.md` for the work plan.
See `release_notes.md` for the v12A1 closure summary.
