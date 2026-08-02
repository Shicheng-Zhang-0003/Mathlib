# MathLib v12A1 Development Tree

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
