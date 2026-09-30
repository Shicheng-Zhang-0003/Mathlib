# Maintainer Notes (v12R2 / V1.2-RC2)

Before accepting changes:

Check:
- API compatibility (`extern "C"`, status codes additive only)
- numerical behavior (exact-binary mpmath truth; report worst ULP + grid)
- compiler cleanliness (`-Werror`, 3 profiles, C++ header check)
- test coverage impact (modular + edge 23/23 + oracle + fuzz + sanitizers)
- documentation updates (tiers only; never quote smoke tolerances as ULP)
- canonical TU list still matches (CMake/Makefile/run_all_tests/edge, 34 TUs)

Numerical changes require justification and validation:
- state the error law (e.g. `pow_err ≈ |y|·log_err`),
- measure before/after ULP on exact-binary grids + oracle,
- hold the ≤5 ULP gate; 0-ULP measurements are not proofs (see `ULP_PUSH.md`).
