# Numerical Policy (v12R2 / V1.2-RC2)

Prioritizes:

1. predictable behavior
2. documented error bounds
3. reproducibility
4. portability

Floating point behavior follows IEEE-754 expectations where applicable.
Exact bit identity across all platforms is not guaranteed (x86-64 LD paths
vs short-LD fallback; guarantees stay at the ≤5 ULP gate).

Accuracy tiers are non-interchangeable (`PRECISION_CONTRACT.md`):
gate ≤5 ULP proven; 0 ULP measured 2026-09-30 on oracle 212 + wide grids
(x86-64 only) is not a proof for all inputs (see `ULP_PUSH.md`).
Ground truth must be exact-binary (integer-ratio) mpmath, never decimal.
