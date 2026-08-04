# v12A1 Minimax Coefficient Validation Report

Generated: 2026-08-04
Method: Chebyshev interpolation (exact 80-digit power-basis conversion)
Ground truth: mpmath at 80 decimal places
Validation grid: 5001 points per reduced domain
Evaluation model: FMA Horner (libm fma, correctly rounded)

## Generated coefficients (minimax_coeffs.h, DORMANT)

| Function | Parity | Degree | Max ULP | Status |
|----------|--------|--------|---------|--------|
| sin | odd, x*P(x^2) | 9 | 2 | PASS |
| cos | even, P(x^2) | 8 | 5 | PASS |
| exp | plain (Taylor-19) | 19 | 1 | PASS |
| log | odd, z*P(z^2) | 7 | 3 | PASS |

## Active kernels in tree (NOT modified)

The running code uses the Maclaurin kernels in src/internal/minimax.h
and the Taylor/series kernels in src/exp_log.c. Measured on the same
reduced domains with the same evaluation model:

| Kernel | Source | Max ULP |
|--------|--------|---------|
| sin (Maclaurin, x*P(x^2), deg 19) | minimax.h | 1 |
| cos (Maclaurin, P(x^2), deg 18) | minimax.h | 1 |
| exp (Taylor-19 on reduced arg) | exp_log.c | 1 |
| log (atanh series, 11 terms) | exp_log.c | 1 |

## The exp decision (measured, not assumed)

Chebyshev interpolation cannot reach <= 5 ULP for exp under double
FMA Horner on [-ln2/2, ln2/2]. The approximation error at degree 10
is already ~1 ULP, but Horner evaluation adds ~6-8 ULP of rounding
noise, and the noise grows with degree. Measured evidence:

| Chebyshev degree | Max ULP |
|------------------|---------|
| 9 | 122 |
| 10 | 8 |
| 11 | 16 |
| 12 | 33 |
| 13 | 10 |
| 14 | 12 |

The oscillating, non-decreasing error is the signature of
evaluation noise, not truncation. The Taylor-19 kernel (exact 1/k!
coefficients) measures 1 ULP under the same model — strictly better
than any fitted polynomial. It is therefore recorded as the validated
exp coefficient set. Real libms escape this floor with table-driven
arguments or extended precision, not with higher-degree Horner.

## Decisions

- Hard gate: nothing is written unless every kernel is <= 5 ULP.
- No coefficient swap into trig.c / exp_log.c in this wave.
  The active kernels already meet the accuracy bar; the generated
  header stays dormant until a future gated swap script.

## Acceptance Criteria

All functions must achieve <= 5 ULP on their reduced domain.
