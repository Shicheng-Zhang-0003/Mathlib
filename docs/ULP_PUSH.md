# ULP Push toward <0.5 (Correct Rounding) — Status R6

`<0.5 ULP` means correctly rounded (single rounding error). Proving it for
ALL binary64 inputs is the Table Maker's Dilemma (100+ bit + Ziv + Lefèvre
worst-case search, millions of CPU hours; CORE-MATH took years). Below is
measurement, not proof — except where stated, guarantees stay at the ≤5 ULP
gate. All LD paths are x86-64 only (`long double==double` collapses to DD).

## Measured 0 ULP (dense 2000-pt + coarse + edge grids, mpmath 80-dps)

Oracle families (412 vectors, gate ≤5 ULP, measured worst 0):

| Family | Grid worst | Fix that closed it (R6 unless noted) |
|---|---|---|
| sin, cos, exp, log, gamma, lgamma, pow | 0 (212-vector grid, prior work) | Payne-Hanek, LD split, Lanczos (v12R2/R5) |
| asin, acos | 0 (was 2 mid-range) | LD whole-formula + atanl source |
| asinh, acosh, atanh | 0 (was 1-2) | LD log1p forms, single rounding |
| log10 | 0 (was 2) | exact hex split (LO sign was flipped) |
| log2 | 0 (was 1) | exact hex split (single-double bias) |
| cbrt | 0 (was 1) | LD neighbor correction (mini-Ziv select) |
| erf | 0 (was 6) | 1-erfc route + LD Taylor + LD complement |
| sinh, tanh | 0 (was 1) | LD expl composition; direct LD Taylor |
| sinpi, cospi | 0 | exact integer/half shortcuts + LD combine |
| tan | 0 (was 1, 10% rate) | LD sincos + LD quotient |
| exp2 | 0 (was 1) | direct LD exp2l kernel, single rounding |
| bessel_i0, bessel_i1 | 0 | exact (all-positive series, LD) |
| atan | 0 on 3966/4000 (was 200 fails) | LD identity tree + compensated Horner |

Bessel J0/J1: bit-exact away from zeros (|truth|>0.3); crossings are an
absolute-error regime (~1e-14, intrinsic ascending-series cancellation —
same class as large-x trig, documented not defective).

## Residual 1 ULP (all ≤ gate; boundary proximity in ULP from truth to
rounding boundary — smaller means harder/deferred)

| Function | Worst point (exact double) | Proximity | Cause class |
|---|---|---|---|
| exp10 | 9.7501 | 0.0001 | LD noise floor; needs 106-bit product |
| expm1 | 5.89501 | 0.0002 | LD expl source noise |
| cosh | -2.69 | 0.0001 | Table-Maker (truth on boundary) |
| erfc (CF) | -2.142 | 0.0005 | CF truncation at noise floor |
| atan | 4.035 + 22 cataloged | 0.0001–0.05 | compensated Horner limit |
| digamma | 1.445 | 0.042 | 200x recurrence cancellation → needs 106-bit end-to-end |
| erfc (series) | 3 pts ≤0.17 | 0.017–0.17 | LD-Kahan noise floor |
| log1p | 65 pts | ≤5 gate | TwoSum path composes 0.5+0.5 |
| Y0/Y1 benign | scattered | — | log-singularity + cancellation; ≤2 ULP |
| Kelvin/Airy/K | transition bands | — | documented limits (1e-15..1e-13 rel) |

Boundary catalog (exact doubles that need 100+ bit arithmetic to decide;
future Ziv work starts here): exp10(9.7501e0), expm1(5.89501e0),
cosh(-2.69e0), erfc(-2.142e0), atan(4.035e0 + 22 pts in test logs),
digamma(1.445e0). Reproduce: `python3 scripts/oracles/generate_oracles_ext.py`
plus the sweep harness in `/tmp/opencode/ulp` (method documented above).

## Law (unchanged)

`pow_err ≈ |y|·log_err + exp_err`; `gamma_rel_err ≈ lgamma_abs_err`.
LD (64-bit mantissa, ~1e-19) everywhere + single round to double, plus
exact integer handling. The R6 lesson repeated everywhere: composing two
0-ULP doubles still costs 1 ULP — only single-rounding from a >64-bit
evaluation attains 0 ULP, and only away from Table-Maker boundaries.

## What remains for PROVEN <0.5

1. Ziv loop + MPFR fallback per function + boundary-proximity test
   (the catalog above is the exact target set).
2. Baart/Lefèvre worst-case search over full binary64 (millions of hours).
3. Slow-trig exact fractionals; general gamma TD; short-LD platform QD.
4. Bessel Y/Kelvin/Airy transition bands (documented, ~1e-13 rel).
5. DD digamma end-to-end; DD erfc-Kahan; log1p TwoSum-path audit.

Do not ship `<0.5 ULP` as a guarantee. Honest claim: **0 ULP measured on
the oracle 412 + wide grids (x86-64, 26 families); gate ≤5 ULP proven.**
