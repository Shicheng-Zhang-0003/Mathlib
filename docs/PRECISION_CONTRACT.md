# Precision Contract

The v12R2 core measurements and current v12R3 working tree use validation
tiers with **different, non-interchangeable**
tolerances. A passing test in one tier says nothing about the guarantees of
another. This table resolves the split.

Gate stays **≤5 ULP** on the oracle (pass/fail). The 2026-10-19 rerun has
406 mpmath reference vectors across 27 families plus 21 large-argument
domain checks: all 427 checks pass, and the worst reference error is 1 ULP
at `gamma(50)`. This is a measurement, not a proof for all binary64 inputs
(Table Maker's Dilemma; extended-precision paths are platform-dependent; see
`docs/ULP_PUSH.md`). Do not ship `<0.5 ULP` as a guarantee until per-function
Ziv + worst-case search land.

| Tier | Ground truth | Gate | What it proves |
| :--- | :--- | :--- | :--- |
| **Oracle** (`tests/test_oracle.c` + `tests/oracle_data.h` + `tests/oracle_data_ext.h`) | mpmath at 50 decimal places (`scripts/oracles/generate_oracles{,_ext}.py: mp.dps=50`; sufficient for double ULP work) | **<= 5 ULP** per call, 406 reference vectors across 27 families (sin, cos, exp, log, gamma, lgamma, pow; asin, acos, asinh, acosh, atanh, atan, log10, log2, cbrt, erf, sinh, tanh, sinpi, cospi, tan, exp2, bessel_i0, bessel_i1, bessel_j0, bessel_j1), plus 21 large-argument domain checks | Numerical accuracy of the validated core. This is the only tier that certifies precision. |
| **Smoke / modular** (`tests/test.c`, `test_core`, `test_trig`, `test_linalg`, `test_dsp`) | Self-agreement / mathematical identities (e.g. `fft_dc`, round-trip solves) | Loose absolute `ASSERT_NEAR` tolerances (`1e-9` … `1e-15`); at magnitude ~1, `1e-12` is ~4500 ULP | No-crash, no-NaN, plausible-output smoke coverage. A pass here is **not** a precision claim. |
| **Edge directed** (`tests/test_edge_*.c`) | Exact IEEE-754 / boundary expectations (signed zero, Inf/NaN guards, remainder identities, domain clamps) | Exact or near-exact match on directed cases | Boundary and special-value regression coverage. It is not an accuracy certificate. |

Rules:

1. Only the oracle tier may be cited for ULP claims: the current run has 406
   reference vectors within ≤5 ULP, plus 21 large-argument domain checks;
   its maximum reference error is 1 ULP at `gamma(50)`.
2. Smoke/modular tolerances exist to catch regressions and crashes, not to
   bound error. Do not tighten them into pseudo-oracles and do not quote
   them as accuracy figures.
3. `tests/test_edge_accuracy_audit.c` is a **regression guard**, not an
   oracle-tier certification. It pins mpmath-derived references with ULP
   tolerances for functions outside the oracle vector set (atan, expm1,
   log1p, erfc, the eight integer-order Bessel functions, Airy, digamma,
   Jacobi, exp10, cosh/sinh). Its tolerances are chosen to catch regressions, not to
   certify precision; the K and Airy transition bands use loose guards
   (~1e-8 relative) because those are documented limits, not defects.
4. New Batch-1 modules (`optim_n`, `ode_sys`, `spectral`, `stats_inv`, `sde`,
    `pde`, `harmonic`, `mcmc`, `manifold`, `info`, `analytic_nt`, `control`,
    `kelvin`), the contest module (`ml_contest.h`, 42 functions: u64
    exact-or-sentinel paths are edge-tier exact like rule 3's directed
    cases; LD series paths are regression guards), the applied modules
    (`ml_geometry.h`, `ml_graph.h`, `ml_modeling.h`, olympiad-NT annex:
    exact-or-sentinel paths are edge-tier exact; iterative solvers are
    regression guards with fail-loud CONVERGENCE), and the doctorate
    modules (`ml_measure.h`, `ml_funcanal.h`, `ml_topology.h`,
    `ml_diffgeom.h`, `ml_algebraic.h`: exact integer paths are edge-tier
    exact; LD/FP paths are regression guards with fail-loud codes)
    have **no oracle coverage** and inherit no precision guarantee until
    oracle vectors exist for them. See `docs/API_STATUS.md` and
    `docs/CONTEST_GUIDE.md`.
   The eleven v12R3 extension modules (`bigint`, `math_extra`, `special`,
   `calculus_extra`, `complexan`, `geom_extra`, `graph_extra`, `linalg_extra`,
   `ode_extra`, `optim_extra`, and `stats_extra`) also have no oracle-tier
   certification. Their directed tests establish only the pinned cases.
5. Despot audit (2026-09-30): targeted regression `/tmp/opencode/mathlib-work/
   despot_check` pins the fixed semantics (exp2 subnormal, crt2 24/30,
   cross_entropy +Inf, ber/bei NaN, sphere unit gates, fixed-point rounding).
   It is a contract guard, not a ULP certification.

Round-6 ULP push: `tests/test_oracle.c` grew 19 families / 200 vectors
(`tests/oracle_data_ext.h` via `scripts/oracles/generate_oracles_ext.py`;
Bessel J vectors are benign-only by generator-enforced |truth| > 0.3,
crossings stay an absolute-error regime per KNOWN_LIMITATIONS).
The 2026-10-19 R7 update adds 15 atan references. The current oracle contains
406 mpmath vectors and 21 large-argument domain checks; 427/427 checks pass,
with a maximum of 1 ULP at `gamma(50)`.

Round-5 (2026-10-10, Edmonton): `tests/test_edge_despot_r5.c` (42 assertions)
pins the R5 programming/robustness fixes at regression-guard tier (rule 3):
CPU dispatch, fixed-point, cross3/transpose aliasing, egcd, status precision,
DoS caps, RNG warmup, ODE arrival, Hurwitz guard + closed-form physics
(SUVAT/projectile/SHM/exp). Oracle stays 212/212 at 0 ULP measured (gate ≤5
ULP). See `docs/DESPOT_R5_AUDIT.md`.

Round-4 (2026-10-07): `tests/test_edge_despot_r4.c` (21 assertions) pins the
Round-4 fixes at regression-guard tier (rule 3): zeta(0), Hurwitz EM 0 ULP,
clog/csqrt/cpow, remainder bit-exact vs system, gamma P/Q at Inf. Oracle stays
212/212 at 0 ULP measured (gate ≤5 ULP). See `docs/DESPOT_R4_AUDIT.md`.

Round-3 (2026-10-01): the K-quadrature bridge (Bessel K 4<=x<16, Airy 2.5<x<8.5) and the Kelvin DLMF full sums are new numeric methods with no oracle vectors; they are pinned by `tests/test_edge_accuracy_audit.c` at 65536 ULP (a regression guard, per rule 3), NOT by the oracle tier. `accuracy_audit` grew to 371 assertions. Measured vs mpmath 80-dps: K 6.8e-15, Ai 2.1e-13, Kelvin ~1e-13.\n\nRound-2 (2026-10-01): zeta/kelvin/prime_pi/acosh/softplus/normal_inv/hurwitz/catalan fixes re-measured 0 ULP on oracle 212 + wide grids (x86-64); gate stays ≤5 ULP. Batch-1 still no oracle.

Round-7 (2026-10-19): `tests/test_edge_despot_r7.c` adds 973 directed and
reference regression assertions for extension APIs, including general-order
Bessel, Airy, Jacobi, Dawson, Struve, spherical functions, and polygamma.
These assertions are regression guards only. The 48-test full gate passed;
the core oracle tier and its ≤5 ULP gate are unchanged. See
`docs/DESPOT_R7_AUDIT.md` for sampled mpmath evidence and documented domains.

Exact foundations (2026-10-19): `ml_mpz_t`, `ml_mpq_t`, `ml_polyq_t`, and
`ml_gfp_ctx_t` use integer algorithms and do not inherit floating-point ULP
claims. Their contracts are exact within available memory; the test suite
cross-checks small results against native integer arithmetic and checks large
division/multiplication identities. `ml_interval_t` encloses basic arithmetic
results by outward widening under its explicit IEEE binary64/round-to-nearest
precondition. Rational-to-interval conversion uses exact comparisons and
returns either a point or adjacent-double enclosure. This is a conditional
arithmetic guarantee, not a proof of the platform compiler's conformance or
of any transcendental function.
