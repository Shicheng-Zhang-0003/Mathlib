# Precision Contract

v12R2 uses three validation tiers with **different, non-interchangeable**
tolerances. A passing test in one tier says nothing about the guarantees of
another. This table resolves the split.

Gate stays **≤5 ULP** on the oracle (pass/fail). Measured accuracy on
2026-09-30 is **worst 0 ULP (212/212 correctly rounded)** on the grid vs
mpmath 50-dps exact-binary truth — a measurement, not a proof for all
binary64 inputs (Table Maker's Dilemma; x86-64 LD paths only, see
`docs/ULP_PUSH.md`). Do not ship `<0.5 ULP` as a guarantee until per-function
Ziv + worst-case search land.

| Tier | Ground truth | Gate | What it proves |
| :--- | :--- | :--- | :--- |
| **Oracle** (`tests/test_oracle.c` + `tests/oracle_data.h`) | mpmath at 50 decimal places (`scripts/oracles/generate_oracles.py: mp.dps=50`; sufficient for double ULP work) | **<= 5 ULP** per call, all 7 families (sin, cos, exp, log, gamma, lgamma, pow), 212 vectors | Numerical accuracy of the validated core. This is the only tier that certifies precision. |
| **Smoke / modular** (`tests/test.c`, `test_core`, `test_trig`, `test_linalg`, `test_dsp`) | Self-agreement / mathematical identities (e.g. `fft_dc`, round-trip solves) | Loose absolute `ASSERT_NEAR` tolerances (`1e-9` … `1e-15`); at magnitude ~1, `1e-12` is ~4500 ULP | No-crash, no-NaN, plausible-output smoke coverage. A pass here is **not** a precision claim. |
| **Edge directed** (`tests/test_edge_*.c`, 23 suites / 700+ assertions) | Exact IEEE-754 / boundary expectations (signed zero, Inf/NaN guards, `fmod` identities, domain clamps) | Exact or near-exact match on directed cases (`1e-15` on identities such as `atan2(±0,-1) == ±π`) | Boundary and special-value correctness. Covers branches the oracle grid never hits. |

Rules:

1. Only the oracle tier may be cited for ULP claims (see `release_notes.md`:
   212 passed, 0 failed, all ≤ 5 ULP).
2. Smoke/modular tolerances exist to catch regressions and crashes, not to
   bound error. Do not tighten them into pseudo-oracles and do not quote
   them as accuracy figures.
3. `tests/test_edge_accuracy_audit.c` (v12R2) is a **regression guard**, not an
   oracle-tier certification. It pins mpmath-derived references with ULP
   tolerances for functions outside the 212-vector oracle grid (erfc, the
   eight integer-order Bessel functions, Airy, digamma, Jacobi, exp10,
   cosh/sinh). Its tolerances are chosen to catch regressions, not to
   certify precision; the K and Airy transition bands use loose guards
   (~1e-8 relative) because those are documented limits, not defects.
4. New Batch-1 modules (`optim_n`, `ode_sys`, `spectral`, `stats_inv`, `sde`,
    `pde`, `harmonic`, `mcmc`, `manifold`, `info`, `analytic_nt`, `control`,
    `kelvin`)
    have **no oracle coverage** and inherit no precision guarantee until
    oracle vectors exist for them. See `docs/API_STATUS.md`.
5. Despot audit (2026-09-30): targeted regression `/tmp/opencode/mathlib-work/
   despot_check` pins the fixed semantics (exp2 subnormal, crt2 24/30,
   cross_entropy +Inf, ber/bei NaN, sphere unit gates, fixed-point rounding).
   It is a contract guard, not a ULP certification.

Round-3 (2026-10-01): the K-quadrature bridge (Bessel K 4<=x<16, Airy 2.5<x<8.5) and the Kelvin DLMF full sums are new numeric methods with no oracle vectors; they are pinned by `tests/test_edge_accuracy_audit.c` at 65536 ULP (a regression guard, per rule 3), NOT by the oracle tier. `accuracy_audit` grew to 371 assertions. Measured vs mpmath 80-dps: K 6.8e-15, Ai 2.1e-13, Kelvin ~1e-13.\n\nRound-2 (2026-10-01): zeta/kelvin/prime_pi/acosh/softplus/normal_inv/hurwitz/catalan fixes re-measured 0 ULP on oracle 212 + wide grids (x86-64); gate stays ≤5 ULP. Batch-1 still no oracle.
