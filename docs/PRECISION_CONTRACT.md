# Precision Contract

v12R2 uses three validation tiers with **different, non-interchangeable**
tolerances. A passing test in one tier says nothing about the guarantees of
another. This table resolves the split.

| Tier | Ground truth | Gate | What it proves |
| :--- | :--- | :--- | :--- |
| **Oracle** (`tests/test_oracle.c` + `tests/oracle_data.h`) | mpmath at 80 decimal places | **<= 5 ULP** per call, all 7 families (sin, cos, exp, log, gamma, lgamma, pow), 212 vectors | Numerical accuracy of the validated core. This is the only tier that certifies precision. |
| **Smoke / modular** (`tests/test.c`, `test_core`, `test_trig`, `test_linalg`, `test_dsp`) | Self-agreement / mathematical identities (e.g. `fft_dc`, round-trip solves) | Loose absolute `ASSERT_NEAR` tolerances (`1e-9` … `1e-15`); at magnitude ~1, `1e-12` is ~4500 ULP | No-crash, no-NaN, plausible-output smoke coverage. A pass here is **not** a precision claim. |
| **Edge directed** (`tests/test_edge_*.c`, 22 suites / 441 assertions) | Exact IEEE-754 / boundary expectations (signed zero, Inf/NaN guards, `fmod` identities, domain clamps) | Exact or near-exact match on directed cases (`1e-15` on identities such as `atan2(±0,-1) == ±π`) | Boundary and special-value correctness. Covers branches the oracle grid never hits. |

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
   `pde`, `harmonic`, `mcmc`, `manifold`, `info`, `analytic_nt`, `control`)
   have **no oracle coverage** and inherit no precision guarantee until
   oracle vectors exist for them. See `docs/API_STATUS.md`.
