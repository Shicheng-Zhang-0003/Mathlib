# Testing Strategy (v12R2 / V1.2-RC2)

Validation layers (tiers are non-interchangeable, see `PRECISION_CONTRACT.md`):

1. Unit / modular (`test_core/trig/linalg/dsp`, `test.c` smoke 30/0)
2. Numerical edge (23/23 suites, 700+ assertions, exact IEEE-754 expectations)
3. Oracle comparison (212/212 at 0 ULP measured vs mpmath 50-dps
   exact-binary truth; gate stays ≤5 ULP — measurement, not proof)
4. Wide grids (pow 100 cases, sin/cos 25 benign + 2000 random, gamma 11-point;
   exact-binary truth; all core 0 ULP except large-x trig near-zero
   absolute-tiny and vs-sys residual 1)
5. Boundary fuzzing (25/0) + god-mode fuzz (61393/0, `--seed 123456789`)
6. Long-running soak (`--soak`, 10k iterations, off by default)
7. Regression checks (`accuracy_audit` pins, `despot_check` 17 contracts,
   `oracle_ulp1` 0-ULP scan)
8. Sanitizers (ASan+UBSan oracle+core clean) + 3-profile `-Werror` builds

Passing normal tests does not imply correctness over unsupported domains
(Batch-1 no oracle; `long double==double` collapses LD paths; large-x trig
crossings are absolute-error regime). Passing the oracle grid does not prove
<0.5 ULP for all inputs (Table Maker's Dilemma; see `ULP_PUSH.md`).

Evidence: `docs/GATE_V12R2.md`. Logs live in `/tmp` by policy (never in repo).

Round-2 (2026-10-01): oracle 212/0 ASan clean (zeta OOB fixed), god 61473/0, boundary 25/0, edge 23/23 sanitizer-clean sample, 3-profile -Werror 34/34. Probes in /tmp/opencode/mathlib-work (acosh/pi/ninv/sp/kelvin/hk) pin each fix.
