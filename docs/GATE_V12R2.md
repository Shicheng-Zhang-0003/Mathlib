# Gate v12R2 / V1.2-RC2 — Evidence (2026-09-30, ULP-push tip)

Commit: ULP-push tip (this tag).
Logs: `/tmp/opencode/mathlib-work/gate-v12R2/` (local) + `/tmp/opencode/mathlib-work/oracle_ulp1` (0-ULP scan); summary below.

## Builds (`-Werror` clean, `/tmp` objects only, repo tree untouched)

- `SCIENTIFIC / GRAPHICS / EMBEDDED`: all 34 TUs compile-only clean.
- `EMBEDDED ml_rsqrt`: link + run OK (was compile failure before despot fix).
- ASan+UBSan CMake configure + build: clean.

## Tests

- modular: `test_core / test_trig / test_linalg / test_dsp` — all passed.
- smoke `test`: 30 passed, 0 failed.
- oracle: **212 passed, 0 failed, worst 0 ULP** (was 5 at push start; n1=0 —
  all 212 correctly rounded on the grid vs mpmath 50-dps exact-binary truth;
  gate stays ≤5 ULP).
- Wide grids: pow 100 cases 0 ULP (was 41); sin/cos benign 25 cases 0 ULP
  (was 1); gamma 11-point grid 0 ULP (6.7 was 21, 0.1 was 3, 0.001 was 5,
  -0.5 was 1).
- ASan oracle: 212/0. ASan core: pass.
- god-mode fuzz `--seed 123456789`: **61473 passed, 0 failed** (Round-2 re-measured; was 61393 at push start).
- boundary gauntlet: **25 passed, 0 failed**.
- edge: **24/24 PASS** (`accuracy_audit, audit_ip1/2/3, combinatorics,
  complex, core, despot_r4, fixed, hyperbolic, integral, linalg, numerical, pow,
  quaternion, redteam_p0_1-5, redteam_p1_1, redteam_p2_logexp, stats, trig`).
  R5: was 23/23 pre-R4; despot_r4 added the 24th suite.
- `despot_check` (17 contract assertions, historical `/tmp` probe): ALL PASS
  at v12R2 time; R5 pins the same contracts in-repo via
  `tests/test_edge_despot_r4.c` (21) + `tests/test_edge_accuracy_audit.c`
  (371). The standalone 17-assertion script was never vendored; do not cite
  it as a current gate tier (see PRECISION_CONTRACT rule 5).

## Profiles / sanitizers

- `SCIENTIFIC`, `GRAPHICS`, `EMBEDDED` compile clean.
- ASan+UBSan oracle + core clean (Round-2: zeta OOB fixed, oracle ASan abort → 212/0 clean).

## Naming

- `v12R2` (internal) == `V1.2-RC2` (public), same commit, same bits.
- `MATHLIB_VERSION_STRING "12.2.0 (v12R2 refinement / V1.2-RC2 public)"`.
