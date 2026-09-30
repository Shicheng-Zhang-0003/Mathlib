# Gate v12R2 / V1.2-RC2 — Evidence (2026-09-30)

Commit: `e649976` + release fixes (this tag).
Logs: `/tmp/opencode/mathlib-work/gate-v12R2/` (local); summary below.

## Builds (`-Werror` clean, `/tmp` objects only, repo tree untouched)

- `SCIENTIFIC / GRAPHICS / EMBEDDED`: all 35 TUs compile-only clean.
- `EMBEDDED ml_rsqrt`: link + run OK (was compile failure before despot fix).
- ASan+UBSan CMake configure + build: clean.

## Tests

- modular: `test_core / test_trig / test_linalg / test_dsp` — all passed.
- smoke `test`: 30 passed, 0 failed.
- oracle: **212 passed, 0 failed**, worst 5 ULP `gamma(1e-3)`.
- ASan oracle: 212/0. ASan core: pass.
- god-mode fuzz `--seed 123456789`: **61393 passed, 0 failed**.
- boundary gauntlet: **25 passed, 0 failed**.
- edge: **23/23 PASS** (`accuracy_audit, audit_ip1/2/3, combinatorics,
  complex, core, fixed, hyperbolic, integral, linalg, numerical, pow,
  quaternion, redteam_p0_1-5, redteam_p1_1, redteam_p2_logexp, stats, trig`).
- `despot_check` (17 contract assertions): ALL PASS.

## Profiles / sanitizers

- `SCIENTIFIC`, `GRAPHICS`, `EMBEDDED` compile clean.
- ASan+UBSan oracle + core clean.

## Naming

- `v12R2` (internal) == `V1.2-RC2` (public), same commit, same bits.
- `MATHLIB_VERSION_STRING "12.2.0 (v12R2 refinement / V1.2-RC2 public)"`.
