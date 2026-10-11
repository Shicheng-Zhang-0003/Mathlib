# Despot Round-5 Audit (2026-10-10, Edmonton) — full math / programming / operational

Method: hostile re-audit of the entire tree against external truth
(mpmath 80-dps, DLMF identities, closed-form physics), strict-C99
program analysis (`-Werror -Wconversion -Wshadow -Wpedantic -Wundef`,
`-Wstack-usage`, ASan+UBSan), and operational cross-check of all four
canonical TU lists + every claimed count. All builds/tests ran in
`/tmp/opencode/despot-r5` (repo tree untouched except the fix itself).
Pinned by `tests/test_edge_despot_r5.c` (42/42).

## 1. Mathematical verdict: NO wrong-answer lie found (verified, not assumed)

Sampled against mpmath 80-dps on this host (gcc 13.3, x86-64):

- Core: sin/cos/exp/log/pow/gamma/lgamma 0 diff at printed precision;
  oracle 212/0 worst 0 ULP holds after R5 (gate stays ≤5 ULP).
- Specials sweep: K0/K1 4..16, Ai -14..16, ber/bei/ker/kei 0.5..18,
  zeta(0/-2/-4/0.5/1.5), digamma(-0.5), gamma(-0.5/-1.1), large-arg
  sin(1e10/1e100), pow(-2,0.5)=NaN — all match mpmath (rel 0 to 1.5e-16).
- Physics outcomes: ODE RK4 exp(1) err 2.2e-10 (dt=0.01, 4th order);
  SUVAT, projectile v0²/g, SHM cos, chi2_inv(0.5,2)=2ln2 (3.8e-15),
  normal_logcdf, entropy ln2, BB mean, CG 2x2, Kalman, sphere norm=1 —
  all exact/plausible. DP5_sys exp(1) 1e-8 pinned in R5.
- Batch-1 (no oracle, openly): CG/GMRES relative tol, SVD/Jacobi
  scale-invariant stops, SDE forward-only dt, BE2 fail-loud, Kelvin DLMF
  full sums, K/Airy quadrature bridges — re-verified, no change needed.

One dead-clause fix (not a wrong answer on valid inputs):
- `ml_hurwitz_zeta` had a dead `s==round(s)&&a==round(a)&&a<=0` limb
  (a<=0 already returned above) + meaningless integer test for
  |v|>=2^53. R5 replaces with explicit |v|>=2^53 early-out + documents
  s<=0,a!=1 NaN stub (needs functional equation). Pinned R5 B5.

## 2. Programming fixes (28 findings, all fixed)

Critical:
- `cpu_dispatch.c`: non-GNUC path called `ml_matmul_avx2()` unguarded
  (compile fail / SIGILL risk); now scalar. `ml_cpu_has_avx2()` returned
  build-flag constant (lies on heterogeneous CPUs); now 0 when undetectable.
- `fixed_point.c`: `1LL<<shift` UB (shift>=63/<0), `-value` UB
  (INT64_MIN) in both asr helpers; now fail-safe (saturate + exact
  -2^(63-shift) for MIN). In-tree shifts 0..16 unaffected.
- LDBL guards: 12× `(LDBL_MANT_DIG>=64)` without `defined()` (silent
  fallback if float.h ordering changes; caused complex-log incident).
  Unified to `defined(__LDBL_MANT_DIG__)&&(__LDBL_MANT_DIG__>=64)`
  (trig.c/minimax.h/payne_hanek.h pattern). `ml_compiler.h` documents it.
- `ml_cross3`: out==a alias read-after-write (wrong y,z); now temporaries.
- `ml_transpose`: in-place square garbaged upper triangle; now correct
  swap path (out-of-place remains contract).
- `ml_egcd`: recursive `(a/b)*y1` signed-overflow UB; now iterative with
  `__builtin_mul/sub_overflow` fail-loud (pedantic-clean, no __int128).
  Sign fixed for negative a (was x=1 for a=-5,b=0; now x=-1).
- `fft.c`: `1<<24` signed edge → plain `16777216` (1U broke -Wsign-compare).
- `ml_compiler.h`: MSVC `__builtin_*` mapping layer added (sqrt/log/cos/
  exp/sqrtl/cosl/sinl/expl/logl/powl/fabsl → libc); precision then follows
  MSVC LD==double (documented; K/Airy exempt).

Robustness/API:
- Status precision: CG/GMRES/power/QR-eig exhaustion, ode_dp5/be2
  non-convergence, calculus Gauss-Legendre residual, control LQR it==99,
  linalg Jacobi off>tol now `ML_ERR_CONVERGENCE` (were SINGULAR/INTERNAL);
  LU overflow-Inf now `ML_ERR_OVERFLOW`; circular_conv overflow now
  OVERFLOW. True singularity keeps SINGULAR (ABI; only test `zero matrix
  is singular` depends on it — preserved).
- DoS caps: `ml_majorizes` n>1M→0 + SIZE_MAX guard; sieve SIZE_MAX guard
  (OOM sentinel UINT64_MAX documented); `ml_ess` n>100k→NaN (was O(n²)
  uncapped); DCT-II/III/DST n>4096 silent-return (were unbounded O(n²));
  parseval n>1M→NaN; Pollard-rho 1M-iteration cap (was unbounded Floyd).
- RNG: SDE euler/milstein ±ctx + MCMC MH now share the 16-round xorshift
  warmup (bridge_sample already did; same seed gave differently-biased
  first Gaussian before). Seed-0→fallback documented.
- ODE: dp5/be2 `t==t1` exact-float arrival → `>=` clamp + snap (was
  spin-until-h-underflows on some steps).
- Thread/stack: `ml_workspace_t` one-per-thread contract in
  `ml_tensor.h`; existing n-caps already bound all large frames
  (pde/calculus/harmonic/linalg/spectral/info/manifold) — documented,
  no silent stack blowup. `ml_tensor_at` checked vs `ML_TENSOR_AT`
  unchecked documented; `ml_matmul` no-alias documented.

## 3. Operational fixes (20 findings, all fixed)

- Counts (honest, host-noted): god-mode is **host-sensitive**
  (`rand()` + conditional CHECKs + AVX paths; +4 AVX-only, ±80 across
  libc). This host (AVX2, gcc 13.3): **61473/0** (`--seed 123456789`);
  other hosts measured 61393/0. Gate is **0 failed**, not an exact count
  (old docs gated the count and called 61473 "never reproducible" —
  reproduced here three times). Updated: README, TESTING, KNOWN_LIMITS,
  PRECISION_CONTRACT R5 notes; GATE_V12R2 left historical + R5 addendum.
- Edge: 25/25 (24 + new `test_edge_despot_r5` 42 pins); was 23/23→24/24
  stale in GATE/release_notes. `accuracy_audit` 371 (was 282 stale in
  KNOWN_LIMITATIONS:59 — fixed). `despot_check` 17 was a `/tmp`-only probe
  never vendored; docs now mark it historical, pinned in-repo by R4+R5.
- `closure_gate.sh`: was 7 steps, missing smoke/god/python-runner,
  `rm -rf build`, hardcoded build, v12A1 banner. Now 10 steps
  (modular+smoke+edge+boundary+god+oracle+ultimate+python-runner),
  honors `MATHLIB_BUILD_DIR`, opt-in clean, v12R2 banner, prefers CMake
  ultimate_fuzzer target.
- `run_all_tests.py`: v12A1 docstring→repo-root/v12R2; honors
  `MATHLIB_BUILD_DIR`, `MATHLIB_EDGE_BUILD_DIR`,
  `MATHLIB_EDGE_SANITIZERS`; banners v12R2.
- `CMakeLists.txt`: new `ultimate_fuzzer` target (was divergent ad-hoc
  gcc lines). `.gitignore`: +`/ultimate_fuzzer` +`/oracle_check`,
  `/test_oracle` kept as legacy alias.
- Root clutter: `bench_graphics`/`bench_scientific` were tracked ELFs
  despite ignore → `git rm --cached`; untracked root `test_*` removed.
  `build/` stays ignored; evidence goes to `/tmp`.
- Soak: `soak_test.sh` + `tests/run_soak.sh` + `run_fuzz_vault.py` wrote
  `tests/regression/artifacts/` in-tree → now
  `${MATHLIB_SOAK_DIR:-/tmp/opencode/mathlib-soak}` +
  `${MATHLIB_BUILD_DIR:-build}`.
- Scripts: `scripts/closure/README` + `scripts/v12a1/README` referenced
  missing bootstrap scripts / stale v12A1 paths → marked HISTORICAL with
  surviving-sequence notes. TU lists re-verified 34/34 in all four
  canonical files.

## 4. Gate on the R5 tree (measured `/tmp/opencode/despot-r5`)

- oracle 212/0 worst 0 ULP (gate ≤5 ULP), modular 4/4, smoke 30/0,
  boundary 25/0, god 61473/0 (`--seed 123456789`, AVX2 host; 0 failed
  is the gate, count is host-sensitive), edge 25/25 (incl. R4 21/21 +
  R5 42/42), accuracy_audit 371/0, 3-profile `-Werror` clean,
  ASan/UBSan clean (edge sanitized run + CMake sanitizer build).

## 5. What remains (documented limits, not defects)

- No universal correct rounding; no cross-arch bit identity; Batch-1
  still no oracle ULP certification (regression-guard tier only).
- `long double==double` (MSVC/ARM): LD paths collapse (K/Airy exempt).
- `ml_hurwitz_zeta` s<=0,a!=1 NaN stub (needs functional equation).
- Y0/Y1 absolute 1.1e-13 at zeros (intrinsic); Airy |x|>1e6 phase loss;
  lgamma near-zero relative caution. See KNOWN_LIMITATIONS.
