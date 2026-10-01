# Despot Audit — 2026-09-30

Full mathematical, programming, and operational audit. Every claim below is
implemented in code and pinned by `despot_check` or the existing gauntlet.

## Fixed (code + math)

- `ml_exp2`: cutoff `-1074.0 -> -1075.0` (`exp_log.c:664`). `2^-1074` is min
  subnormal, not zero. Dead `(void)n` removed.
- `ml_remainder`: quotient in `long double` (`core.c:250`). Double-rounded
  `x/y` within ~1 ULP of a half tie no longer picks the wrong neighbor.
- `ml_erfinv`: Acklam starter documented; `+-1 -> NaN` recorded as fail-loud
  contract (limit `+-Inf`); dead `w` + stale Winitzki comment removed.
- `ml_coth`: dead `(void)e` removed. `cplx_sqrt`: Annex-G `(-0,+0)->(+0,+0)`;
  dead `(void)m` removed. `ieee754.h`: marked DEPRECATED (naive series,
  subnormal flush; nothing includes it).
- `ml_kronecker`: `INT64_MIN` UB fixed via `ml_abs_u64` magnitude domain
  (`numbertheory.c:63`). Unsigned loop throughout.
- `ml_mobius`: `2` documented as UNRESOLVED sentinel (header + code).
- `ml_mult_order`: gcd-guarded residual reduction; composite residual cannot
  over-divide. Returns the exact order, not a multiple.
- `ml_crt2`: `%m2 -> %q` fix (`numbertheory.c:196`). `0 mod6, 4 mod10`
  now `24 mod30` (was `6`). `p%q` inversion, `INT64_MAX` egcd guard.
- `ml_covariance`: `n<=1 -> NaN` (was `0.0` for `n==1`), consistent with
  `ml_variance_s`.
- `ml_cross_entropy`: `H+Inf -> +Inf` (was `NaN`); NaN only for bad input.
  `mi_discrete` dead pointer removed.
- `ml_mh_sample_ctx`: non-finite proposal is a REJECTION with accumulation
  (was `continue` shrinking the sample); `log(u+1e-300) -> log(u)` with
  `u==0 -> -Inf` (never accept, correct MH).
- `ml_sde_*`: backward `dt<0` rejected (Wiener `sqrt(|dt|)` vs drift sign
  inconsistency); `sqrt(dt)` explicit; NEW `_ctx` variants thread user state
  (legacy passes NULL). Bridge mean vs sampler documented; OU exact kept.
- `ml_ode_dp5_sys`: `k1` finite check before `yt` formation.
- `ml_ode_be2_sys`: Newton failure forces step rejection (`h*=0.5`, retry),
  never acceptance of bogus `yf/ym`. Residual re-check added.
- `ml_cg/gmres`: relative `tol*||b||` (was absolute). GMRES finite pre-check
  added. `ml_svd_jacobi`: scale-invariant `(eps*||A||_F)^2` stop (was
  absolute `1e-30`); dead `UU` removed. `ml_qr_iter_eig` + GMRES Givens via
  `ml_hypot` (overflow-safe at 1e200); dead loop left intact as documentation.
- `ml_jacobi_eigen`: SINGULAR on `off>tol` exhaustion (was silent SUCCESS).
  `ml_determinant`: finite-break on overflow product.
- `ml_kelvin ber/bei(Inf)`: `0.0 -> NaN` (unbounded oscillatory). `XMAX`
  macro used (was literal `20.0` drift).
- `ML_FIXED_HALF_PI/TWO_PI`: `102943/411774 -> 102944/411775`
  (round, not trunc; 1-LSB bias removed). `1<<shift -> 1LL<<shift`.
- `ml_dst2`: documented as DST-I with `ml_dst1` alias; DCT unscaled
  normalization documented; `parseval_energy` naming clarified.
- `ml_haar_fwt`: buffer contract `det[n-1]` documented; dead `tmp` self-copy
  removed; `ml_haar_iwt` validates finiteness (NaN propagation, not garbage).
- `ml_exp_sphere`: `|x|=1` enforced + scale-aware tangency. `ml_sphere_dist`:
  `2*asin(|x-y|/2)` small-angle stable + unit-input gate with acos fallback.
- `ml_lyapunov_2x2_trace`: documented as Hurwitz margin, not Lyapunov trace.
- `ml_polynomial_derivative degree==0`: vacuous return documented.
  `ml_carlson_rf`: dead `e2` removed. `ml_cubic` return `0/1/2/3` documented.
- `profiles.h` EMBEDDED: `ml_rsqrt` fallback added (was missing -> compile
  failure). `ml_types.h`: `OVERFLOW/CONVERGENCE/UNSUPPORTED` added
  (ABI-additive). All 45 public headers: `extern "C"` guards (C++ safe,
  verified `g++ -std=c++17`).
- Build coherence: kelvin in CMake + Makefile + run_all_tests.py +
  run_edge_tests.sh (34 TUs everywhere; was CMake-only drift).

## Validation

- Oracle: 212 passed, 0 failed, worst 0 ULP on the grid (was 5 at push
  start; gate stays ≤5 ULP — measurement, not proof, see `docs/ULP_PUSH.md`).
- Smoke/modular/linalg/dsp: all passed. Fuzz god 61393/0, boundary 25/0.
- `despot_check` (17 assertions in `/tmp/opencode/mathlib-work/`): ALL PASS.
- Full per-TU `-Werror` compile: 34/34 clean × 3 profiles. C++ header check clean.
- Edge full sweep: 23/23 PASS (gate-v12R2 logs).

## Remaining (honest, not hidden)

- Batch-1 still no oracle ULP certification (bisection inverses, Krylov
  without preconditioning, xorshift RNG, Dirichlet-only PDEs, 2x2 control).
- `long double==double` platforms (MSVC/ARM) collapse the 80-bit accuracy
  arguments; documented, not fixed.
- Large stacks (harmonic 64KB, GMRES 25KB, PDE 16KB) thread-safe but heavy
  for embedded; heap-workspace API deferred.
- Absolute `1e-12/1e-9` geometric tolerances (stewart/ceva) remain
  arbitrary-but-documented.

## Round-2 — 2026-10-01 (despot returns)

Full mathematical/programming/operational re-audit in `/tmp/opencode/mathlib-work/`.
Probes built with `gcc -std=c99 -O2 -fno-fast-math -ffp-contract=off` vs
mpmath 80-dps / system libm. All fixes below are implemented and re-verified
(ASan+UBSan oracle 212/0, 3-profile -Werror clean).

### P0 — wrong answer / sanitizer abort (fixed)

- `integral.c:306-332 zeta[]` — missing ζ(12) shifted k>=12 by one slot and
  `zeta[25]` was OOB (`long double[25]` with kmax=25). ASan
  `global-buffer-overflow` + UBSan `index 25 OOB` on `oracle_check`.
  Fixed: 26-entry table `[0..25]` with mpmath 80-dps constants, ζ(12)=
  `1.000246086553308048298637998047739670960L` restored, full precision
  for k>=12 (was 15-digit truncation). Masked by `u^k/k` (`|u|<=0.15`)
  but UB nonetheless — now 0 ULP on lgamma zeros grid.
- `kelvin.c:79 asym phase` — single `+pi/8` for all four functions.
  DLMF 10.67.3-4: ber/bei use `-pi/8`, ker/kei use `+pi/8`.
  Was sign flip: `ber(20)=-47186` vs true `+47489`, `ber(20.1)` rel -2.4.
  Fixed: split `phb/phk`. Residual 2-term truncation ~0.5% at 20
  (EXPERIMENTAL, documented, no oracle).
- `numbertheory.c:231 prime_pi` — `size=nn/2` missed odd `nn` itself.
  `pi(3)=1` (true 2), `pi(5)=2` (3), `pi(7)=3` (4). Fixed: `(nn+1)/2`.

### P1 — systematic bias (fixed)

- `exp_log.c:616 acosh` — `log(x+sqrt((x-1)(x+1)))` cancels at 1+
  (458k ULP at `1+1e-12`). Fixed: `log1p((x-1)+sqrt(...))` → 0 ULP vs sys.
- `info.c:75 softplus` — `x>20 return x` drops `log1p(exp(-x))`
  (57k ULP at 22). Fixed: cutoff 36 (`exp(-36)~2.3e-16`), `x+log1p(exp(-x))`
  for 20<x<=36 → exact.
- `statistics.c:310 normal_inv` — Newton `e/(sigma*pdf)` scales step by
  `1/sigma` (dCDF/dz is `phi`, not `sigma*phi`). Stalled 10 iters for
  sigma!=1 (`inv(.99,5,.001)=5.248` vs `5.002`). Fixed: `e/pdf` → exact,
  ratio `inv(.975,0,2)/inv(.975,0,1)=2.0`.
- `analytic_nt.c:14 hurwitz` — tail was integral+1/2 only, missing
  `s/(12 N^{s+1})` (~10k ULP at s=1.1: `-2.28e-05`). Fixed: B2 term added
  → `9.5844484658669` vs `9.5844484649508` true (9e-10).
- `combinatorics.c:310 catalan` — cap `n>34` falsely overflowed C35/C36
  (both <2^64) and `ml_ncr(2n,n)` overflowed intermediate. Fixed:
  recurrence `C_{k+1}=C_k*2(2k+1)/(k+2)` with GCD cancellation (strict C99,
  no `__int128` for `-Wpedantic`); C0..C36 exact, C37+ MAX.
- `exp_log.c:710 log1p lc[14]` — `0.069` truncated (true `2/29=
  0.06896551724137931`). Masked (`z^29~1e-21`) but fixed for hex-exact intent.
- `integral.c:803-816 Airy/K` — Taylor was double (24M ULP at 5.5),
  now LD Kahan (`1e-22` stop) → ~2.7e-10 at 5.0. XK=9 kept (series 5e-10
  at 8.8 beats asym 2e-9; lowering to 7 worsens). Residual 5-6 hole and
  8<x<10 K hole need Temme uniform expansion — documented, not hidden.

### Programming hardening (fixed)

- `spectral.c:226 qr_iter_eig` — SINGULAR left `evals_*` untouched (stale
  garbage on ignored status). Fixed: poison to NaN on both exhaustion and
  post-check `off>tol` paths; post-check added (was silent SUCCESS on
  unconverged real 2x2 block).
- `control.c:99` — added honestly-named `ml_hurwitz_margin_2x2` alias;
  legacy `ml_lyapunov_2x2_trace` kept for ABI, doc points to alias.
- Headers: 7 use `LIBMATHC_` prefix vs 38 `MATHLIB_` — both guarded,
  no missing guard; standardized in docs (no ABI break).
- `__int128` rejected: strict C99 `-Wpedantic` forbids it (build break);
  recurrence above is the portable fix. `malloc` in `prime_pi`/`majorizes`
  is heap-per-call (thread-safe), not core zero-alloc violation — documented.
- Thread-safety re-verified: `grep static.*\[` shows only `static const`
  tables + stateless helpers; no per-call mutable state.

### Operational (fixed)

- `GATE_V12R2.md:8`, `BUILD_AND_INSTALL.md:10`: `35 TUs` → `34 TUs`
  (canonical list is 34 everywhere: CMake/Makefile/run_all_tests/edge).
- `run_all_tests.py:586` `Path("src/core.c")` false-positive in naive
  `grep -c src/` counts (35/36) — canonical `LIB_SOURCES[70:103]` is 34.
- Validation Round-2: oracle 212/0 worst 0 ULP ASan clean (was abort),
  modular 4/4, smoke 30/0, boundary 25/0, god 61473/0 seed 123456789,
  edge 23/23 (sanitizer sample clean), 3-profile -Werror clean.

## Remaining (honest, not hidden) — updated

- Kelvin 2-term truncation ~0.5% at 20; crossover 20 kept (series holds to
  19 to 6 digits, asym fixed). Uniform `I0/K0` continuation deferred.
- Airy 5-6 hole (~1e-9..1e-8) + K 8-10 (~1e-9) need Temme/Amos; LD Taylor
  halves the hole but does not close it. Pinned by `accuracy_audit`.
- Batch-1 still no oracle ULP (bisection inverses, Krylov unpreconditioned,
  xorshift RNG, Dirichlet PDEs, 2x2 control) — see PRECISION_CONTRACT.
- `long double==double` (MSVC/ARM) collapses LD paths; guarantees stay ≤5 ULP.
- Large stacks (harmonic 64KB, GMRES 25KB) thread-safe but heavy; heap API deferred.
