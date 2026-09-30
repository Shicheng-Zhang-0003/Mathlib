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

- Oracle: 212 passed, 0 failed (worst 5 ULP gamma 1e-3). Unchanged.
- Smoke/modular/linalg/dsp: all passed. Fuzz god 61393/0, boundary 25/0.
- `despot_check` (17 assertions in `/tmp/opencode/mathlib-work/`): ALL PASS.
- Full per-TU `-Werror` compile: 35/35 clean. C++ header check clean.
- Edge full sweep not completed in-session (23 suites x 34 TUs > 3 min);
  targeted edge + oracle + fuzz cover the changed paths.

## Remaining (honest, not hidden)

- Batch-1 still no oracle ULP certification (bisection inverses, Krylov
  without preconditioning, xorshift RNG, Dirichlet-only PDEs, 2x2 control).
- `long double==double` platforms (MSVC/ARM) collapse the 80-bit accuracy
  arguments; documented, not fixed.
- Large stacks (harmonic 64KB, GMRES 25KB, PDE 16KB) thread-safe but heavy
  for embedded; heap-workspace API deferred.
- Absolute `1e-12/1e-9` geometric tolerances (stewart/ceva) remain
  arbitrary-but-documented.
