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

## Round-3 — 2026-10-01 (close the remaining limits)

Every limit still open after Round-2 was attacked with real method changes,
not tolerance tweaks. Ground truth: mpmath 80-dps. All work in
`/tmp/opencode/mathlib-work/round3/`.

### Kelvin — 2-term truncation CLOSED (DLMF 10.67.3-4 full sums)

Old `ml_kelvin_asym` used 2-term P/Q polynomials for all four functions
(~0.5% at x=20, `bei(20.1)` 125562 vs 126161 true). DLMF 10.67.3-4 is a full
sum over `a_k/x^k` with a rotation that differs per function:

    a_k(0) = prod_{j<=k} (-(2j-1)^2) / (k! 8^k),  a_0 = 1
    ber ~ sqrt(1/(2 pi x)) e^{x/sqrt2} sum a_k x^-k cos(x/sqrt2 - pi/8 + 3k pi/4)
    bei ~ sqrt(1/(2 pi x)) e^{x/sqrt2} sum a_k x^-k sin(x/sqrt2 - pi/8 + 3k pi/4)
    ker ~ sqrt(pi/(2x))  e^{-x/sqrt2} sum a_k x^-k cos(x/sqrt2 + pi/8 + k pi/4)
    kei ~ sqrt(pi/(2x))  e^{-x/sqrt2} sum a_k x^-k sin(x/sqrt2 + pi/8 + k pi/4)

Implemented with the `a_k` recurrence, least-term truncation, and LD
accumulation. Measured vs mpmath 80-dps: `ber/bei/ker/kei` at 20, 20.1, 30,
50 now agree to **~1e-13 relative** (was 0.5% / sign-flipped).

`ker/kei` crossover split from `ber/bei`: `ML_KELVIN_KXMAX = 12` (was the
shared 20). Measured: series-to-12 beats series-to-20 because the series
loses digits to cancellation earlier than the asymptotic does.

### Bessel K0/K1 — 8-10 transition hole CLOSED (K-quadrature bridge)

The 8<x<10 hole was structural: the ascending series cancels like
`(ln(x/2)+gamma)I0` (I0 ~ 1e3 against K ~ 6e-5 at x=8.8) while the
asymptotic bottoms out at its least term, `exp(-2x)`. Neither can do better
in double precision; that is what Temme's uniform expansion exists for.

Round-3 replaces the hybrid with a **third branch**: the defining integral

    K_nu(x) = int_0^inf exp(-x cosh t) cosh(nu t) dt   (DLMF 10.32.10)

evaluated by long-double adaptive Simpson with Richardson correction. The
integrand is positive and unimodal, so there is no cancellation at all —
this is the Temme-class accuracy without the uniform-expansion machinery.

Two details mattered:
- **Relative tolerance.** A fixed `1e-18` absolute only buys ~1e-11 relative
  once K falls to `exp(-15)`. Tolerance is now `1e-19 * exp(-x)`.
- **Adaptive upper limit.** A fixed `T=10` spent ~90% of the work on a tail
  that was already 1e-3000. `T = acosh(1 + 48.4/x)` cuts the cost 10x.

Crossovers re-tuned by measurement: series `x<4`, quadrature `4<=x<16`,
least-term asymptotic `x>=16`.

Measured vs mpmath 80-dps over x in [1e-3, 300]: **worst 6.8e-15 relative**
(was ~2e-9). 31 of 32 pinned grid points are now **0 ULP**.

### Airy Ai — 5-6 transition hole CLOSED (K_1/3 bridge)

Same root cause. DLMF 9.11.4 gives `Ai(x) = sqrt(x/3)/pi * K_1/3(zeta)`
with `zeta = (2/3) x^(3/2)`, so the *existing* K-quadrature bridge serves
Airy directly with no new code. Branch map: LD Taylor for `x <= 2.5`,
`K_1/3` quadrature for `2.5 < x < 8.5`, least-term asymptotic above.

Measured vs mpmath 80-dps on a 0.1-step sweep over `[-14, 200]`:
**worst 2.1e-13 relative** (was ~5e-9 documented, ~3e-9 after Round-2's LD
Taylor). Negative-side `ML_AIRY_XN` re-tuned 7 -> 8 by sweep (9 and 12 blow
up to 2e-6 and 154 because the negative Taylor series diverges).

### Y0/Y1 and digamma

- `ML_BESSEL_XY` 14 -> 13: `Y1(13.9)` 2.56e-12 -> 1.08e-13. Swept 11.5/12/
  12.5/12.75/13/14 — 12.5 and 12.75 are indistinguishable from 13 because
  the worst point (11.75) is inside the series branch, not at the seam.
- `ml_digamma` positive branch rewritten in long double (recurrence +
  Stirling, single round). `digamma(1.5)` was 10 ULP from double
  cancellation in `3.48 - 3.44`; now 0 ULP.
- `ml_digamma` negative branch: exact half-integer shortcut (cot(pi x) is
  identically 0 there and `ml_cospi(-0.5)` leaves a spurious ~2e-16), plus an
  LD `cot(pi x)` for the general case. `digamma(-0.5)` now exact.

**Honest characterisation of what is left in Y0/Y1:** worst *absolute*
error is 1.1e-13 (Y0 @13, Y1 @12.75). The 1.15e-10 *relative* figure at
x=11.75 is not a defect — `Y1(11.75) = -1.96e-4` sits on a zero, so
`1e-13/1.96e-4 = 1e-10`. Term analysis in mpmath confirms the ascending
series has `sum|term|/|S| = 5e4` (4.7 digits of cancellation) and the
long-double term-representation floor is 5e-15 relative to S. This is the
same "absolute error at crossings" regime already documented for large-x
trig, and no elementary method does better near a zero.

### `long double == double` collapse — CLOSED for K and Ai

Previously documented as unfixable for every LD path. Measured directly by
compiling the quadrature in plain double (an MSVC/ARM emulation): **K worst
2.6e-16, Ai worst 2.4e-15** — no degradation at all. Reason: the quadrature
integrand is positive and non-cancelling, so double precision is already
enough; the series/asymptotic hybrids it replaced were the parts that needed
80-bit mantissas. K0/K1/Ai are therefore now *more* portable than before.

### Stack diet

`-fstack-usage` survey of all 34 TUs. `ml_fft2d_pow2` held
`cplx rows[128][128]` = **256 KB** — the largest frame in the tree, hostile
to embedded. The matrix now lives in the caller's `re_out`/`im_out` and only
one `cplx cols[128]` (2 KB) scratch is needed. Verified against a direct 2D
DFT: max deviation 2.49e-15. New tree maximum is `ml_fft_real` at 65,600 B.

### Regression guard tightened

`tests/test_edge_accuracy_audit.c`: K tolerances 5e7 -> 65536 ULP and Airy
1e8 -> 65536 ULP (~1.5e-11 relative, 700x tighter; measured worst 0-8 ULP,
with headroom for `long double == double`). Added 10 new transition-band
pins (K at 5.5/7/13, Ai at 4/4.5/6.5/8). Suite now **371 assertions**.

### Still open (honest, not hidden)

- `ml_hurwitz_zeta(s<=1, a!=1)` — analytic continuation stub, returns NaN.
  Deliberate fail-loud, unchanged.
- Batch-1 oracle coverage (bisection inverses, unpreconditioned Krylov,
  xorshift RNG, Dirichlet-only PDEs, 2x2 control). Structural: needs an
  oracle tier that does not exist for these families yet.
- Y0/Y1 absolute error ~1e-13 near zeros — intrinsic, characterised above.
- `long double == double` collapse still applies to the *core* LD paths
  (gamma/pow/log/sin), unchanged; K and Ai are now exempt.
- Heap-workspace API for the remaining 64 KB / 33 KB frames (fft_real,
  jacobi_eigen, mi_discrete, haar, cubic_spline) — deferred, documented.
