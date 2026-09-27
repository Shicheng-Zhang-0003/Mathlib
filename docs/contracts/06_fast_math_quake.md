# Math Diary: Fast Math (Quake III Inverse Sqrt & Integer-Float Isomorphism)

> **STATUS: OWNED (2026-09-27)** — verified against source. Fast paths are
> approximate-by-contract and have no oracle coverage per
> `docs/PRECISION_CONTRACT.md`. v11S closure: domain-guarded, UB-hardened,
> not correctly-rounded libm replacements.

> **Orchestrator's Mission (done):** reverse-engineered from
> `include/mathlib/fast_math.h`, `simd_batch.h`, `simd_bare_metal.h`,
> `src/cpu_dispatch.c`, `include/mathlib/profiles.h`,
> `validation/fma_audit.s`, `tests/fuzz_god_mode.c`, `benchmarks/bench.c`.

---

## 1. The Mathematical Identity
*Approximate-by-contract — never oracle tier.*

- **Core identity:** `y ≈ 1/sqrt(x)`, i.e. `y` solves `1/y² - x = 0` via Newton iteration
  `y ← y * (1.5 - x*0.5*y*y)`.
- **Scalar kernel `ml_fast_rsqrt` (`include/mathlib/fast_math.h:32-53`):**
  integer initial guess `bits = 0x5fe6ec85e7de30daULL - (bits >> 1)` then
  **two** Newton refinements (`y = y*(1.5-(number*0.5*y*y))` twice, lines 51-52).
- **`fast_log2` (`ml_fast_log2`, `fast_math.h:68-90`):** frexp isomorphism
  `x = m*2^e`, `m ∈ [0.5,1.0)`, `f = 2m-1 ∈ [0,1)`,
  `log2(x) = (e-1) + log2(1+f)` with 3rd-degree poly
  `p = f*(1.442695 + f*(-0.721347 + f*0.278652))`.
- **`fast_exp2` (`ml_fast_exp2`, `fast_math.h:105-152`):** split `x = xi + frac`
  (floor for negatives), 5th-degree poly for `2^frac` with `ML_LN2`,
  `0.2402265069591007, 0.05550410866482158, 0.009618129107628477,`
  `0.00133335581464249`, `mant_approx = 1+p`, exponent inserted as
  `bits = biased<<52`, subnormal/underflow routed to `ml_ldexp_pure`.
- **Tier status:** per `docs/PRECISION_CONTRACT.md` only the oracle tier
  (`tests/test_oracle.c`, mpmath 80 digits, ≤5 ULP, 7 families, 212 vectors)
  certifies precision. Fast-math paths are smoke/edge tier at best —
  approximate by design, never cited for ULP claims.

## 2. The Algorithm (Execution Path)
*Real names, real files.*

1. **Phase 1: Scalar rsqrt guard + guess — `ml_fast_rsqrt` (`fast_math.h:32-53`)**
   - `ml_isnan → return input; number<0 → ml_make_nan(); ml_isinf → 0.0;`
     `number==0 → ml_copysign(Inf,number)` (so `+0→+Inf`, `-0→-Inf`).
   - Positive subnormal (`ml_is_subnormal`) bails to exact path
     `1.0/ml_sqrt(number)` (line 40) — no bit-hack on denormals.
   - Else `memcpy` bits, magic subtract, `memcpy` back, two Newton steps.
2. **Phase 2: Batched rsqrt — `ml_simd_batch_rsqrt` (`simd_batch.h:23-74`, AVX2;
   scalar fallback lines 86-94)**
   - Exception-lane pre-scan (lines 40-48): if any lane is NaN/Inf/subnormal
     or `!(in[i]>0.0)` (covers ±0, negatives), the **whole batch of 4** falls
     back lane-by-lane to `ml_fast_rsqrt` — batch exceptional behavior equals
     scalar contract exactly (`+0→+Inf`, `-0→-Inf`, neg→NaN, NaN→NaN,
     `+Inf→+0`, subnormal→finite fallback).
   - Clean all-positive-normal batches only: `_mm256_loadu_pd`,
     `_mm256_castpd_si256`, broadcast magic
     `_mm256_set1_epi64x(0x5fe6ec85e7de30daLL)`, `srli 1`, `sub_epi64`,
     cast back, then **two** vector Newton rounds with `half=0.5`,
     `three_half=1.5` (`mul/sub` sequence lines 62-71), `_mm256_storeu_pd`.
3. **Phase 3: Bare-metal rsqrt — `ml_avx2_fast_rsqrt(__m256d v)`
   (`simd_bare_metal.h:10-28`, `__AVX2__` only)**
   - Same integer guess in-register (zero-cost casts, 1-cycle integer
     shift+subtract, no scalar division). Only **one** Newton round
     (lines 24-27) — hence its looser fuzz tolerance (1e-3) vs the 2-step
     batch kernel (1e-5). No exception handling inside; callers must
     pre-screen (the batch wrapper does).
4. **Phase 4: IEEE plumbing for log2/exp2 — `ml_frexp_pure` / `ml_ldexp_pure`
   (`src/core.c:63-91` / `11-61`, via `internal/ieee_exact.h`)**
   - `ml_fast_log2` routes through `ml_frexp_pure` so subnormals are
     normalized before the exponent is read (no raw-exponent-field bug).
   - `ml_fast_exp2` bounds input first (`NaN→NaN`, `x≥1024→+Inf`,
     `x<-1074→0`, avoiding UB long-long cast), then `biased = xi+1023`;
     `biased≤0` routes to `ml_ldexp_pure(mant_approx,xi)` for gradual
     underflow; `biased≥0x7FF → +Inf`.
5. **Phase 5: Vector matmul dispatch — `ml_matmul` / `ml_matmul_avx2` /
   `ml_matmul_scalar` (`src/cpu_dispatch.c:55/26/10`, decl
   `include/mathlib/cpu_dispatch.h:15`)**
   - Naive O(N³) `i,j,k` nest. AVX2 path: per `(i,k)` broadcast
     `_mm256_broadcast_sd(&A[i*n+k])`, stream `_mm256_loadu_pd(&B[k*n+j])`,
     accumulate `_mm256_fmadd_pd`, tail `j` scalar; `size_t` indexing,
     NULL/`N≤0` early return, `n > SIZE_MAX/n` overflow guard, no-alias
     `ML_RESTRICT` contract.
   - Dispatch: `#if ML_COMPILE_TIME_AVX2` + **runtime**
     `__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")`
     (line 69) selects AVX2 else scalar — 2026-09-27 fix: was a
     compile-macro (`__FMA__`) test that reflected build flags, not host
     capability (see `ml_cpu_has_fma` comment lines 92-93). Without it an
     AVX2 binary on a non-AVX2 CPU would SIGILL.

## 3. The "Dark Arts" (Hardware & IEEE-754 Sympathy)
- **Magic constant `0x5fe6ec85e7de30da`:** 64-bit-double variant of the Quake III
  `0x5f3759df` float hack (Chris Lomont / Charles McEniry derivation). It is
  `floor(0x3FE6EC85E7DE30DA-ish)` tuned so `(magic - bits>>1)` reinterpreted
  as double lands within a few percent of `1/sqrt(x)` on the log-exponent
  line; Newton polishes the mantissa. Same literal in `fast_math.h:46`,
  `simd_batch.h:53`, `simd_bare_metal.h:14` (verified by grep).
- **Why 2 Newton steps (scalar + batch) vs 1 (bare-metal):** one step leaves
  ~1e-3 relative error (exactly the `fuzz_god_mode.c:377` bare-metal gate);
  the second step buys ~1e-5–1e-6 (the `fuzz_god_mode.c:367`
  `|expected|*1e-5+1e-12` batch gate and `test_edge_redteam_p1_1.c:97-98`
  comment). Third step costs a full dependent mul-chain for sub-ULP gains
  nobody in the fast tier pays for — GRAPHICS workloads are raster-bound.
- **`memcpy` bit-twiddle vs `union` UB:** header comment `fast_math.h:18-19`
  is explicit — `memcpy(&bits,&number,8)` / `memcpy(&y,&bits,8)` guarantees
  ISO C99 strict-aliasing compliance; `union` type-punning is UB in this
  codebase's reading and trips optimizers at `-O2/-O3`.
- **`vfmadd213sd` proof in `validation/fma_audit.s`:** lines 38/40/43 show the
  Horner loop of `audit_exp` compiled to contracted FMA
  (`vfmadd213sd (%rax),%xmm0,%xmm1` ×2 per unroll + remainder) — evidence the
  toolchain emits single-rounding FMA where the audit expects it, not
  separate mul+add. Scalar `ml_matmul_avx2` likewise uses `_mm256_fmadd_pd`.
- **GRAPHICS profile routing (`include/mathlib/profiles.h:8-18`):**
  `MATHLIB_PROFILE_GRAPHICS → ml_rsqrt(x) = ml_fast_rsqrt(x)`;
  default/SCIENTIFIC `→ 1.0/ml_sqrt(x)`; `MATHLIB_PROFILE_EMBEDDED` hides
  FPU fast math entirely. Inline functions (not macros) per the
  macro-poisoning fix — no double-evaluation.

## 4. Error Budget & Failure Modes
- **Maximum error (fuzz-measured, not ULP-certified):** SIMD batch rsqrt gate
  is `tol = |1/sqrt(x)|*1e-5 + 1e-12` over 1024 random `[1,100]` inputs
  (`fuzz_god_mode.c:358-369`); bare-metal single-step gate is `1e-3`
  (`:376-377`); directed edge gate is `1e-4` relative
  (`test_edge_redteam_p1_1.c:107-115`); core smoke `1e-4`
  (`test_core.c:15`). No MPFR ULP figure exists or is claimed — fast math
  is outside the oracle tier by contract.
- **Catastrophic/exceptional lanes:** never vectorized blindly. Any
  NaN/Inf/±0/negative/subnormal lane forces whole-batch scalar fallback
  (`simd_batch.h:40-48`) to the `ml_fast_rsqrt` contract
  (NaN→NaN, neg finite→NaN incl. `-Inf`, `+Inf→+0`, `+0→+Inf`, `-0→-Inf`,
  subnormal→`1/ml_sqrt`, verified in `test_edge_redteam_p1_1.c:25-89` and
  `test_edge_audit_ip1.c:51-65`).
- **Overflow/underflow boundaries:** `ml_fast_exp2`: `≥1024→+Inf`,
  `<-1074→0`, exponent-insert overflow `biased≥0x7FF→+Inf`, underflow via
  `ml_ldexp_pure` gradual path. `ml_fast_log2`: neg→NaN, `0→-Inf`,
  `+Inf→+Inf`. `ml_fast_rsqrt`: `+Inf→0`, `0→signed Inf` (see §2.1).
- **NaN propagation:** payload-preserving `return x` on NaN in all three fast
  kernels; negative finite rsqrt returns fresh `ml_make_nan()`.
- **Benchmark-only carve-out:** `benchmarks/bench.c:12-26` isolates non-portable
  `rdtsc` (`__asm__ volatile("rdtsc")`) under
  `#if defined(__x86_64__)||defined(__i386__)` with `BENCH_UNITS "cycles"`,
  else `clock()` fallback with `"ticks"`. Excluded from the lib build —
  the only sanctioned inline-asm cycle counter in-tree.
- **Documented gap — no physics kernels vectorized:** grep for `_mm256_|immintrin`
  hits only `src/cpu_dispatch.c` (matmul). ODE/ODE_SYS/PDE/SDE
  (`src/ode.c`, `ode_sys.c`, `pde.c`, `sde.c`, `spectral.c`, …) are scalar;
  no `ml_vec4`/batch-rsqrt/FMA fast path exists for integrators or field
  solvers. Vectorizing them is future work, not a silent claim.

## 5. Orchestrator's Reverse-Engineering Log
- **2026-09-27 — Dispatch runtime fix owned:** `src/cpu_dispatch.c:69` guards
  `ml_matmul_avx2` with runtime `__builtin_cpu_supports("avx2")&&("fma")`
  plus scalar fallback; `ml_cpu_has_avx2/fma/sse41` use the same runtime
  query. Rationale in `ml_cpu_has_fma` comment: `__FMA__` compile macro
  reflects build flags, not host CPU — old macro test risked SIGILL.
  Verified with grep; `cpu_dispatch.h` header documents compile-time +
  runtime split.
- **2026-09-27 — Bench rdtsc carve-out owned:** `benchmarks/bench.c:12-26`
  x86 guard + `clock()` fallback, `BENCH_UNITS` cycles/ticks, excluded from
  lib build. Only inline asm in the fast-math story; library code stays
  zero-dependency C.
- **2026-09-27 — Profile routing owned:** `profiles.h:8-18` GRAPHICS routes
  `ml_rsqrt→ml_fast_rsqrt`, default routes `→1.0/ml_sqrt`, EMBEDDED hides
  fast math; inline-function (not macro) form per double-evaluation fix.
  Matches `bench.c:35-41` active-profile print.
- **2026-09-27 — Batch fallback + tolerances owned:** `simd_batch.h:40-48`
  whole-batch scalar fallback verified against
  `test_edge_redteam_p1_1.c` (signed-zero/Inf/NaN/subnormal lanes) and
  `fuzz_god_mode.c:352-381` (1e-5 batch, 1e-3 bare-metal); `fma_audit.s:38`
  `vfmadd213sd` confirms contracted-FMA codegen for the Horner/FMA claims.
