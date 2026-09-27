# Math Diary: Exponentials & Logarithms (Horner's & Cody-Waite)

> **STATUS: OWNED (2026-09-27)**

> History: v11S closure hardened Cody-Waite/Horner pow + hyperbolics; v12A1 validated Taylor-19/atanh-11 kernels (report); v12R2 fixed pow limit, FMA fallback, Dekker split, hex constants.

---

## 1. The Mathematical Identity
*What fundamental mathematical truth or definition does this module enforce?*
- **Core Equation (exp):** `exp(x) = 2^n * exp(r)` with `n = round(x / ln2)`, `r = x - n*ln2` in `[-ln2/2, ln2/2]`. Implemented in `ml_exp` (`src/exp_log.c:43`).
- **Core Equation (log):** `log(x) = e*ln2 + 2*atanh(z)`, `z = (m-1)/(m+1)` after `ml_frexp_pure` mantissa `m` in `[0.7071, 1.4142)` (adjust `m *= 1+adjust; e -= adjust`). The atanh series `2*(z + z^3/3 + ... + z^21/11)` is 11 terms (`lc[0..10]`, `src/exp_log.c:139-144`), the exact kernel measured at 1 ULP in `docs/V12A1_MINIMAX_REPORT.md`.
- **Core Equation (pow):** general path `pow(x,y) = exp(y*log(x))` computed in double-double (`ml_log_split` + DD product + `ml_exp` + second-order `elo` correction, `src/exp_log.c:320-338`) PLUS an exact integer path: if `ml_is_integer_double(y)` (`src/internal/pow_util.h:11`) and `|y| <= 1023`, binary exponentiation (`result *= base` / squaring loop, `src/exp_log.c:285-296`) is exact — `pow(2,10)=1024`, `pow(10,3)=1000` with no log/exp roundtrip. Negative integer bases recurse on `ml_pow(-x, y)` with odd/even sign via `ml_is_odd_integer_double` (`src/exp_log.c:304-314`).
- **Algebraic Invariants:** `exp(log(x)) == x` (x > 0), `log(exp(x)) == x` (finite), `pow(x,0) == 1` even for `x = 0/NaN` (except `1^NaN == 1` per IEEE), `exp2(x) = 2^x` via exact `ml_ldexp_pure`, `log2/log10` via DD division by `ML_LN2`/`LN10_HI+LO`, `ml_logb(x,b) = ml_log(x)/ml_log(b)`.
- **Hyperbolic thresholds:** `ml_sinh` returns `x` for `|x| < 1e-8` (cubic term below 1 ULP); `ml_tanh`/`ml_asinh`/`ml_atanh` return `x` for `|x| < 1.5e-8` (release_notes.md item 4: old `1e-4` cost ~3.3M ULP at `x=1e-4`).

## 2. The Algorithm (Execution Path)
*Trace the data flow from the raw `double` input to the final `double` output. What are the distinct phases?*
1. **Phase 1: Range Reduction / Pre-processing**
   - *Code location:* `ml_exp` reduction (`src/exp_log.c:67-87`); `ml_log`/`ml_log_split` frexp+adjust (`src/exp_log.c:128-136,194-198`); `ml_pow` IEEE guards + integer path (`src/exp_log.c:229-296`).
   - *Why we do it:* the raw input cannot go into the polynomial: `exp` Taylor only converges fast on `[-ln2/2, ln2/2]`, so `n = ml_round(x / ML_LN2)` then the FMA Cody-Waite residual `r = ML_FMA(-n, ML_LN2_HI, x); r = ML_FMA(-n, ML_LN2_LO, r)` keeps ~106 bits (old `x - n*hi - n*lo` lost 2 ULP). `log` needs `m` near 1 so `|z| <= 0.06` and the atanh series converges in 11 terms; `pow` must dispatch NaN/Inf/zero first because `exp(y*log(x))` is undefined there.
2. **Phase 2: The Core Approximation**
   - *Code location:* `ml_exp` Taylor-19 DD Horner (`src/exp_log.c:89-107`); `ml_log` atanh DD Horner (`src/exp_log.c:145-149`); `ml_expm1` Taylor+Kahan / `ml_log1p` atanh (`src/exp_log.c:498-564`); `ml_sinh` Taylor to x^19/19! for `|x| < 0.5` (`src/exp_log.c:369-383`).
   - *Why we do it:* `exp`: `inv_fact[0..19]` (`1/k!`) evaluated as `ml_ddx_from_d(inv_fact[19])` then `ml_ddx_mul_d(acc, r)` / `ml_ddx_add_d(acc, inv_fact[i])` down to `i=1`, times `r`, plus 1.0 — DD Horner, <0.5 ULP eval. `log`: `lc[0..10]` (`2, 2/3, 2/5, ..., 2/21`) Horner in `t = z*z` via the same DD helpers (was 2-rounding mul+add). `ml_expm1` for `|x| < 0.5` uses Taylor with Kahan summation and the quadratic return-x threshold `|x| < 1.11e-16` (1.5e-8 is cubic-only, cost 2.4M ULP at `x=1e-9` — see comment at `src/exp_log.c:503-506`); `ml_log1p` for `|x| < 0.5` uses `z = x/(2+x)`, `z*P(z^2)` DD Horner. `ml_tanh` uses `(1-e)/(1+e)` with `e = ml_exp(-2*|x|)` and early `|x|>20 -> +/-1`, `|x|<1.5e-8 -> x`; `ml_asinh`/`ml_atanh` share the `1.5e-8` return-x cutoff.
3. **Phase 3: Reconstruction / Post-processing**
   - *Code location:* `ml_exp` ldexp split (`src/exp_log.c:108-116`); `ml_log` DD `e*ln2 + z*P` (`src/exp_log.c:152-164`); `ml_pow` DD product + `ml_exp` (`src/exp_log.c:320-338`); `ml_exp2`/`ml_sinh`/`ml_cosh` shifted overflow paths.
   - *Why we do it:* `exp`: `p = ml_ddx_to_d(acc)` then `ml_ldexp_pure(p, ni)`; the `ni > 1023` split (`ml_ldexp_pure(p,1023) * ml_ldexp_pure(1.0, ni-1023)`) avoids spurious overflow at `n == 1024` (`x ~ 709.78`) where `ldexp(p,1024) = Inf` but true `exp(x) <= DBL_MAX`. `log`: `ehi = ed*ML_LN2_HI; elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed*ML_LN2_LO`, `ml_ddx_renorm` + `ml_two_sum` merge with `zp` — compensated reconstruction. `pow`: `p = y*log_hi; e = ML_FMA(y, log_hi, -p) + y*log_lo`, `PE = ml_ddx_renorm(p,e)`, `g = ml_exp(PE.hi)`, `elo = ML_FMA(PE.lo, PE.lo*0.5, PE.lo)+1.0`, result `ML_FMA(g, elo, 0.0)`. `ml_sinh`/`ml_cosh` shift by `ML_LN2` above 700.0 (`ml_exp(ax-ML_LN2)`) and saturate above `ML_LOG_HYP_OVERFLOW`.

## 3. The "Dark Arts" (Hardware & IEEE-754 Sympathy)
*Why did the AI write it this way instead of the naive textbook way?*
- **Bitwise / Memory Tricks:** [e.g., Why use `memcpy` instead of `union`? Why use Horner's Method?]
- **Branchless Logic:** [Are there masks or bitwise ops used to avoid CPU pipeline stalls?]
- **Constants:** [Where do the magic numbers (e.g., `0x5fe6ec85...`) come from?]

## 4. Error Budget & Failure Modes
*Where does the physical limit of the silicon break the mathematical ideal?*
- **Maximum ULP Error:** [ ] ULPs (Units in the Last Place) vs MPFR ground-truth.
- **Catastrophic Cancellation:** [Identify inputs where subtracting two nearly equal numbers destroys precision]
- **Overflow/Underflow Boundaries:** [At what exact input value does the result become `Inf` or `0.0`?]
- **NaN Propagation:** [How does the code handle `NaN` and `Inf` inputs?]

## 5. Orchestrator's Reverse-Engineering Log
*Paste insights, proofs, and "Aha!" moments from your Red/Green team consultations here.*
- [Date] - Insight 1: ...
- [Date] - Insight 2: ...
