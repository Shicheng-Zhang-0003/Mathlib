# Math Diary: Signal Processing (Cooley-Tukey FFT & Twiddle Refresh)

> **STATUS: OWNED (2026-09-27)**
>
> Note (condensed): v11S closure still applies — power-of-two only,
> unsupported lengths are safe silent no-ops, per-stage twiddle refresh bounds
> drift. Validation tier is smoke/edge only: FFT has no oracle coverage per
> `docs/PRECISION_CONTRACT.md`, so no ULP claim may be cited from FFT tests.

---

## 1. The Mathematical Identity
*The DFT and its exact inverse, plus the energy identity the smoke test checks.*

- **Core Equation (DFT):** `X[k] = sum_{n=0}^{N-1} x[n] W_N^{kn}`,
  `W_N = exp(-2πi/N)`, enforced by `ml_fft_execute` (`src/fft.c:15-72`).
  Forward transform is unscaled — no normalization flag exists.
- **Inverse (conjugate-forward-conjugate + 1/N):**
  `x[n] = (1/N) sum_{k=0}^{N-1} X[k] W_N^{-kn}`, implemented literally in
  `ml_ifft_execute` (`src/fft.c:74-95`): negate all `imag` (`82-84`), call
  `ml_fft_execute` (`86`), negate `imag` again and multiply both parts by
  `inv = 1.0/(double)ns` (`88-94`).
- **Algebraic Invariants:** `FFT(IFFT(x)) = x` up to rounding; DC bin
  `X[0] = sum x[n]` (`tests/test_dsp.c:10-13` checks impulse `->` all-ones
  spectrum); Parseval energy `sum_n |x[n]|² = (1/N) sum_k |X[k]|²`, checked by
  `test_fft_parseval` (`tests/fuzz_boundary_gauntlet.c:100-121`, n=64,
  `|E_time − E_freq/N| < E_time·1e-13`). Time-side energy helper is
  `ml_parseval_energy` (`src/transforms.c:48-59`, long-double accumulation,
  NaN on non-finite input or non-finite result).

---

## 2. The Algorithm (Execution Path)
*Trace the data flow from the raw `double` input to the final `double` output. What are the distinct phases?*

1. **Phase 1: Pre-check & Early Exit** (`ml_fft_execute:15-23`)
   - NULL pointer → silent no-op
   - Non-power-of-two `n` → silent no-op (safe, no UB)
   - `n > 2^24` (256MB) → silent no-op (DoS cap)

2. **Phase 2: Bit-Reversal Permutation** (`ml_fft_execute:27-42`)
   - Iterative bit-reversal: `for (i=1, j=0; i<ns; i++) { bit=ns>>1; while (j&bit) j^=bit, bit>>=1; j^=bit; if (i<j) swap(x[i],x[j]); }`
   - In-place, no extra memory

3. **Phase 3: Cooley-Tukey Iterative Radix-2** (`ml_fft_execute:44-71`)
   - Outer loop: `len = 2, 4, 8, ..., ns`
   - Twiddle base: `ang = -2π/len`, `wlen = {cos(ang), sin(ang)}`
   - Inner loop: stride `len`, process `half = len/2` butterflies
   - **Twiddle refresh every 64 steps** (`j & 63 == 0`): recompute `w = {cos(theta), sin(theta)}` from `theta = ang*j` instead of accumulating `w *= wlen`
   - Butterfly: `u = x[i+j]`, `v = x[i+j+half] * w`, `x[i+j] = u+v`, `x[i+j+half] = u-v`
   - `w = w * wlen` for next `j` (complex multiply via `ml_cplx_mul`)
   - Early break when `len > ns/2` (prevents unnecessary final shift)

4. **Phase 4: Inverse Transform** (`ml_ifft_execute:74-95`)
   - Conjugate input: `x[i].imag = -x[i].imag`
   - Forward FFT
   - Conjugate output: `x[i].imag = -x[i].imag`
   - Scale by `1/N`: `x[i].real *= inv`, `x[i].imag *= inv`

5. **Phase 5: Real FFT Helper** (`ml_fft_real:76-84` in `src/harmonic.c`)
   - Pack real input into complex buffer (imag=0)
   - Call `ml_fft_execute`
   - Unpack to separate `re[]`, `im[]` arrays

6. **Phase 6: 2D FFT** (`ml_fft2d_pow2:85-100` in `src/harmonic.c`)
   - Row-major layout, `rows[128][128]` stack buffer (n≤128)
   - FFT each row, then FFT each column (transpose via `cols[128]` scratch)
   - Dirichlet-only (no Neumann/Robin)

---

## 3. The "Dark Arts" (Hardware & IEEE-754 Sympathy)
*Why did the AI write it this way instead of the naive textbook way?*

- **Bitwise / Memory Tricks:**
  - Power-of-two check: `(n & (n-1)) == 0` — single instruction, no loops
  - `size_t` loop variables (`ns`, `i`, `j`, `len`, `half`) — `size_t` overflow impossible for `n <= 2^24`
  - In-place transform — no allocation, cache-friendly
  - `cplx` struct `{double real, imag}` — matches `complex.h` layout for potential SIMD

- **Branchless Logic:**
  - Twiddle refresh at fixed stride (`j & 63 == 0`) — bounds drift from accumulated `w *= wlen` (rounding error grows as O(√N) vs O(N) without refresh)
  - No `if` inside innermost butterfly loop — pure FMA-friendly complex arithmetic

- **Constants:**
  - `ML_PI` for `2π/len` angle computation — exact hex `0x1.921fb54442d18p+1` in header
  - Twiddle refresh period 64: empirically chosen — drift < 1e-15 at N=4096, cost ~1.5% extra cos/sin calls

---

## 4. Error Budget & Failure Modes
*Where does the physical limit of the silicon break the mathematical ideal?*

- **Maximum ULP Error:** No oracle coverage (smoke/edge tier only). Parseval energy test passes at `1e-13` relative. No ULP claim may be cited.
- **Catastrophic Cancellation:**
  - Butterfly `u-v` for nearly-equal inputs — inherent to DFT structure
  - Twiddle refresh limits phase error accumulation
  - Long-double energy accumulation in `ml_parseval_energy` avoids sum-of-squares overflow
- **Overflow/Underflow Boundaries:**
  - DC bin `X[0] = sum x[n]` can overflow for large `N` and large inputs — no scaling in forward transform
  - Inverse scaling `1/N` can underflow for large `N` — `inv = 1.0/N` computed in double
  - Input `n > 2^24` rejected before any compute (DoS prevention)
- **NaN Propagation:**
  - Non-finite input → NaN output (complex multiply `u*v` propagates NaN)
  - Silent no-op on invalid `n` — no NaN returned, input unchanged
  - `ml_parseval_energy` returns NaN on any non-finite input or result

---

## 5. Orchestrator's Reverse-Engineering Log
*Paste insights, proofs, and "Aha!" moments from your Red/Green team consultations here.*

- 2026-09-27 - Twiddle refresh rationale: without refresh, `w *= wlen` accumulates rounding error; at N=4096, phase drift ~1e-12 rad → Parseval error ~1e-9. Refresh every 64 steps bounds drift to <1e-15.
- 2026-09-27 - Silent no-op contract: `ml_fft_execute` does not return error code for invalid `n` — it's a void function. Callers MUST check `ml_fft_is_supported(n)` first. This is a deliberate v11S closure choice: zero-overhead hot path, validity is caller's responsibility.
- 2026-09-27 - Inverse uses conjugate trick: avoids separate IFFT code path, guarantees `FFT(IFFT(x)) = x` identity up to rounding (same code path both ways).
- 2026-09-27 - 2D FFT uses stack buffer `rows[128][128]` (max 128x128=16384 complex = 256KB). Thread-safe (stack-local), no heap. Larger n rejected (`ML_ERR_INVALID_ARG`).
- 2026-09-27 - `ml_fft_real` and `ml_fft2d_pow2` live in `harmonic.c` (Batch-1), not `fft.c` (core). They reuse core `ml_fft_execute`.
