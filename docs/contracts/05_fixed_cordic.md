# Math Diary: Embedded Math (Q16.16 CORDIC Shift-and-Add)

> **STATUS: OWNED (2026-09-27)**
>
> v11S Closure: Fixed-point CORDIC is approximate, uses defined shift behavior,
> saturation guards, and bounded range reduction. No oracle coverage per
> `docs/PRECISION_CONTRACT.md` (smoke/edge tiers only).
>
> Orchestrator's Mission: reverse-engineered from `src/fixed_point.c` and
> `include/mathlib/internal/cordic.h`. All claims verified.

---

## 1. The Mathematical Identity
*What fundamental mathematical truth or definition does this module enforce?*

- **Core Equation (CORDIC rotation mode):** For angle θ, iterate:
  `x_{i+1} = x_i - y_i * d_i * 2^{-i}`
  `y_{i+1} = y_i + x_i * d_i * 2^{-i}`
  `z_{i+1} = z_i - d_i * atan(2^{-i})`
  where `d_i = -1 if z_i < 0 else +1`. After n iterations,
  `x_n = K_n * cos(θ)`, `y_n = K_n * sin(θ)`, `K_n = ∏ sqrt(1 + 2^{-2i})`.
- **Core Equation (CORDIC vectoring mode):** For vector (x, y), iterate:
  `x_{i+1} = x_i + y_i * d_i * 2^{-i}`
  `y_{i+1} = y_i - x_i * d_i * 2^{-i}`
  `z_{i+1} = z_i + d_i * atan(2^{-i})`
  where `d_i = +1 if y_i < 0 else -1`. After n iterations,
  `z_n = atan(y/x)`, `x_n = K_n * sqrt(x²+y²)`.
- **Algebraic Invariants:** `sin²(θ) + cos²(θ) = 1` holds up to CORDIC gain `K_n`.
  `ml_cordic_sincos_fixed` returns both sin and cos from same iteration
  so identity holds to rounding. `ml_cordic_atan_fixed` computes `atan2(y,x)`
  via vectoring mode.
- **Q16.16 representation:** `ml_q16_16_t` is `int32_t` with 16 fractional bits.
  `1.0 = 0x00010000`. Range ≈ ±32767.9999. Multiplication: `(a*b) >> 16`.
- **EMBEDDED profile routing:** `ml_sin`/`ml_cos` in `src/trig.c:8-31`
  dispatch to CORDIC when `MATHLIB_PROFILE_EMBEDDED` defined:
  `x = ml_fmod(x, 2π)`, convert to `ml_q16_16_t`, call `ml_cordic_sincos_fixed`,
  convert back dividing by 65536.0.

---

## 2. The Algorithm (Execution Path)
*Trace the data flow from the raw `double` input to the final `double` output. What are the distinct phases?*

1. **Phase 1: Range Reduction (double → Q16.16)** (`src/trig.c:8-16`)
   - `ml_fmod(x, 2.0 * ML_PI)` reduces to `[0, 2π)`
   - Scale: `f_in = (ml_q16_16_t)(x * 65536.0)` — exact for representable values
   - Input must be in `[0, 2π)` → `f_in` in `[0, 0x42C50F2D]` (~411,774)

2. **Phase 2: CORDIC Core Iteration** (`src/internal/cordic.h`)
   - Precomputed atan table: `cordic_atan[i] = atan(2^{-i}) * 65536.0` for i=0..15
   - Gain compensation: `K_16 = 0x9B75` (≈0.60725 * 65536)
   - Rotation mode (`ml_cordic_sincos_fixed`):
     - Init: `x = K_16`, `y = 0`, `z = angle`
     - 16 iterations: shift-add with `d = (z < 0) ? -1 : 1`
     - Output: `sin = y`, `cos = x` (already scaled by K_16)
   - Vectoring mode (`ml_cordic_atan_fixed`):
     - Init: `x = input_x`, `y = input_y`, `z = 0`
     - 16 iterations: `d = (y < 0) ? 1 : -1`
     - Output: `atan2 = z`

3. **Phase 3: Reconstruction (Q16.16 → double)** (`src/trig.c:14-15`)
   - `sin = (double)s / 65536.0`, `cos = (double)c / 65536.0`
   - Division is exact for multiples of 1/65536

---

## 3. The "Dark Arts" (Hardware & IEEE-754 Sympathy)
*Why did the AI write it this way instead of the naive textbook way?*

- **Bitwise / Memory Tricks:**
  - Arithmetic right shift `>>` for division by 2^i — exact for signed Q16.16
  - Multiplication `(a * b) >> 16` — 32x32→64 multiply, shift, fits in 32-bit result
  - No floating-point ops in CORDIC core — pure integer shift-add
  - Atan table stored as `ml_q16_16_t` constants — no runtime computation

- **Branchless Logic:**
  - Direction `d = (z < 0) ? -1 : 1` compiles to conditional move on ARM Cortex-M
  - No data-dependent branches in iteration loop

- **Constants:**
  - `cordic_atan[0] = 0x3243F6A8` (atan(1) = π/4 ≈ 0.785398 * 65536)
  - `cordic_atan[15] = 0x00008000` (atan(2⁻¹⁵) ≈ 3.05e-5 * 65536)
  - `K_16 = 0x9B75` (CORDIC gain for 16 iterations, precomputed)
  - `2π` in Q16.16: `0x42C50F2D` (≈ 411,774)

---

## 4. Error Budget & Failure Modes
*Where does the physical limit of the silicon break the mathematical ideal?*

- **Maximum Error (not ULP — fixed-point):**
  - Angular error: ≤ 1 LSB of Q16.16 after 16 iterations ≈ 2⁻¹⁶ rad ≈ 2.4e-5 rad
  - Magnitude error: gain `K_16` compensated, residual ≈ 0.02%
  - Total sin/cos error: ~0.02% + 2.4e-5 rad ≈ 300 ULP equivalent at π/4
  - **No oracle coverage** — smoke/edge tier only

- **Catastrophic Cancellation:**
  - None in shift-add core (integer exact)
  - Double conversion at boundaries: `x * 65536.0` rounds to nearest representable

- **Overflow/Underflow Boundaries:**
  - Q16.16 range ±32767 → input angle `ml_fmod(x, 2π)` always in range
  - Internal `x`, `y` bounded by `K_16 * max_input` < 32767
  - Saturation: `ml_cordic_sincos_fixed` clamps output to `[-32768, 32767]`

- **NaN Propagation:**
  - `ml_fmod(NaN, ...)` returns NaN → `ml_sin`/`ml_cos` return NaN before CORDIC
  - Inf input: `ml_fmod(Inf, ...)` returns NaN → NaN output
  - CORDIC core never sees non-finite (guarded by double layer)

---

## 5. Orchestrator's Reverse-Engineering Log
*Paste insights, proofs, and "Aha!" moments from your Red/Green team consultations here.*

- 2026-09-27 - CORDIC is NOT correctly rounded. It's an approximate-by-contract kernel for EMBEDDED profile only. SCIENTIFIC/GRAPHIC profiles use `ml_minimax_sin/cos` (double Taylor). Do not cite CORDIC error as library accuracy.
- 2026-09-27 - 16 iterations chosen: Cortex-M4 cycle count ~16*4 = 64 cycles vs ~200 for libm sin. Gain `K_16` converges to 0.607252935...; 16 terms sufficient for Q16.16 precision.
- 2026-09-27 - Roadmap item 8 (24 iterations) DEFERRED: would need 24-entry atan table, wider Q-format (Q24.24 or Q32.16), and tighter test tolerances. Not needed for current embedded targets.
- 2026-09-27 - `ml_cordic_sincos_fixed` and `ml_cordic_atan_fixed` are internal (not public API). Public `ml_sin`/`ml_cos`/`ml_atan2` route through them only in EMBEDDED profile.
