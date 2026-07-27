# MathLib v11S API Status (Closure Candidate)
<!-- v11S CLOSURE IP-22: docs alignment -->

This document defines the public interface status for MathLib `11.0.0-rc1` (v11S closure candidate).

**No new features or signature changes are permitted beyond the v11S closure boundary.**

---

## Module Status

| Module / Header | Status | Notes |
| :--- | :--- | :--- |
| `ml_core.h` | **STABLE** | Bitwise IEEE-754 helpers, `ml_sqrt`, `ml_fmod`, `ml_round` |
| `ml_trig.h` | **STABLE** | Maclaurin kernels, bounded Cody-Waite reduction, NaN/Inf guards |
| `ml_exp_log.h` | **STABLE** | Cody-Waite reduction, Horner evaluation, hyperbolic edge hardening |
| `ml_complex.h` | **STABLE** | Overflow-safe abs, atan2-based arg, NaN guards |
| `fft.h` | **STABLE** | Power-of-two radix-2 Cooley-Tukey only |
| `ml_linalg.h` | **STABLE** | Zero-alloc LU solve, relative singularity threshold |
| `ml_tensor.h` | **STABLE** | Workspace bump allocator with fixed ABI and canary |
| `ml_statistics.h` | **STABLE** | Invalid-argument guards, finite-input validation |
| `ml_combinatorics.h` | **STABLE** | Portable overflow detection with `UINT64_MAX` sentinel |
| `ml_numerical.h` | **STABLE** | Root-finding, derivative, Simpson integration guards |
| `ml_optimization.h` | **STABLE** | Golden-section and gradient-descent guards |
| `ml_ode.h` | **STABLE** | Euler / RK4 guards |
| `ml_polynomial.h` | **STABLE** | Horner evaluation and Newton guards |
| `ml_quadratics.h` | **STABLE** | Citardauq-style stable quadratic roots |
| `ml_integral.h` | **EXPERIMENTAL** | Positive-domain gamma only; traditional integrator is experimental |
| `ml_fixed_point.h` | **STABLE** | Q16.16 CORDIC approximate trig with defined shifts |
| `ml_quaternion.h` | **STABLE** | Quaternion algebra and hardened slerp |
| `fast_math.h` | **STABLE** | Approximate fast paths with explicit domain guards |
| `profiles.h` | **STABLE** | Compile-time profile routing |
| `simd.h` / `simd_batch.h` | **STABLE** | SIMD and scalar fallback paths |
| `cpu_dispatch.h` | **STABLE** | Compile-time capability queries |
| `version.h` | **STABLE** | Version metadata |

---

## Error Handling Contract

* Core math functions return IEEE-754 `double`:
  - `NaN` for domain errors
  - `Inf` for overflow where appropriate

* Structural functions return `ml_status_t`:
  - solver failures
  - workspace exhaustion
  - invalid arguments
  - non-finite input rejection

---

## Removed / Legacy

* `compat.h` is not part of the v11S core.
* `legacy/` modules are not built by default and are outside the v11S stability boundary.
