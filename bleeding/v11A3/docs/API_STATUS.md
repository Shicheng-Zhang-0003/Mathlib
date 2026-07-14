# MathLib v11A2 API Status (Freeze)

This document defines the permanent public interface for MathLib v1.1.
**No new features or signature changes are permitted beyond this point.**

| Module / Header | Status | Notes |
| :--- | :--- | :--- |
| `ml_core.h` | **STABLE** | Bitwise IEEE-754, pure C fmod/round, hardware sqrt. |
| `ml_trig.h` | **STABLE** | Horner's Minimax, Payne-Hanek, Polymorphic API (`_scientific`, `_fixed`). |
| `ml_exp_log.h` | **STABLE** | Horner's Method, Cody-Waite split, hyperbolic functions. |
| `ml_linalg.h` | **EXPERIMENTAL** | Zero-alloc LU decomposition. API may shift in v1.2. |
| `ml_tensor.h` | **EXPERIMENTAL** | Workspace bump allocator. API may shift in v1.2. |
| `fft.h` | **STABLE** | Cooley-Tukey iterative, O(N) twiddle refresh. |
| `ml_complex.h` | **STABLE** | Pure C99 complex arithmetic. |
| `compat.h` | **LEGACY** | V1.0 shims. Throws `#warning`. Will be removed in v2.0. |
| `legacy/` | **REMOVED** | Quarantined malloc-heavy modules. Not built by default. |

## Error Handling Contract
* Core math (`ml_sin`, `ml_exp`) returns IEEE-754 `double` (`NaN`/`Inf` for domain errors).
* Complex operations (`ml_solve_v10`) return `ml_status_t` enum.
