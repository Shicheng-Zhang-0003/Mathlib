# MathLib Design Contract (v11A1)

This document defines the strict architectural and operational boundaries of MathLib.
Any future contributions must adhere to these policies.

## 1. Memory Policy
* **Zero Internal Allocation:** Core APIs (`ml_solve_v10`, `ml_fft_execute`, etc.) will **never** call `malloc`, `calloc`, or `free`.
* **Client-Provided Scratchpads:** Any operation requiring temporary memory must accept a `ml_workspace_t` bump-allocator from the caller.
* **Legacy Isolation:** Modules requiring heap allocation are strictly quarantined in the `legacy/` directory and are not part of the core static library.

## 2. Threading Policy
* **Stateless & Thread-Safe:** All core math functions are pure, stateless, and rely only on local stack variables and caller-provided workspaces.
* **No Global State:** There are no global variables, no hidden caches, and no thread-local storage. MathLib is inherently thread-safe.

## 3. Determinism Policy
* **Bit-Identical Across Profiles:** The `SCIENTIFIC` profile guarantees deterministic output regardless of the underlying SIMD execution path (AVX2 vs Scalar fallback).
* **No Fast-Math Compiler Flags:** The build system strictly enforces `-fno-fast-math` and `-ffp-contract=off` to prevent the compiler from reordering IEEE-754 operations.

## 4. Error Handling Policy
* **Core Math:** Returns standard IEEE-754 `double` (using `NaN` and `Inf` for domain errors, matching `glibc libm`).
* **Complex Operations:** Returns `ml_status_t` enum (e.g., `ML_ERR_SINGULAR`, `ML_ERR_WORKSPACE`) to explicitly signal structural failures.

## 5. Precision Policy
* **Software Bounds:** Pure C99 software implementations (without hardware FMA) guarantee `<= 5 ULP` deviation from ground-truth `glibc libm`.
* **Exact Operations:** Bitwise parsers (`ml_isnan`, `ml_fabs`) and Error-Free Transformations (`ml_two_sum`) are 100% IEEE-754 exact with zero branching.
