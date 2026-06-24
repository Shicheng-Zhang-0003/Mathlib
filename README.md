# libmathc ⚡
**A bare-metal, hardware-optimized C99 scientific computing engine.**

![C99](https://img.shields.io/badge/C-C99-blue)
![CMake](https://img.shields.io/badge/build-CMake-green)
![Tests](https://img.shields.io/badge/tests-1300%2B%20Fuzzed-brightgreen)
![License](https://img.shields.io/badge/license-MIT-yellow)

`libmathc` is a comprehensive, from-scratch mathematical library written in C. It covers everything from basic combinatorics and $O(1)$ trigonometry to hardware-level AVX SIMD intrinsics, IEEE 754 bit-manipulation, and advanced NxN matrix decompositions (LU/QR/Eigenvalues).

## 🚨 Transparency Notice: The "Vibe Coding" Journey
**I wrote v1 entirely by hand.**

`v2` through `v8` were **100% vibe-coded**. I acted as the Prompt Engineer, while an AI pair-programmer wrote, compiled, and debugged the actual C syntax, Python fuzzing scripts, and CMake configurations.

I directed the mathematical constraints, discovered the edge-case failures, and designed the simpler algorithmic architecture.

---

## 🚀 Key Features

### 1. Pure Mathematics & Calculus
*   **$O(1)$ Trigonometry:** Replaced slow `while` loops with mathematical identities and `fmod`.
*   **Base-2 Split Exponentials:** $e^x = 2^n \cdot e^r$ for rapid convergence.
*   **Fast Logarithms:** $O(1)$ exponent extraction via `frexp` and ultra-fast $z = \frac{m-1}{m+1}$ Taylor series.
*   **Numerical Methods:** Newton-Raphson, Bisection, Simpson's Rule integration, and Runge-Kutta 4 (ODEs).

### 2. Advanced Linear Algebra (NxN)
*   **Struct-based API:** Beautiful, pass-by-value C99 structs (`vec2`, `vec3`, `mat3x3`).
*   **LU Decomposition:** With partial pivoting for numerically stable $Ax=b$ solving.
*   **Householder QR:** For least-squares regression and orthogonal factorization.
*   **Eigenvalue Extraction:** Power Iteration and the **Wilkinson-Shifted QR Algorithm** for finding all roots of the characteristic polynomial.

### 3. Signal Processing & Optimization
*   **Cooley-Tukey Radix-2 FFT:** In-place Fast Fourier Transform and Inverse FFT for signal processing.
*   **Optimization Solvers:** Golden Section Search and Gradient Descent.

### 4. Silicon & Bitwise Math (v8)
*   **CORDIC Engine:** Shift-and-add trigonometry (no multiplication in the inner loop) with full 360-degree quadrant mapping.
*   **Pure IEEE 754 Bit-Masking:** Custom `ml_frexp_pure` and `ml_ldexp_pure` using 64-bit integer unions, completely bypassing `<math.h>`.
*   **Minimax Polynomials:** Remez-generated coefficients evaluated via Horner's Method.
*   **Payne-Hanek Reduction:** 3-part Cody-Waite extended precision for zero-loss trigonometric range reduction on massive angles.
*   **AVX SIMD:** Raw `_mm256_mul_pd` intrinsics for 4-wide double-precision vector math.
*   **Error-Free Transformations:** Dekker/Knuth algorithms and Software FMA to capture microscopic IEEE 754 rounding errors.

---

## 📊 Benchmarks
During Phase 3, the trigonometry engine was upgraded from $O(N)$ loops to $O(1)$ algorithmic reduction.

**Task:** Calculate `sine(x)` 100,000 times on massive angles (`x > 10,000.0`).
*   **v1 / v2 (while loops):** ~450.00 ms
*   **v3+ (O(1) fmod):** ~1.05 ms *(400x speedup)*

---

## 🧪 The Chaos Monkey Test Suite
Math libraries live and die by their edge cases. `test.c` includes an **Ultimate Fuzzing Suite** that dynamically generates **~1,300 randomized tests** via Python and C, verifying the library against IEEE 754 standards.

It explicitly tests:
*   Astronomical angles (e.g., $10^6$ radians) to verify range reduction.
*   Subnormal numbers, `NaN`, and `Infinity` propagation.
*   Singular matrices and exact eigenvalue convergence.
*   Software FMA capturing swallowed `1e-16` rounding errors.

---

## 🛠️ Building and Installation

`libmathc` uses CMake to compile into a static library (`libmathc.a`).

### Prerequisites
*   GCC or Clang with C99 support
*   CMake (3.10+)
*   A CPU with AVX support (for the SIMD module)

### Build Steps
```bash
mkdir build && cd build
cmake ..
make
