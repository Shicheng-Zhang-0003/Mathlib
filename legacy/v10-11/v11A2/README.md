# MathLib V1.0
![C99 Standard](https://img.shields.io/badge/C-99-blue.svg)
![Zero Allocation](https://img.shields.io/badge/Heap-Zero%20Alloc-brightgreen.svg)
![Invariant Audited](https://img.shields.io/badge/Fuzzer-65,000%2B%20Assertions-orange.svg)
![Bare Metal](https://img.shields.io/badge/SIMD-AVX2%20Optimized-red.svg)

**MathLib V1.0** is a high-performance, zero-dependency, bare-metal C99 scientific computing engine. 

Engineered to eliminate `<math.h>` dependencies, eradicate heap allocations in hot loops, and provide cycle-accurate determinism for high-throughput pipelines.

## 🚀 Quickstart (One-Command Install)

```bash
git clone https://github.com/yourname/mathlib.git
cd mathlib
cmake -B build -DMATHLIB_PROFILE=SCIENTIFIC
cmake --build build
```

To integrate into your own CMake project:
```cmake
find_package(mathlib REQUIRED)
target_link_libraries(your_app PRIVATE mathlib::mathc)
```

## 📊 API Precision Contracts
MathLib avoids overclaiming "strict IEEE-754" across the board. Instead, we guarantee specific mathematical contracts per module:

| Module | Guarantee | Implementation |
| :--- | :--- | :--- |
| **Bitwise Parsers** (`ml_isnan`, `ml_fabs`) | **100% IEEE-754 Exact** | Pure 64-bit integer bitmasking. Zero branching. |
| **Transcendentals** (`ml_sin`, `ml_exp`) | **< 1 ULP Error** | 19th-degree Maclaurin series + Payne-Hanek reduction. |
| **Fast Math** (`ml_fast_rsqrt`) | **< 10^-4 Relative Error** | Quake III integer-cast bit-hack + 2x Newton-Raphson. |
| **Linear Algebra** (`ml_solve_v10`) | **Zero Heap Allocation** | Client-provided scratchpad bump-allocator (`ml_workspace_t`). |

## 🏗️ Architecture & Profiles
MathLib uses compile-time hardware profiling to route math functions to the optimal implementation:

*   **`SCIENTIFIC` (Default):** Strict precision, 19th-degree polynomials, Payne-Hanek range reduction.
*   **`GRAPHICS`:** Blazing fast AVX2 SIMD batch processing and Quake III `rsqrt` bit-hacks.
*   **`EMBEDDED`:** True Q16.16 Fixed-Point CORDIC. Zero floats, zero FPU dependencies.

## 🧪 Testing & Invariant Fuzzing
MathLib ships with a deterministic smoke test and two dedicated property-based invariant fuzzers.

```bash
make test           # Runs the lean CI/CD smoke test
make fuzz_god_mode  # Unleashes 65,000+ randomized algebraic & IEEE-754 assertions
make fuzz_boundary  # Tests exact mathematical boundaries and energy conservation
```

## 🔄 Migration Guide (Linker-Level Wrapping)
If you have an existing codebase heavily reliant on standard `libm`, you do not need to rewrite your code. Use the GCC linker wrap flag to seamlessly route standard calls into MathLib:

```bash
gcc your_app.c -Wl,--wrap=sin -Wl,--wrap=cos -lmathc -lm
```

## 📜 License
MIT License. Free for commercial and open-source use.