# Build and Installation (v12R3 development, 59 TUs)

Recommended build:

cmake -S . -B build -DMATHLIB_PROFILE=SCIENTIFIC
cmake --build build

Profiles: `SCIENTIFIC` (accurate `1/sqrt`), `GRAPHICS` (`ml_fast_rsqrt`),
`EMBEDDED` (accurate fallback; fixed-point CORDIC path for trig).
All three profiles build the 58 listed translation units. GCC and Clang are
the validated toolchain families; the CMake configuration fails early for
unconfigured compiler flag sets. CMake 3.14 or later is required.

Strict flags (CMake + Makefile aligned):

-Wall -Wextra -Wconversion -Wshadow -Wpedantic -Werror -std=c99
-fno-fast-math -ffp-contract=off -Iinclude/mathlib -Isrc

Canonical TU list (58) must match across `CMakeLists.txt`, `Makefile`,
`run_all_tests.py`, `tests/run_edge_tests.sh` (kelvin drift fixed).

Sanitizers: `-DMATHLIB_SANITIZERS=ON` (ASan+UBSan, compile+link).
C++: all public headers carry `extern "C"` guards.

Install and consume the exported static target:

```bash
cmake --install build --prefix /path/to/prefix
```

Downstream CMake projects can use `find_package(mathlib 12.3 REQUIRED)` and
link `mathlib::mathc`. On Unix, the exported target carries its `libm`
dependency so static-library consumers do not need to repeat it.

Before deployment run the closure gate (`closure_gate.sh`: ASan configure,
build, modular, edge with sanitizers, boundary, oracle, ultimate fuzzer).
ULP work needs x86-64 long double (`LDBL_MANT_DIG>=64`); MSVC/ARM falls back
to DD (guarantees stay at ≤5 ULP gate).
