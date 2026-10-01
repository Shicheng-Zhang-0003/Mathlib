# Build and Installation (v12R2 / V1.2-RC2, 34 TUs)

Recommended build:

cmake -S . -B build -DMATHLIB_PROFILE=SCIENTIFIC
cmake --build build

Profiles: `SCIENTIFIC` (accurate `1/sqrt`), `GRAPHICS` (`ml_fast_rsqrt`),
`EMBEDDED` (accurate fallback; fixed-point CORDIC path for trig).
All three compile `-Werror` clean (34 TUs incl. kelvin).

Strict flags (CMake + Makefile aligned):

-Wall -Wextra -Wconversion -Wshadow -Wpedantic -Werror -std=c99
-fno-fast-math -ffp-contract=off -Iinclude/mathlib -Isrc

Canonical TU list (34) must match across `CMakeLists.txt`, `Makefile`,
`run_all_tests.py`, `tests/run_edge_tests.sh` (kelvin drift fixed).

Sanitizers: `-DMATHLIB_SANITIZERS=ON` (ASan+UBSan, compile+link).
C++: all public headers carry `extern "C"` (verified `g++ -std=c++17`).

Before deployment run the closure gate (`closure_gate.sh`: ASan configure,
build, modular, edge with sanitizers, boundary, oracle, ultimate fuzzer).
ULP work needs x86-64 long double (`LDBL_MANT_DIG>=64`); MSVC/ARM falls back
to DD (guarantees stay at ≤5 ULP gate).
