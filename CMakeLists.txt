cmake_minimum_required(VERSION 3.14)
project(mathlib_v12R3 VERSION 12.3.0 LANGUAGES C)

# Strict C99 enforcement (No GNU extensions like typeof or statement expressions)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Profile Validation & Definition
set(MATHLIB_PROFILE "SCIENTIFIC" CACHE STRING "Mathlib Execution Profile")
if(NOT MATHLIB_PROFILE MATCHES "^(SCIENTIFIC|GRAPHICS|EMBEDDED)$")
    message(FATAL_ERROR "MATHLIB_PROFILE must be SCIENTIFIC, GRAPHICS, or EMBEDDED. Got: ${MATHLIB_PROFILE}")
endif()

if(MATHLIB_PROFILE STREQUAL "GRAPHICS")
    add_compile_definitions(MATHLIB_PROFILE_GRAPHICS)
elseif(MATHLIB_PROFILE STREQUAL "EMBEDDED")
    add_compile_definitions(MATHLIB_PROFILE_EMBEDDED)
else()
    add_compile_definitions(MATHLIB_PROFILE_SCIENTIFIC)
endif()

# Strict Compiler Flags (Matching the hardened Makefile)
# NOTE: -Wfloat-equal is intentionally omitted (The Float-Equal Paradox).
# In numerical computing (libm), exact comparisons to 0.0, 1.0,
# and NaN checks (x != x) are mathematically intentional and required.
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-O3 -Wall -Wextra -Wconversion -Wshadow -Wpedantic
                        -Wundef -Werror -fno-fast-math -ffp-contract=off)
elseif(MSVC)
    add_compile_options(/O2 /W4 /WX /fp:strict)
else()
    message(FATAL_ERROR "Strict compiler flags are configured for GCC, Clang, or MSVC; got ${CMAKE_C_COMPILER_ID}")
endif()

option(MATHLIB_NATIVE "Enable -march=native (Warning: Breaks binary portability)" OFF)
if(MATHLIB_NATIVE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        add_compile_options(-march=native)
    else()
        message(FATAL_ERROR "MATHLIB_NATIVE is only configured for GCC/Clang toolchains")
    endif()
endif()

# Sanitizers (Fixed: Added to link options as well)
option(MATHLIB_SANITIZERS "Enable Address and Undefined Behavior Sanitizers" OFF)
if(MATHLIB_SANITIZERS)
    if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        message(FATAL_ERROR "MATHLIB_SANITIZERS requires GCC or Clang")
    endif()
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    message(STATUS "Sanitizers enabled")
endif()

# The Real v11S Translation Units
set(MATHLIB_SOURCES
    src/core.c
    src/trig.c
    src/exp_log.c
    src/complex.c
    src/linalg.c
    src/fft.c
    src/cpu_dispatch.c
    src/combinatorics.c
    src/quadratics.c
    src/polynomial.c
    src/numerical.c
    src/statistics.c
    src/integral.c
    src/ode.c
    src/optimization.c
    src/quaternion.c
    src/fixed_point.c
    src/orthogonal.c
    src/calculus.c
    src/numbertheory.c
    src/transforms.c
    src/optim_n.c
    src/ode_sys.c
    src/spectral.c
    src/stats_inv.c
    src/sde.c
    src/pde.c
    src/harmonic.c
    src/mcmc.c
    src/manifold.c
    src/info.c
    src/analytic_nt.c
    src/control.c
    src/kelvin.c
    src/contest.c
    src/geometry.c
    src/graph.c
    src/modeling.c
    src/measure.c
    src/funcanal.c
    src/topology.c
    src/diffgeom.c
    src/algebraic.c
    src/math_extra.c
    src/bigint.c
    src/exact.c
    src/interval.c
    src/polyq.c
    src/gfp.c
    src/special.c
    src/linalg_extra.c
    src/calculus_extra.c
    src/stats_extra.c
    src/graph_extra.c
    src/optim_extra.c
    src/geom_extra.c
    src/complexan.c
    src/ode_extra.c
    src/school.c
)

add_library(mathc STATIC ${MATHLIB_SOURCES})
# A sanitizer-instrumented static archive needs the runtime at the final link.
# add_link_options() above covers this project, while the interface option
# carries that requirement through the exported target to package consumers.
if(MATHLIB_SANITIZERS)
    target_link_options(mathc INTERFACE -fsanitize=address,undefined
                        -fno-omit-frame-pointer)
endif()
# The static archive leaves libm symbols for the final executable to resolve.
# Publish that dependency so installed/imported consumers link correctly too.
if(UNIX)
    target_link_libraries(mathc PUBLIC m)
endif()
target_include_directories(mathc PRIVATE src)
target_include_directories(mathc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/mathlib>
    $<INSTALL_INTERFACE:include>
)

# Core Test Executables
add_executable(test tests/test.c)
target_link_libraries(test mathc)
target_include_directories(test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(bench benchmarks/bench.c)
target_link_libraries(bench mathc)
target_include_directories(bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(fuzz_god_mode tests/fuzz_god_mode.c)
target_link_libraries(fuzz_god_mode mathc)
target_include_directories(fuzz_god_mode PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(fuzz_boundary tests/fuzz_boundary_gauntlet.c)
target_link_libraries(fuzz_boundary mathc)
target_include_directories(fuzz_boundary PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# R5 despot: first-class ultimate_fuzzer target (was ad-hoc gcc lines in
# closure_gate.sh vs run_all_tests.py that diverged: -Lbuild -lmathc vs
# direct LIB_SOURCES compile). Single source of truth now.
add_executable(ultimate_fuzzer tests/ultimate_fuzzer.c)
target_link_libraries(ultimate_fuzzer mathc)
target_include_directories(ultimate_fuzzer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Modular CI/CD Tests
add_executable(test_core tests/test_core.c)
target_link_libraries(test_core mathc)
target_include_directories(test_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_trig tests/test_trig.c)
target_link_libraries(test_trig mathc)
target_include_directories(test_trig PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_linalg tests/test_linalg.c)
target_link_libraries(test_linalg mathc)
target_include_directories(test_linalg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_dsp tests/test_dsp.c)
target_link_libraries(test_dsp mathc)
target_include_directories(test_dsp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Examples & Validation (Hidden behind options to prevent CI bloat)
option(MATHLIB_BUILD_EXAMPLES "Build example programs" OFF)
if(MATHLIB_BUILD_EXAMPLES)
    add_executable(hello_tensor examples/hello_tensor.c)
    target_link_libraries(hello_tensor mathc)
target_include_directories(hello_tensor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

option(MATHLIB_BUILD_VALIDATION "Build validation tools" OFF)
if(MATHLIB_BUILD_VALIDATION)
    add_executable(ulp_validation validation/ulp_error.c)
    target_link_libraries(ulp_validation mathc)
target_include_directories(ulp_validation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

# Installation and Export (Enables find_package(mathlib) for downstream consumers)
include(GNUInstallDirs)
install(TARGETS mathc EXPORT mathlibTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/mathlib/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT mathlibTargets
    FILE mathlibTargets.cmake
    NAMESPACE mathlib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mathlib
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/mathlibConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/mathlibConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mathlib
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/mathlibConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mathlibConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/mathlibConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mathlib
)

# MATHLIB_V12A1_ORACLE_BUILD_FIX
# oracle_check is a first-class target, always built.
# It inherits sanitizer link flags from add_link_options() when
# MATHLIB_SANITIZERS=ON, eliminating the manual link lines that
# produced undefined __asan_* / __ubsan_* references.
add_executable(oracle_check tests/test_oracle.c)
target_compile_definitions(oracle_check PRIVATE MATHLIB_HAS_ORACLE_DATA)
target_include_directories(oracle_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(oracle_check mathc)
