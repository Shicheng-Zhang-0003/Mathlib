cmake_minimum_required(VERSION 3.10)
project(mathlib_v11S VERSION 11.0.0 LANGUAGES C)

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
add_compile_options(-O3 -Wall -Wextra -Wconversion -Wshadow -Wpedantic -Werror -fno-fast-math -ffp-contract=off)

option(MATHLIB_NATIVE "Enable -march=native (Warning: Breaks binary portability)" OFF)
if(MATHLIB_NATIVE)
    add_compile_options(-march=native)
endif()

# Sanitizers (Fixed: Added to link options as well)
option(MATHLIB_SANITIZERS "Enable Address and Undefined Behavior Sanitizers" OFF)
if(MATHLIB_SANITIZERS)
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
)

add_library(mathc STATIC ${MATHLIB_SOURCES})
target_include_directories(mathc PRIVATE src)
target_include_directories(mathc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/mathlib>
    $<INSTALL_INTERFACE:include>
)

# Core Test Executables
add_executable(test tests/test.c)
target_link_libraries(test mathc m)
target_include_directories(test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(bench benchmarks/bench.c)
target_link_libraries(bench mathc m)
target_include_directories(bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(fuzz_god_mode tests/fuzz_god_mode.c)
target_link_libraries(fuzz_god_mode mathc m)
target_include_directories(fuzz_god_mode PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(fuzz_boundary tests/fuzz_boundary_gauntlet.c)
target_link_libraries(fuzz_boundary mathc m)
target_include_directories(fuzz_boundary PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Modular CI/CD Tests
add_executable(test_core tests/test_core.c)
target_link_libraries(test_core mathc m)
target_include_directories(test_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_trig tests/test_trig.c)
target_link_libraries(test_trig mathc m)
target_include_directories(test_trig PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_linalg tests/test_linalg.c)
target_link_libraries(test_linalg mathc m)
target_include_directories(test_linalg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_dsp tests/test_dsp.c)
target_link_libraries(test_dsp mathc m)
target_include_directories(test_dsp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Examples & Validation (Hidden behind options to prevent CI bloat)
option(MATHLIB_BUILD_EXAMPLES "Build example programs" OFF)
if(MATHLIB_BUILD_EXAMPLES)
    add_executable(hello_tensor examples/hello_tensor.c)
    target_link_libraries(hello_tensor mathc m)
target_include_directories(hello_tensor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

option(MATHLIB_BUILD_VALIDATION "Build validation tools" OFF)
if(MATHLIB_BUILD_VALIDATION)
    add_executable(ulp_validation validation/ulp_error.c)
    target_link_libraries(ulp_validation mathc m)
target_include_directories(ulp_validation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

# Installation and Export (Enables find_package(mathlib) for downstream consumers)
include(GNUInstallDirs)
install(TARGETS mathc EXPORT mathlibTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/mathlib DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
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
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mathlibConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mathlib
)

option(MATHLIB_ORACLE_TESTS "Enable mpmath oracle validation" OFF)
if(MATHLIB_ORACLE_TESTS)
    add_executable(test_oracle tests/test_oracle.c)
    target_compile_definitions(test_oracle PRIVATE MATHLIB_HAS_ORACLE_DATA)
    target_link_libraries(test_oracle mathc m)
endif()
