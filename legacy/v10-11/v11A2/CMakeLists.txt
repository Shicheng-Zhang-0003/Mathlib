cmake_minimum_required(VERSION 3.10)
project(mathlib_v1_0 C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(MATHLIB_PROFILE "SCIENTIFIC" CACHE STRING "Mathlib Execution Profile")
if(MATHLIB_PROFILE STREQUAL "GRAPHICS")
    add_compile_definitions(MATHLIB_PROFILE_GRAPHICS)
elseif(MATHLIB_PROFILE STREQUAL "EMBEDDED")
    add_compile_definitions(MATHLIB_PROFILE_EMBEDDED)
endif()

option(MATHLIB_NATIVE "Enable -march=native (Warning: Breaks binary portability)" OFF)
if(MATHLIB_NATIVE)
    add_compile_options(-O3 -march=native -Wall -Wextra -fno-fast-math -ffp-contract=off)
else()
    add_compile_options(-O3 -Wall -Wextra -fno-fast-math -ffp-contract=off)
endif()


option(MATHLIB_HARDENED "Enable memory canaries and strict input sanitization" OFF)
if(MATHLIB_HARDENED)
    add_compile_definitions(MATHLIB_PROFILE_HARDENED)
    message(STATUS "Hardened perimeter enabled")
endif()
add_library(mathc STATIC src/mathlib.c src/cpu_dispatch.c)
target_include_directories(mathc PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/mathlib>
    $<INSTALL_INTERFACE:include/mathlib>
)

add_executable(test tests/test.c)
target_link_libraries(test mathc m)

add_executable(bench benchmarks/bench.c)
target_link_libraries(bench mathc m)

add_executable(fuzz_god_mode tests/fuzz_god_mode.c)
target_link_libraries(fuzz_god_mode mathc m)

add_executable(fuzz_boundary tests/fuzz_boundary_gauntlet.c)
target_link_libraries(fuzz_boundary mathc m)

include(GNUInstallDirs)
install(TARGETS mathc EXPORT mathlibTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/mathlib DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT mathlibTargets FILE mathlibTargets.cmake NAMESPACE mathlib:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mathlib)
option(MATHLIB_ENABLE_LEGACY "Enable legacy malloc-based matrix and FFT modules" OFF)
if(MATHLIB_ENABLE_LEGACY)
    message(STATUS "Legacy modules enabled")
    endif()


# Modular CI/CD Tests
add_executable(test_core tests/test_core.c)
target_link_libraries(test_core mathc m)

add_executable(test_trig tests/test_trig.c)
target_link_libraries(test_trig mathc m)

add_executable(test_linalg tests/test_linalg.c)
target_link_libraries(test_linalg mathc m)

add_executable(test_dsp tests/test_dsp.c)
target_link_libraries(test_dsp mathc m)

option(MATHLIB_SANITIZERS "Enable Address and Undefined Behavior Sanitizers" OFF)
if(MATHLIB_SANITIZERS)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
    message(STATUS "Sanitizers enabled")
endif()

add_executable(ulp_validation validation/ulp_error.c)
target_link_libraries(ulp_validation mathc m)

add_executable(hello_tensor examples/hello_tensor.c)
target_link_libraries(hello_tensor mathc m)
