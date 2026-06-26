cmake_minimum_required(VERSION 3.10)
project(mathlib_v9 C)

set(CMAKE_C_STANDARD 99)
set(MATHLIB_PROFILE "SCIENTIFIC" CACHE STRING "Mathlib Execution Profile (GRAPHICS, SCIENTIFIC, EMBEDDED)")

if(MATHLIB_PROFILE STREQUAL "GRAPHICS")
    add_compile_definitions(MATHLIB_PROFILE_GRAPHICS)
    message(STATUS "Building with GRAPHICS profile (Speed Optimized)")
elseif(MATHLIB_PROFILE STREQUAL "EMBEDDED")
    add_compile_definitions(MATHLIB_PROFILE_EMBEDDED)
    message(STATUS "Building with EMBEDDED profile (No FPU / CORDIC)")
else()
    message(STATUS "Building with SCIENTIFIC profile (Strict IEEE-754)")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -mavx -Wall -Wextra -fgnu89-inline")

add_library(mathc STATIC mathlib.c)
target_include_directories(mathc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test test.c)
target_link_libraries(test mathc m)

add_executable(bench bench.c)
target_link_libraries(bench mathc m)
