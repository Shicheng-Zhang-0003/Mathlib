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

add_compile_options(-O3 -march=native -Wall -Wextra)

add_library(mathc STATIC src/mathlib.c)
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