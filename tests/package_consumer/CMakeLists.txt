cmake_minimum_required(VERSION 3.14)
project(mathlib_package_consumer C)

find_package(mathlib 12.3 CONFIG REQUIRED)

add_executable(mathlib_package_consumer main.c)
target_link_libraries(mathlib_package_consumer PRIVATE mathlib::mathc)
