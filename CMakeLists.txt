cmake_minimum_required(VERSION 3.24)

set(CMAKE_Fortran_COMPILER "mpif90")

project(KAMEL LANGUAGES C CXX Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion")

enable_testing()

include(${CMAKE_SOURCE_DIR}/cmake/Dependencies.cmake)

set(CMAKE_BUILD_TYPE Release)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Unified build: add subprojects for KiLCA, KIM, and QL-Balance
add_subdirectory(common/math)
add_subdirectory(KiLCA)
add_subdirectory(KIM)
add_subdirectory(QL-Balance)
