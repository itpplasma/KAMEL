cmake_minimum_required(VERSION 3.24)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_Fortran_COMPILER "mpif90")
set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS}"
    "-O3"
    "-Wall"
    "-fbacktrace"
    "-fcheck=all"
    "-ffree-line-length-none"
    "-fopenmp"
    "-g"
    "-mtune=generic"
)
string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

project(KIM VERSION 0.1
            DESCRIPTION "KiLCA Integral Model"
            LANGUAGES Fortran C CXX)
enable_language(Fortran)
enable_testing()

include(FetchZeal)
add_subdirectory(src)
add_subdirectory(src/math/slatec)
add_subdirectory(src/math/libcerf-main)
add_subdirectory(src/math/ddeabm)
add_subdirectory(src/tests)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
