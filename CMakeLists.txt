cmake_minimum_required(VERSION 3.18)

project(KIM VERSION 0.1
            DESCRIPTION "KiLCA Integral Model"
            LANGUAGES Fortran C CXX)
enable_language(Fortran)

add_subdirectory(src)
add_subdirectory(src/math/slatec)
add_subdirectory(src/math/lapack-3.11.0)
add_subdirectory(src/math/libcerf-main)
add_subdirectory(src/math/ddeabm)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set (CMAKE_Fortran_FLAGS "-O2 -mtune=generic")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_FLAGS "")
endif()

