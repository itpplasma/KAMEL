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

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#find_package(OpenMP_Fortran)
find_package(OpenMP_Fortran REQUIRED)
IF(OPENMP_FORTRAN_FOUND)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()

message("${CMAKE_Fortran_COMPILER}")
message("${CMAKE_Fortran_FLAGS}")
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message("set gfortran compiler, will use flags")
    set(CMAKE_Fortran_FLAGS "-O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_FLAGS "")
endif()

