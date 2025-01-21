cmake_minimum_required(VERSION 3.24)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(KIM VERSION 0.1
            DESCRIPTION "KiLCA Integral Model"
            LANGUAGES Fortran C CXX)
enable_language(Fortran)

add_subdirectory(src)
add_subdirectory(src/math/slatec)
#add_subdirectory(src/math/lapack-3.11.0)
add_subdirectory(src/math/libcerf-main)
add_subdirectory(src/math/ddeabm)
#add_subdirectory(src/hdf5_tools)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(FetchZeal)

#find_package(OpenMP_Fortran)
#find_package(OpenMP REQUIRED Fortran)
#message("Test this 2")
#IF(OPENMP_FORTRAN_FOUND)
#message("Test this ")
    #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#ENDIF()

message("${CMAKE_Fortran_COMPILER}")
message("${CMAKE_Fortran_FLAGS}")
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message("set gfortran compiler, will use flags")
    set(CMAKE_Fortran_FLAGS "-O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_FLAGS "")
endif()

