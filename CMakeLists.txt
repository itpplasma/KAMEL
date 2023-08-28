cmake_minimum_required(VERSION 3.18)

project(KIM VERSION 0.1
            DESCRIPTION "KiLCA Integral Model"
            LANGUAGES Fortran C CXX)
enable_language(Fortran)


set(exe "${CMAKE_SOURCE_DIR}/exe/")

add_subdirectory(src)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set (CMAKE_Fortran_FLAGS "-O2 -mtune=generic")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set (CMAKE_Fortran_FLAGS "")
endif()

