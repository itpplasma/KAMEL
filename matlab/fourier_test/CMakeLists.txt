cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(Fouriermodes)

enable_language(Fortran)

if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
    include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()

#find_program(BASH_PROGRAM bash)
#find_package(MPI REQUIRED)

add_library(mc_efit SHARED
 src/field_divB0.f90
 src/field_eq_mod.f90
 src/bdivfree.f90
 src/bdivfree_mod.f90
 src/spline5_RZ.f90
 src/rk4dc.f90
 src/mag.f90
 src/rhs_flt.f90
 src/input_files.f90
 src/period_mod.f90
 src/field_c_mod.f90
 src/field_mod.f90
 src/amn_mod.f90
 src/theta_rz_mod.f90
 src/extract_fluxcoord_mod.f90
 src/inthecore_mod.f90
 src/fouriermodes.f90)

add_executable(fouriermodes.x src/fouriermodes.f90)
target_link_libraries(fouriermodes.x mc_efit)