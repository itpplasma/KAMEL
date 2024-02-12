set(KIM_mods modules/config_mod.f90 
             modules/constants_mod.f90
             modules/setup_mod.f90
             modules/grid_mod.f90
             modules/back_quants_mod.f90
             modules/kernel_mod.f90
             modules/resonances_mod.f90)

set(KIM_subs plasma_parameter.f90
             plag_coeff.f90
             #kernel_phi.f90
             #kernel_B.f90
             calculate_backs.f90
             calculate_equil.f90
             gengrid.f90
             plasma_disp_func.f90
             spline_function.f90
             #basis_transform_kernel.f90
             solve_poisson.f90
             #integrands.f90
             loading_bar.f90
             #integrate_krp_kr_rg.f90
             #integration.f90
             #fft_kernel.f90
             cut_off_integration.f90
             kernel_functions.f90
             math/suitesparse/umf4_f77zwrapper.c
             math/suitesparse/umf4_f77wrapper.c
             math/superlu/c_fortran_dgssv.c
             math/superlu/c_fortran_zgssv.c
             math/bessel.f90
             math/bessk.f90
             math/use_libcerf_mod.f90
             math/sparse_mod.f90
             math/adaptive_int.f90
             hdf5_tools/hdf5_tools.f90)


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_library(KIM_lib "${KIM_mods}" "${KIM_subs}")
set_target_properties(KIM_lib PROPERTIES
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/")

find_package(SuiteSparse REQUIRED)
if("${SuiteSparse_FOUND}")
    message(" ====== SuiteSparse was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${SUITESPARSE_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${SUITESPARSE_INCLUDE_DIRS}")
    message("${SUITESPARSE_INCLUDE_DIRS}")
endif()

find_package(Umfpack REQUIRED)
if("${UMFPACK_FOUND}")
    message(" ====== Umfpack was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${UMFPACK_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${UMFPACK_INCLUDES}")
    message("${UMFPACK_INCLUDES}")
endif()

#find_package(SUNDIALS REQUIRED)
#if("${SUNDIALS_FOUND}")
#    message(" ====== Sundials was found! =====")
#    target_link_libraries(KIM_lib PUBLIC "${SUNDIALS_LIBRARIES}")
#    target_include_directories(KIM_lib PUBLIC "${SUNDIALS_INCLUDE_DIRS}")
#    message("${SUNDIALS_INCLUDE_DIRS}")
#endif()

find_package(SuperLU REQUIRED)
if("${SUPERLU_FOUND}")
    message(" ====== SuperLU was found! =====")
    target_link_libraries(KIM_lib PUBLIC "${SUPERLU_LIBRARIES}")
    message("${SUPERLU_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${SUPERLU_INCLUDES}")
    message("${SUPERLU_INCLUDES}")
endif()

find_library(UMFPACK_LIBRARY umfpack)

find_package(OpenMP_Fortran REQUIRED)
IF(OPENMP_FORTRAN_FOUND)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()

#find_package(OpenMP REQUIRED)

#if(OpenMP_C_FOUND)
#    message(" ====== OPENMP_C was found! =====")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    target_include_directories(KIM_lib PUBLIC "${OpenMP_C_INCLUDE_DIRS}")
#endif()

#find_package(FFTW3 COMPONENTS fftw3 REQUIRED)
#IF("${FFTW3_FOUND}")
#    message(" ====== FFTW3 was found! =====")
#    target_link_libraries(KIM_lib PUBLIC "${FFTW3_fftw3_LIBRARY}")
#    message("${FFTW3_fftw3_LIBRARY}")
#    target_include_directories(KIM_lib PUBLIC "${FFTW3_INCLUDE_DIRS}")
#    message("${FFTW3_INCLUDE_DIRS}")
#endif()

file(
    DOWNLOAD
    https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.38.3/CPM.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
    EXPECTED_HASH SHA256=cc155ce02e7945e7b8967ddfaff0b050e958a723ef7aad3766d368940cb15494
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

#CPMAddPackage(
#    NAME Finufft
#    GIT_REPOSITORY https://github.com/flatironinstitute/finufft.git
#    GIT_TAG master
#    GIT_SHALLOW Yes
#    GIT_PROGRESS Yes
#    EXCLUDE_FROM_ALL Yes
#    SYSTEM
#)
#target_link_libraries(KIM_lib PUBLIC finufft_static)

# or for shared linking
#target_link_library(KIM_lib [PUBLIC|PRIVATE|INTERFACE] finufft)

add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM_exe
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/")
target_link_libraries(KIM_exe KIM_lib 
                              bessel 
                              lapack 
                              cerf 
                              ddeabm)
                              #finufft_static)
