set(KIM_mods modules/config_mod.f90 
             modules/constants_mod.f90
             modules/setup_mod.f90
             modules/grid_mod.f90
             modules/back_quants_mod.f90
             modules/kernel_mod.f90
             modules/resonances_mod.f90)

set(KIM_subs plasma_parameter.f90
             plag_coeff.f90
             calculate_backs.f90
             calculate_equil.f90
             gengrid.f90
             plasma_disp_func.f90
             spline_function.f90
             solve_poisson.f90
             loading_bar.f90
             cut_off_integration.f90
             kernel_functions.f90
             math/suitesparse/umf4_f77zwrapper.c
             math/suitesparse/umf4_f77wrapper.c
             math/superlu/c_fortran_dgssv.c
             math/superlu/c_fortran_zgssv.c
             math/bessel.f90
             math/bessk.f90
             math/use_libcerf_mod.f90
             math/sparse_mod.f90)


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_library(KIM_lib STATIC "${KIM_mods}" "${KIM_subs}")
#set_target_properties(KIM_lib PROPERTIES
#                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/")

find_package(LAPACK REQUIRED)
if("${LAPACK_FOUND}")
    message(" ====== LAPACK was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${LAPACK_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${LAPACK_INCLUDE_DIRS}")
    message("${LAPACK_INCLUDE_DIRS}")
endif()

#find_package(SLATEC REQUIRED)
#if("${SLATEC_FOUND}")
#    message(" ====== SLATEC was found! ======")
#    target_link_libraries(KIM_lib PUBLIC "${SLATEC_LIBRARIES}")
#    target_include_directories(KIM_lib PUBLIC "${SLATEC_INCLUDE_DIRS}")
#    message("${SLATEC_INCLUDE_DIRS}")
#endif()

find_package(SuiteSparse REQUIRED)
if("${SuiteSparse_FOUND}")
    message(" ====== SuiteSparse was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${SuiteSparse_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${SuiteSparse_INCLUDE_DIRS}")
    message("${SUITESPARSE_INCLUDE_DIRS}")
endif()

find_package(Umfpack REQUIRED)
if("${UMFPACK_FOUND}")
    message(" ====== Umfpack was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${UMFPACK_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${UMFPACK_INCLUDES}")
    message("${UMFPACK_INCLUDES}")
endif()

find_library(UMFPACK_LIBRARY umfpack)

find_package(SuperLU REQUIRED)
if("${SUPERLU_FOUND}")
    message(" ====== SuperLU was found! =====")
    target_link_libraries(KIM_lib PUBLIC "${SUPERLU_LIBRARIES}")
    message("${SUPERLU_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${SUPERLU_INCLUDES}")
    message("${SUPERLU_INCLUDES}")
endif()


find_package(OpenMP_Fortran REQUIRED)
IF(OPENMP_FORTRAN_FOUND)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()

file(
    DOWNLOAD
    https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.38.3/CPM.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
    EXPECTED_HASH SHA256=cc155ce02e7945e7b8967ddfaff0b050e958a723ef7aad3766d368940cb15494
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)


add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM_exe
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/")
target_link_libraries(KIM_exe KIM_lib 
                              lapack 
                              cerf 
                              slatec 
                              ddeabm)
                              #finufft_static)
