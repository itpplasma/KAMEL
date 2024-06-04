set(KIM_setup setup/config_mod.f90 
              setup/setup_mod.f90
)

set(KIM_poisson poisson_equation/solve_poisson.f90
)

set(KIM_math math/suitesparse/umf4_f77zwrapper.c
             math/suitesparse/umf4_f77wrapper.c
             math/superlu/c_fortran_dgssv.c
             math/superlu/c_fortran_zgssv.c
             math/bessel.f90
             math/bessk.f90
             math/use_libcerf_mod.f90
             math/sparse_mod.f90
             math/plasma_disp_func.f90
)

set(KIM_debye debye_case/debye_kernel.f90
              debye_case/fill_spline_kernel_debye.f90
)

set(KIM_grid grid/generate_k_space_grid.f90
             grid/gengrid.f90
             grid/kr_space_adjustments.f90
             grid/spline_function.f90
             grid/kr_grid_mod.f90
             grid/grid_mod.f90
             grid/resonances_mod.f90
)

set(KIM_util util/findIndex.f90
             util/loading_bar.f90
             util/plag_coeff.f90
             util/unit_tests.f90
             util/dense_to_sparse.f90
             util/write_config_to_hdf5.f90
             util/constants_mod.f90
)

set(KIM_kernels kernels/integrands.f90
                #kernels/kernel_functions.f90
                kernels/cut_off_integration.f90
                kernels/kernel_mod.f90
)

set(KIM_background_equilibrium background_equilibrium/calculate_backs.f90
                               background_equilibrium/calculate_equil.f90
                               background_equilibrium/plasma_parameter.f90
                               background_equilibrium/back_quants_mod.f90
)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_library(KIM_lib STATIC "${KIM_setup}" 
                           "${KIM_poisson}" 
                           "${KIM_debye}" 
                           "${KIM_grid}" 
                           "${KIM_util}"
                           "${KIM_kernels}"
                           "${KIM_background_equilibrium}"
                           "${KIM_math}"
)

find_package(LAPACK REQUIRED)
if("${LAPACK_FOUND}")
    message(" ====== LAPACK was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${LAPACK_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${LAPACK_INCLUDE_DIRS}")
    message("${LAPACK_INCLUDE_DIRS}")
endif()

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
