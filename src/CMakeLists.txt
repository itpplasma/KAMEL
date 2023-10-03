set(KIM_mods modules/plas_parameter_mod.f90 
             modules/config_mod.f90 
             modules/constants_mod.f90
             modules/setup_mod.f90
             modules/grid_mod.f90
             modules/back_quants_mod.f90
             modules/kernel_mod.f90
             modules/equil_mod.f90
             modules/resonances_mod.f90)

set(KIM_subs read_profiles.f90
             plag_coeff.f90
             kernel_phi.f90
             kernel_B.f90
             calculate_backs.f90
             calculate_equil.f90
             gengrid.f90
             plasma_disp_func.f90
             spline_function.f90
             basis_transform_kernel.f90
             solve_poisson.f90
             math/suitesparse/umf4_f77zwrapper.c
             math/suitesparse/umf4_f77wrapper.c
             math/superlu/c_fortran_dgssv.c
             math/superlu/c_fortran_zgssv.c
             math/bessel.f90
             math/bessk.f90
             math/use_libcerf_mod.f90
             math/sparse_mod.f90
             math/adaptive_int.f90)


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


add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM_exe
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/")
target_link_libraries(KIM_exe KIM_lib 
                              bessel 
                              lapack 
                              cerf 
                              ddeabm)
