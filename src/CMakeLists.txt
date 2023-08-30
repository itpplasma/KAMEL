set(KIM_mods plas_parameter_mod.f90 
             config_mod.f90 
             constants_mod.f90
             setup_mod.f90
             grid_mod.f90
             back_quants_mod.f90)

set(KIM_subs read_profiles.f90
             plag_coeff.f90
             kernel_rho_phi.f90
             calc_backs.f90
             math/bessel.f90
             math/bessk.f90
             math/use_libcerf_mod.f90)

add_library(KIM_lib "${KIM_mods}" "${KIM_subs}")
set_target_properties(KIM_lib PROPERTIES
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/")

add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM_exe
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/")
target_link_libraries(KIM_exe PRIVATE KIM_lib bessel lapack cerf)