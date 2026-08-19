# OpenMP runtime (GOMP_critical_* from libneo's OMP critical sections) is
# carried transitively by kilca_lib's PUBLIC link interface, so no per-target
# OpenMP linkage is needed here.

# Optional: build QUADPACK/RKF45 integration comparison test
add_executable(test_integration_methods ${CMAKE_SOURCE_DIR}/KIM/tests/test_integration_methods.f90)
set_target_properties(test_integration_methods PROPERTIES
                        OUTPUT_NAME test_integration_methods.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_integration_methods KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_integration_methods
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_integration_methods.x)


add_executable(test_root_finding ${CMAKE_SOURCE_DIR}/KIM/tests/test_root_finding.f90)
set_target_properties(test_root_finding PROPERTIES
                        OUTPUT_NAME test_root_finding.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_root_finding KIM_lib)
add_test(NAME test_root_finding
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_root_finding.x)

add_executable(test_region_roots_vs_muller ${CMAKE_SOURCE_DIR}/KIM/tests/test_region_roots_vs_muller.f90)
set_target_properties(test_region_roots_vs_muller PROPERTIES
                        OUTPUT_NAME test_region_roots_vs_muller.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_region_roots_vs_muller KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_region_roots_vs_muller
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_region_roots_vs_muller.x)

add_executable(test_equil_ode_ddeabm ${CMAKE_SOURCE_DIR}/KIM/tests/test_equil_ode_ddeabm.f90)
set_target_properties(test_equil_ode_ddeabm PROPERTIES
                        OUTPUT_NAME test_equil_ode_ddeabm.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_equil_ode_ddeabm KIM_lib)
add_test(NAME test_equil_ode_ddeabm
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_equil_ode_ddeabm.x)

add_executable(test_wavenumber_geometry
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_wavenumber_geometry.f90)
set_target_properties(test_wavenumber_geometry PROPERTIES
                      OUTPUT_NAME test_wavenumber_geometry.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_wavenumber_geometry KIM_lib)
add_test(NAME test_wavenumber_geometry
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_wavenumber_geometry.x)

add_executable(test_profile_input ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input.f90)
set_target_properties(test_profile_input PROPERTIES
                        OUTPUT_NAME test_profile_input.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
# Links the full math set: using profile_input_m pulls KIM objects whose
# transitive deps (W2_arr) need the hypergeometric routines from KiLCA.
target_link_libraries(test_profile_input KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_profile_input
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input.x)

add_executable(test_profile_input_integration ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input_integration.f90)
set_target_properties(test_profile_input_integration PROPERTIES
                        OUTPUT_NAME test_profile_input_integration.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_profile_input_integration KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_profile_input_integration
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input_integration.x)

add_executable(test_kim_diagnostics ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_diagnostics.f90)
set_target_properties(test_kim_diagnostics PROPERTIES
                        OUTPUT_NAME test_kim_diagnostics.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_diagnostics KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_kim_diagnostics
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_diagnostics.x)

add_executable(test_periodic_ion_tensor
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_ion_tensor.f90)
set_target_properties(test_periodic_ion_tensor PROPERTIES
                      OUTPUT_NAME test_periodic_ion_tensor.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_ion_tensor KIM_lib kilca_lib lapack cerf fortnum fortnum_amos_compat)
add_test(NAME test_periodic_ion_tensor
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_ion_tensor.x)

add_executable(test_quasilinear_integral_algebra
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_quasilinear_integral_algebra.f90)
file(SHA256
  ${CMAKE_SOURCE_DIR}/KIM/mathematica/verify_quasilinear_integral_transport.wl
  QL_INTEGRAL_ORACLE_SHA256)
foreach(oracle_fixture IN ITEMS
    quasilinear_transverse_oracle.dat
    quasilinear_channel_oracle.dat)
  file(STRINGS ${CMAKE_SOURCE_DIR}/KIM/tests/test_data/${oracle_fixture}
       oracle_header LIMIT_COUNT 2)
  list(GET oracle_header 1 oracle_hash_header)
  if(NOT oracle_hash_header STREQUAL
      "# generator_sha256 ${QL_INTEGRAL_ORACLE_SHA256}")
    message(FATAL_ERROR
      "Stale ${oracle_fixture}; regenerate it with "
      "KIM/mathematica/verify_quasilinear_integral_transport.wl")
  endif()
endforeach()
configure_file(
  ${CMAKE_SOURCE_DIR}/KIM/tests/test_data/quasilinear_transverse_oracle.dat
  ${CMAKE_BINARY_DIR}/tests/quasilinear_transverse_oracle.dat COPYONLY)
configure_file(
  ${CMAKE_SOURCE_DIR}/KIM/tests/test_data/quasilinear_channel_oracle.dat
  ${CMAKE_BINARY_DIR}/tests/quasilinear_channel_oracle.dat COPYONLY)
set_target_properties(test_quasilinear_integral_algebra PROPERTIES
                      OUTPUT_NAME test_quasilinear_integral_algebra.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_quasilinear_integral_algebra
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_quasilinear_integral_algebra
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_quasilinear_integral_algebra.x)
set_tests_properties(test_quasilinear_integral_algebra PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

add_executable(test_collisionless_ion_kernel ${CMAKE_SOURCE_DIR}/KIM/tests/test_collisionless_ion_kernel.f90)
set_target_properties(test_collisionless_ion_kernel PROPERTIES
                        OUTPUT_NAME test_collisionless_ion_kernel.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_collisionless_ion_kernel KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_collisionless_ion_kernel
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_collisionless_ion_kernel.x)

# m_phi=0 collisionless ion charge/current kernels in the periodic Fourier basis.
add_executable(test_collisionless_fourier_kernel ${CMAKE_SOURCE_DIR}/KIM/tests/test_collisionless_fourier_kernel.f90)
set_target_properties(test_collisionless_fourier_kernel PROPERTIES
                        OUTPUT_NAME test_collisionless_fourier_kernel.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_collisionless_fourier_kernel KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_collisionless_fourier_kernel
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_collisionless_fourier_kernel.x)

add_executable(test_collisionless_ion_assembly ${CMAKE_SOURCE_DIR}/KIM/tests/test_collisionless_ion_assembly.f90)
set_target_properties(test_collisionless_ion_assembly PROPERTIES
                        OUTPUT_NAME test_collisionless_ion_assembly.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_collisionless_ion_assembly KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_collisionless_ion_assembly
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_collisionless_ion_assembly.x)
set_tests_properties(test_collisionless_ion_assembly PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

add_executable(test_susceptibility_collisionless_limit
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_susceptibility_collisionless_limit.f90)
set_target_properties(test_susceptibility_collisionless_limit PROPERTIES
                      OUTPUT_NAME test_susceptibility_collisionless_limit.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_susceptibility_collisionless_limit
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_susceptibility_collisionless_limit
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_susceptibility_collisionless_limit.x)

add_executable(test_ifunc_conservation_models
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_ifunc_conservation_models.f90)
set_target_properties(test_ifunc_conservation_models PROPERTIES
                      OUTPUT_NAME test_ifunc_conservation_models.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ifunc_conservation_models
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_ifunc_conservation_models
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ifunc_conservation_models.x)

add_executable(test_ifunc_model_config
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_ifunc_model_config.f90)
set_target_properties(test_ifunc_model_config PROPERTIES
                      OUTPUT_NAME test_ifunc_model_config.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ifunc_model_config
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_ifunc_model_config
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ifunc_model_config.x)

add_executable(test_ion_temperature_gradient_model
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_ion_temperature_gradient_model.f90)
set_target_properties(test_ion_temperature_gradient_model PROPERTIES
                      OUTPUT_NAME test_ion_temperature_gradient_model.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ion_temperature_gradient_model
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_ion_temperature_gradient_model
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ion_temperature_gradient_model.x)

add_executable(test_flr2_fourier_kernel ${CMAKE_SOURCE_DIR}/KIM/tests/test_flr2_fourier_kernel.f90)
set_target_properties(test_flr2_fourier_kernel PROPERTIES
                        OUTPUT_NAME test_flr2_fourier_kernel.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_flr2_fourier_kernel KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_flr2_fourier_kernel
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_flr2_fourier_kernel.x)

add_executable(test_flr2_tridiagonal
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_flr2_tridiagonal.f90)
set_target_properties(test_flr2_tridiagonal PROPERTIES
                      OUTPUT_NAME test_flr2_tridiagonal.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_flr2_tridiagonal KIM_lib kilca_lib lapack cerf
                      fortnum_amos_compat)
add_test(NAME test_flr2_tridiagonal
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_flr2_tridiagonal.x)

add_executable(test_flr2_response
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_flr2_response.f90)
set_target_properties(test_flr2_response PROPERTIES
                      OUTPUT_NAME test_flr2_response.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_flr2_response KIM_lib kilca_lib lapack cerf
                      fortnum_amos_compat)
add_test(NAME test_flr2_response
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_flr2_response.x)

add_executable(test_flr2_factory
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_flr2_factory.f90)
set_target_properties(test_flr2_factory PROPERTIES
                      OUTPUT_NAME test_flr2_factory.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_flr2_factory KIM_lib kilca_lib lapack cerf
                      fortnum_amos_compat)
add_test(NAME test_flr2_factory
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_flr2_factory.x)

add_executable(test_kim_solver_flr2
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver_flr2.f90)
set_target_properties(test_kim_solver_flr2 PROPERTIES
                      OUTPUT_NAME test_kim_solver_flr2.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver_flr2 KIM_lib kilca_lib lapack cerf
                      fortnum_amos_compat)
configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/KIM_config_flr2_small.nml
               ${CMAKE_BINARY_DIR}/tests/KIM_config_flr2_small.nml COPYONLY)
add_test(NAME test_kim_solver_flr2
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver_flr2.x)
set_tests_properties(test_kim_solver_flr2 PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

add_executable(test_ampere_matrices ${CMAKE_SOURCE_DIR}/KIM/tests/test_ampere_matrices.f90)
set_target_properties(test_ampere_matrices PROPERTIES
                        OUTPUT_NAME test_ampere_matrices.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ampere_matrices KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_ampere_matrices
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ampere_matrices.x)

# Pure-utility unit tests (Lagrange coefficients, binary search, index lookup,
# equidistant grid generation). Full math link set: pulling these KIM objects
# drags transitive deps (W2_arr) needing hypergeometric routines from KiLCA.
add_executable(test_plag_coeff ${CMAKE_SOURCE_DIR}/KIM/tests/test_plag_coeff.f90)
set_target_properties(test_plag_coeff PROPERTIES
                        OUTPUT_NAME test_plag_coeff.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_plag_coeff KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_plag_coeff
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_plag_coeff.x)

add_executable(test_findIndex ${CMAKE_SOURCE_DIR}/KIM/tests/test_findIndex.f90)
set_target_properties(test_findIndex PROPERTIES
                        OUTPUT_NAME test_findIndex.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_findIndex KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_findIndex
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_findIndex.x)

add_executable(test_grid_equidistant ${CMAKE_SOURCE_DIR}/KIM/tests/test_grid_equidistant.f90)
set_target_properties(test_grid_equidistant PROPERTIES
                        OUTPUT_NAME test_grid_equidistant.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_grid_equidistant KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_grid_equidistant
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_grid_equidistant.x)

# Derived Krook argument z0 recompute + collision-frequency scale.
add_executable(test_collision_scale_z0 ${CMAKE_SOURCE_DIR}/KIM/tests/test_collision_scale_z0.f90)
set_target_properties(test_collision_scale_z0 PROPERTIES
                        OUTPUT_NAME test_collision_scale_z0.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_collision_scale_z0 KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_collision_scale_z0
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_collision_scale_z0.x)

add_executable(test_fp_harmonic_prefactors
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_fp_harmonic_prefactors.f90)
set_target_properties(test_fp_harmonic_prefactors PROPERTIES
                      OUTPUT_NAME test_fp_harmonic_prefactors.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_fp_harmonic_prefactors
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_fp_harmonic_prefactors
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_fp_harmonic_prefactors.x)

add_executable(test_fp_kernel_prefactors
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_fp_kernel_prefactors.f90)
set_target_properties(test_fp_kernel_prefactors PROPERTIES
                      OUTPUT_NAME test_fp_kernel_prefactors.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_fp_kernel_prefactors
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_fp_kernel_prefactors
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_fp_kernel_prefactors.x)

add_executable(test_ion_flr_scale_factor
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_ion_flr_scale_factor.f90)
set_target_properties(test_ion_flr_scale_factor PROPERTIES
                      OUTPUT_NAME test_ion_flr_scale_factor.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ion_flr_scale_factor
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_ion_flr_scale_factor
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ion_flr_scale_factor.x)

# Regression test for the vendored periodization utility.
add_executable(test_periodization ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodization.f90)
set_target_properties(test_periodization PROPERTIES
                        OUTPUT_NAME test_periodization.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodization KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodization
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodization.x)

# Periodic sampling of KIM background primitives (aperfuns callback + wrapper).
add_executable(test_periodic_background ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_background.f90)
set_target_properties(test_periodic_background PROPERTIES
                        OUTPUT_NAME test_periodic_background.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_background KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodic_background
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_background.x)

# Build a periodic background plasma on a redirected window grid (Phase-2.3).
add_executable(test_periodic_background_build ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_background_build.f90)
set_target_properties(test_periodic_background_build PROPERTIES
                        OUTPUT_NAME test_periodic_background_build.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_background_build KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodic_background_build
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_background_build.x)
set_tests_properties(test_periodic_background_build PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Assemble the dense periodic Fourier matrices K_{m,m'} (Phase-2.4).
add_executable(test_periodic_assembly ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_assembly.f90)
set_target_properties(test_periodic_assembly PROPERTIES
                        OUTPUT_NAME test_periodic_assembly.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_assembly KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodic_assembly
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_assembly.x)
set_tests_properties(test_periodic_assembly PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Dense zgesv solve + inverse DFT for the periodic potential (Phase-2.5).
add_executable(test_periodic_solve ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_solve.f90)
set_target_properties(test_periodic_solve PROPERTIES
                        OUTPUT_NAME test_periodic_solve.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_solve KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodic_solve
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_solve.x)
set_tests_properties(test_periodic_solve PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Lifecycle/contract tests for the kim_solver_m API seam.
add_executable(test_kim_solver ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver.f90)
set_target_properties(test_kim_solver PROPERTIES
                        OUTPUT_NAME test_kim_solver.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_kim_solver
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver.x)

# End-to-end keystone: a real electromagnetic solve through the kim_solver seam.
add_executable(test_kim_solver_em ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver_em.f90)
set_target_properties(test_kim_solver_em PROPERTIES
                        OUTPUT_NAME test_kim_solver_em.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver_em KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/KIM_config_em_small.nml
               ${CMAKE_BINARY_DIR}/tests/KIM_config_em_small.nml COPYONLY)
add_test(NAME test_kim_solver_em
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver_em.x)
set_tests_properties(test_kim_solver_em PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# End-to-end: the electrostatic_periodic run-type through the factory (Phase-2.6).
add_executable(test_kim_solver_periodic ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver_periodic.f90)
set_target_properties(test_kim_solver_periodic PROPERTIES
                        OUTPUT_NAME test_kim_solver_periodic.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver_periodic KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_kim_solver_periodic
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver_periodic.x)
set_tests_properties(test_kim_solver_periodic PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# End-to-end public result contract for the periodic KIM local response.
add_executable(test_kim_solver_periodic_contract
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver_periodic_contract.f90)
set_target_properties(test_kim_solver_periodic_contract PROPERTIES
                      OUTPUT_NAME test_kim_solver_periodic_contract.x
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver_periodic_contract
                      KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/KIM_config_em_small.nml
               ${CMAKE_BINARY_DIR}/tests/KIM_config_em_small.nml COPYONLY)
add_test(NAME test_kim_solver_periodic_contract
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver_periodic_contract.x)
set_tests_properties(test_kim_solver_periodic_contract PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Phase-2.7 cross-check: periodic vs global electrostatic delta-Phi on the same
# (m,n) = (-6,2) resonant case. Diagnostic (approx-vs-approx): logs the relative
# difference; soft gate (only on missing/all-zero/NaN).
add_executable(test_periodic_vs_global ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_vs_global.f90)
set_target_properties(test_periodic_vs_global PROPERTIES
                        OUTPUT_NAME test_periodic_vs_global.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_vs_global KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodic_vs_global
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_vs_global.x)
set_tests_properties(test_periodic_vs_global PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Phase-3.1: reusable periodic core compute_periodic_delta_phi evaluated on a
# FIXED diagnostic grid (basis for the convergence harness P3.2-3.4).
add_executable(test_periodic_convergence ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_convergence.f90)
set_target_properties(test_periodic_convergence PROPERTIES
                        OUTPUT_NAME test_periodic_convergence.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_convergence KIM_lib kilca_lib lapack cerf fortnum_amos_compat)
add_test(NAME test_periodic_convergence
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_convergence.x)
set_tests_properties(test_periodic_convergence PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)
