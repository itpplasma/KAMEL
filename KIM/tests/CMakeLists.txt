# OpenMP runtime (GOMP_critical_* from libneo's OMP critical sections) is
# carried transitively by kilca_lib's PUBLIC link interface, so no per-target
# OpenMP linkage is needed here.

# Optional: build QUADPACK/RKF45 integration comparison test
add_executable(test_integration_methods ${CMAKE_SOURCE_DIR}/KIM/tests/test_integration_methods.f90)
set_target_properties(test_integration_methods PROPERTIES
                        OUTPUT_NAME test_integration_methods.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_integration_methods KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_integration_methods
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_integration_methods.x)

# QUADPACK DQAGI sanity test (infinite integral)
add_executable(test_quadpack_qagi ${CMAKE_SOURCE_DIR}/KIM/tests/test_quadpack_qagi.f90)
set_target_properties(test_quadpack_qagi PROPERTIES
                        OUTPUT_NAME test_quadpack_qagi.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_quadpack_qagi KIM_lib lapack cerf ddeabm slatec)
add_test(NAME test_quadpack_qagi
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_quadpack_qagi.x)


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
target_link_libraries(test_region_roots_vs_muller KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_region_roots_vs_muller
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_region_roots_vs_muller.x)

add_executable(test_profile_input ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input.f90)
set_target_properties(test_profile_input PROPERTIES
                        OUTPUT_NAME test_profile_input.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
# Links the full math set: using profile_input_m pulls KIM objects whose
# transitive deps (W2_arr) need the hypergeometric routines from KiLCA/slatec.
target_link_libraries(test_profile_input KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_profile_input
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input.x)

add_executable(test_profile_input_integration ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input_integration.f90)
set_target_properties(test_profile_input_integration PROPERTIES
                        OUTPUT_NAME test_profile_input_integration.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_profile_input_integration KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_profile_input_integration
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input_integration.x)

add_executable(test_kim_diagnostics ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_diagnostics.f90)
set_target_properties(test_kim_diagnostics PROPERTIES
                        OUTPUT_NAME test_kim_diagnostics.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_diagnostics KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_kim_diagnostics
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_diagnostics.x)

add_executable(test_flr2_fourier_kernel ${CMAKE_SOURCE_DIR}/KIM/tests/test_flr2_fourier_kernel.f90)
set_target_properties(test_flr2_fourier_kernel PROPERTIES
                        OUTPUT_NAME test_flr2_fourier_kernel.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_flr2_fourier_kernel KIM_lib kilca_lib lapack cerf ddeabm slatec)
configure_file(${CMAKE_SOURCE_DIR}/verification/oracles/rho_phi_kernels.dat
               ${CMAKE_BINARY_DIR}/tests/rho_phi_kernels.dat COPYONLY)
add_test(NAME test_flr2_fourier_kernel
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_flr2_fourier_kernel.x)
set_tests_properties(test_flr2_fourier_kernel PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

add_executable(test_ampere_matrices ${CMAKE_SOURCE_DIR}/KIM/tests/test_ampere_matrices.f90)
set_target_properties(test_ampere_matrices PROPERTIES
                        OUTPUT_NAME test_ampere_matrices.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ampere_matrices KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_ampere_matrices
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ampere_matrices.x)

# Pure-utility unit tests (Lagrange coefficients, binary search, index lookup,
# equidistant grid generation). Full math link set: pulling these KIM objects
# drags transitive deps (W2_arr) needing hypergeometric routines from KiLCA/slatec.
add_executable(test_plag_coeff ${CMAKE_SOURCE_DIR}/KIM/tests/test_plag_coeff.f90)
set_target_properties(test_plag_coeff PROPERTIES
                        OUTPUT_NAME test_plag_coeff.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_plag_coeff KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_plag_coeff
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_plag_coeff.x)

add_executable(test_findIndex ${CMAKE_SOURCE_DIR}/KIM/tests/test_findIndex.f90)
set_target_properties(test_findIndex PROPERTIES
                        OUTPUT_NAME test_findIndex.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_findIndex KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_findIndex
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_findIndex.x)

add_executable(test_grid_equidistant ${CMAKE_SOURCE_DIR}/KIM/tests/test_grid_equidistant.f90)
set_target_properties(test_grid_equidistant PROPERTIES
                        OUTPUT_NAME test_grid_equidistant.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_grid_equidistant KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_grid_equidistant
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_grid_equidistant.x)

# Regression test for the vendored periodization utility (make_periodic/localizer).
add_executable(test_periodization ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodization.f90)
set_target_properties(test_periodization PROPERTIES
                        OUTPUT_NAME test_periodization.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodization KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_periodization
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodization.x)

# Periodic sampling of KIM background primitives (aperfuns callback + wrapper).
add_executable(test_periodic_background ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_background.f90)
set_target_properties(test_periodic_background PROPERTIES
                        OUTPUT_NAME test_periodic_background.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_background KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_periodic_background
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_background.x)

# Build a periodic background plasma on a redirected window grid (Phase-2.3).
add_executable(test_periodic_background_build ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_background_build.f90)
set_target_properties(test_periodic_background_build PROPERTIES
                        OUTPUT_NAME test_periodic_background_build.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_background_build KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_periodic_background_build
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_background_build.x)
set_tests_properties(test_periodic_background_build PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Equilibrium parallel-flow contribution to harmonic-resolved FP x2, including
# the global and forced-periodic background paths.
add_executable(test_fp_parallel_flow ${CMAKE_SOURCE_DIR}/KIM/tests/test_fp_parallel_flow.f90)
set_target_properties(test_fp_parallel_flow PROPERTIES
                        OUTPUT_NAME test_fp_parallel_flow.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_fp_parallel_flow KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_fp_parallel_flow
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_fp_parallel_flow.x)
set_tests_properties(test_fp_parallel_flow PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Independent arbitrary-precision Mathematica oracle for every FP
# susceptibility consumed by KIM.
add_executable(test_fp_susceptibility_oracle
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_fp_susceptibility_oracle.f90)
set_target_properties(test_fp_susceptibility_oracle PROPERTIES
                        OUTPUT_NAME test_fp_susceptibility_oracle.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_fp_susceptibility_oracle
                      KIM_lib kilca_lib lapack cerf ddeabm slatec)
configure_file(${CMAKE_SOURCE_DIR}/verification/oracles/fp_susceptibilities.dat
               ${CMAKE_BINARY_DIR}/tests/fp_susceptibilities.dat COPYONLY)
add_test(NAME test_fp_susceptibility_oracle
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_fp_susceptibility_oracle.x)
set_tests_properties(test_fp_susceptibility_oracle PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Independent equilibrium/profile oracle, including multispecies collision
# frequencies, quasineutral array profiles, and cell-centre reconstruction.
add_executable(test_plasma_background_oracle
               ${CMAKE_SOURCE_DIR}/KIM/tests/test_plasma_background_oracle.f90)
set_target_properties(test_plasma_background_oracle PROPERTIES
                        OUTPUT_NAME test_plasma_background_oracle.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_plasma_background_oracle
                      KIM_lib kilca_lib lapack cerf ddeabm slatec)
configure_file(${CMAKE_SOURCE_DIR}/verification/oracles/plasma_background.dat
               ${CMAKE_BINARY_DIR}/tests/plasma_background.dat COPYONLY)
add_test(NAME test_plasma_background_oracle
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_plasma_background_oracle.x)
set_tests_properties(test_plasma_background_oracle PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Assemble the dense periodic Fourier matrices K_{m,m'} (Phase-2.4).
add_executable(test_periodic_assembly ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_assembly.f90)
set_target_properties(test_periodic_assembly PROPERTIES
                        OUTPUT_NAME test_periodic_assembly.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_assembly KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_periodic_assembly
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_assembly.x)
set_tests_properties(test_periodic_assembly PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Dense zgesv solve + inverse DFT for the periodic potential (Phase-2.5).
add_executable(test_periodic_solve ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_solve.f90)
set_target_properties(test_periodic_solve PROPERTIES
                        OUTPUT_NAME test_periodic_solve.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_solve KIM_lib kilca_lib lapack cerf ddeabm slatec)
configure_file(${CMAKE_SOURCE_DIR}/verification/oracles/field_operator_matrices.dat
               ${CMAKE_BINARY_DIR}/tests/field_operator_matrices.dat COPYONLY)
add_test(NAME test_periodic_solve
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_solve.x)
set_tests_properties(test_periodic_solve PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Lifecycle/contract tests for the kim_solver_m API seam.
add_executable(test_kim_solver ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver.f90)
set_target_properties(test_kim_solver PROPERTIES
                        OUTPUT_NAME test_kim_solver.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_kim_solver
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver.x)

# End-to-end keystone: a real electromagnetic solve through the kim_solver seam.
add_executable(test_kim_solver_em ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver_em.f90)
set_target_properties(test_kim_solver_em PROPERTIES
                        OUTPUT_NAME test_kim_solver_em.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver_em KIM_lib kilca_lib lapack cerf ddeabm slatec)
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
target_link_libraries(test_kim_solver_periodic KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_kim_solver_periodic
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver_periodic.x)
set_tests_properties(test_kim_solver_periodic PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

# Phase-2.7 cross-check: periodic vs global electrostatic delta-Phi on the same
# (m,n) = (-6,2) resonant case. Diagnostic (approx-vs-approx): logs the relative
# difference; soft gate (only on missing/all-zero/NaN).
add_executable(test_periodic_vs_global ${CMAKE_SOURCE_DIR}/KIM/tests/test_periodic_vs_global.f90)
set_target_properties(test_periodic_vs_global PROPERTIES
                        OUTPUT_NAME test_periodic_vs_global.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_periodic_vs_global KIM_lib kilca_lib lapack cerf ddeabm slatec)
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
target_link_libraries(test_periodic_convergence KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_periodic_convergence
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_periodic_convergence.x)
set_tests_properties(test_periodic_convergence PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)
