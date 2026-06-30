# OpenMP runtime: KIM/libneo objects use OMP critical sections, so test
# executables that pull them need libgomp (GOMP_critical_*) at link time.
find_package(OpenMP REQUIRED COMPONENTS C Fortran)

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
target_link_libraries(test_region_roots_vs_muller KIM_lib lapack)
add_test(NAME test_region_roots_vs_muller
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_region_roots_vs_muller.x)

add_executable(test_equil_ode_dop853 ${CMAKE_SOURCE_DIR}/KIM/tests/test_equil_ode_dop853.f90)
set_target_properties(test_equil_ode_dop853 PROPERTIES
                        OUTPUT_NAME test_equil_ode_dop853.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_equil_ode_dop853 KIM_lib)
add_test(NAME test_equil_ode_dop853
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_equil_ode_dop853.x)

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
target_link_libraries(test_kim_diagnostics KIM_lib kilca_lib lapack cerf fortnum_amos_compat OpenMP::OpenMP_Fortran)
add_test(NAME test_kim_diagnostics
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_diagnostics.x)

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

# Lifecycle/contract tests for the kim_solver_m API seam.
add_executable(test_kim_solver ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver.f90)
set_target_properties(test_kim_solver PROPERTIES
                        OUTPUT_NAME test_kim_solver.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver KIM_lib kilca_lib lapack cerf fortnum_amos_compat OpenMP::OpenMP_Fortran)
add_test(NAME test_kim_solver
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver.x)

# End-to-end keystone: a real electromagnetic solve through the kim_solver seam.
add_executable(test_kim_solver_em ${CMAKE_SOURCE_DIR}/KIM/tests/test_kim_solver_em.f90)
set_target_properties(test_kim_solver_em PROPERTIES
                        OUTPUT_NAME test_kim_solver_em.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_kim_solver_em KIM_lib kilca_lib lapack cerf fortnum_amos_compat OpenMP::OpenMP_Fortran)
configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/KIM_config_em_small.nml
               ${CMAKE_BINARY_DIR}/tests/KIM_config_em_small.nml COPYONLY)
add_test(NAME test_kim_solver_em
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_kim_solver_em.x)
set_tests_properties(test_kim_solver_em PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)
