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

add_executable(test_profile_input ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input.f90)
set_target_properties(test_profile_input PROPERTIES
                        OUTPUT_NAME test_profile_input.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
# Links the full math set: using profile_input_m pulls KIM objects whose
# transitive deps (W2_arr) need the hypergeometric routines from KiLCA/slatec.
target_link_libraries(test_profile_input KIM_lib "${KiLCA_lib_path}" lapack cerf ddeabm slatec)
add_test(NAME test_profile_input
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input.x)

add_executable(test_profile_input_integration ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input_integration.f90)
set_target_properties(test_profile_input_integration PROPERTIES
                        OUTPUT_NAME test_profile_input_integration.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_profile_input_integration KIM_lib "${KiLCA_lib_path}" lapack cerf ddeabm slatec)
add_test(NAME test_profile_input_integration
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input_integration.x)

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
target_link_libraries(test_plag_coeff KIM_lib "${KiLCA_lib_path}" lapack cerf ddeabm slatec)
add_test(NAME test_plag_coeff
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_plag_coeff.x)

add_executable(test_findIndex ${CMAKE_SOURCE_DIR}/KIM/tests/test_findIndex.f90)
set_target_properties(test_findIndex PROPERTIES
                        OUTPUT_NAME test_findIndex.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_findIndex KIM_lib "${KiLCA_lib_path}" lapack cerf ddeabm slatec)
add_test(NAME test_findIndex
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_findIndex.x)

add_executable(test_grid_equidistant ${CMAKE_SOURCE_DIR}/KIM/tests/test_grid_equidistant.f90)
set_target_properties(test_grid_equidistant PROPERTIES
                        OUTPUT_NAME test_grid_equidistant.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_grid_equidistant KIM_lib "${KiLCA_lib_path}" lapack cerf ddeabm slatec)
add_test(NAME test_grid_equidistant
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_grid_equidistant.x)
