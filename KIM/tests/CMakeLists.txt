find_package(GSL REQUIRED)

# Optional: build QUADPACK/RKF45 integration comparison test
add_executable(test_integration_methods ${CMAKE_SOURCE_DIR}/KIM/tests/test_integration_methods.f90)
set_target_properties(test_integration_methods PROPERTIES
                        OUTPUT_NAME test_integration_methods.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_integration_methods KIM_lib lapack cerf slatec ddeabm)
add_test(NAME test_integration_methods
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_integration_methods.x)

# QUADPACK DQAGI sanity test (infinite integral)
add_executable(test_quadpack_qagi ${CMAKE_SOURCE_DIR}/KIM/tests/test_quadpack_qagi.f90)
set_target_properties(test_quadpack_qagi PROPERTIES
                        OUTPUT_NAME test_quadpack_qagi.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_quadpack_qagi KIM_lib lapack cerf slatec ddeabm)
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
target_link_libraries(test_profile_input KIM_lib)
add_test(NAME test_profile_input
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input.x)

add_executable(test_profile_input_integration ${CMAKE_SOURCE_DIR}/KIM/tests/test_profile_input_integration.f90)
set_target_properties(test_profile_input_integration PROPERTIES
                        OUTPUT_NAME test_profile_input_integration.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_profile_input_integration KIM_lib)
add_test(NAME test_profile_input_integration
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_profile_input_integration.x)

add_executable(test_ampere_matrices ${CMAKE_SOURCE_DIR}/KIM/tests/test_ampere_matrices.f90)
set_target_properties(test_ampere_matrices PROPERTIES
                        OUTPUT_NAME test_ampere_matrices.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ampere_matrices KIM_lib lapack cerf slatec ddeabm GSL::gsl GSL::gslcblas)
add_test(NAME test_ampere_matrices
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ampere_matrices.x)
