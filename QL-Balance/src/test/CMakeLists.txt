set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)


function (add_fortran_test name)
    add_test(NAME ${name} COMMAND ${name}.x)
    set_tests_properties(${name} PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
endfunction()

add_executable(test_sparse.x test_sparse.f90)
add_dependencies(test_sparse.x sparse)
target_link_libraries(test_sparse.x PUBLIC sparse SuiteSparse::umfpack_wrappers SuiteSparse::umfpack)
target_include_directories(test_sparse.x PUBLIC "${UMFPACK_INCLUDE_DIR}")
add_fortran_test(test_sparse)

# THE FOLLOWING TEST IS BROKEN, AS 'balance_conf.nml' IS NOT FOUND CORRECTLY
# add_executable(test_D_one_over_nu.x test_D_one_over_nu.f90)
# target_link_libraries(test_D_one_over_nu.x PUBLIC ql-balance_lib)
# add_fortran_test(test_D_one_over_nu)

# RHS Balance module tests
add_executable(test_rhs_balance.x test_rhs_balance.f90)
target_link_libraries(test_rhs_balance.x PUBLIC ql-balance_lib)
add_fortran_test(test_rhs_balance)

# Parameter-scan profile rescaling
add_executable(test_paramscan.x test_paramscan.f90)
target_link_libraries(test_paramscan.x PUBLIC ql-balance_lib)
add_fortran_test(test_paramscan)

# Integrated toroidal torque normalization
add_executable(test_toroidal_torque.x test_toroidal_torque.f90)
target_link_libraries(test_toroidal_torque.x PUBLIC ql-balance_lib)
add_fortran_test(test_toroidal_torque)

# KIM adapter module tests (unit tests for adapter helpers and linkage)
add_executable(test_kim_adapter.x test_kim_adapter.f90)
target_link_libraries(test_kim_adapter.x PUBLIC ql-balance_lib)
add_fortran_test(test_kim_adapter)
configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/KIM_config_em_small.nml
               ${CMAKE_BINARY_DIR}/tests/KIM_config_em_small.nml COPYONLY)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/profiles)
foreach(profile_file n.dat q.dat Te.dat Ti.dat)
    configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/${profile_file}
                   ${CMAKE_BINARY_DIR}/tests/profiles/${profile_file} COPYONLY)
endforeach()
set_tests_properties(test_kim_adapter PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

add_executable(test_periodic_embedding.x test_periodic_embedding.f90)
target_link_libraries(test_periodic_embedding.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_embedding)

add_executable(test_periodic_response_contract.x test_periodic_response_contract.f90)
target_link_libraries(test_periodic_response_contract.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_response_contract)

add_executable(test_periodic_kim_selector.x test_periodic_kim_selector.f90)
target_link_libraries(test_periodic_kim_selector.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_kim_selector)

add_executable(test_periodic_current_normalization.x test_periodic_current_normalization.f90)
target_link_libraries(test_periodic_current_normalization.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_current_normalization)

add_executable(test_periodic_amplitude_state.x test_periodic_amplitude_state.f90)
target_link_libraries(test_periodic_amplitude_state.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_amplitude_state)

add_executable(test_periodic_workflow_validation.x test_periodic_workflow_validation.f90)
target_link_libraries(test_periodic_workflow_validation.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_workflow_validation)
