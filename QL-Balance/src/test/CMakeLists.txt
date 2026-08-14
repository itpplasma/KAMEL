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

add_executable(test_sparse_oracle.x test_sparse_oracle.f90)
target_link_libraries(test_sparse_oracle.x PRIVATE
    sparse SuiteSparse::umfpack_wrappers SuiteSparse::umfpack)
add_fortran_test(test_sparse_oracle)

# Singular system must stop with an UMFPACK error (ERROR STOP, nonzero exit).
# A normal return prints NO_ERROR_STOP and fails; a crash prints nothing and fails.
add_executable(test_sparse_singular.x test_sparse_singular.f90)
target_link_libraries(test_sparse_singular.x PRIVATE
    sparse SuiteSparse::umfpack_wrappers SuiteSparse::umfpack)
add_test(NAME test_sparse_singular COMMAND test_sparse_singular.x)
set_tests_properties(test_sparse_singular PROPERTIES
    PASS_REGULAR_EXPRESSION "UMFPACK singular matrix"
    FAIL_REGULAR_EXPRESSION "NO_ERROR_STOP")

# THE FOLLOWING TEST IS BROKEN, AS 'balance_conf.nml' IS NOT FOUND CORRECTLY
# add_executable(test_D_one_over_nu.x test_D_one_over_nu.f90)
# target_link_libraries(test_D_one_over_nu.x PUBLIC ql-balance_lib)
# add_fortran_test(test_D_one_over_nu)

# RHS Balance module tests
add_executable(test_rhs_balance.x test_rhs_balance.f90)
target_link_libraries(test_rhs_balance.x PUBLIC ql-balance_lib)
add_fortran_test(test_rhs_balance)

# A(y)*y + q(y) of the implicit step must count the RMP torque once
add_executable(test_rhs_source_consistency.x test_rhs_source_consistency.f90)
target_link_libraries(test_rhs_source_consistency.x PUBLIC ql-balance_lib)
add_fortran_test(test_rhs_source_consistency)

# One implicit step must keep a stationary state (accuracy of the sparse solve)
add_executable(test_evolvestep_stationary.x test_evolvestep_stationary.f90)
target_link_libraries(test_evolvestep_stationary.x PUBLIC ql-balance_lib)
add_fortran_test(test_evolvestep_stationary)

# Parameter-scan profile rescaling
add_executable(test_paramscan.x test_paramscan.f90)
target_link_libraries(test_paramscan.x PUBLIC ql-balance_lib)
add_fortran_test(test_paramscan)

# Da interpolation from profiles/Da.dat onto the balance grid
add_executable(test_diff_coeffs.x test_diff_coeffs.f90)
target_link_libraries(test_diff_coeffs.x PUBLIC ql-balance_lib)
add_fortran_test(test_diff_coeffs)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diff_coeffs/profiles)
set_tests_properties(test_diff_coeffs PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diff_coeffs)

# Namelist text stored as /run_info/balancenml
add_executable(test_run_info.x test_run_info.f90)
target_link_libraries(test_run_info.x PUBLIC ql-balance_lib)
add_fortran_test(test_run_info)

# Integrated toroidal torque normalization
add_executable(test_toroidal_torque.x test_toroidal_torque.f90)
target_link_libraries(test_toroidal_torque.x PUBLIC ql-balance_lib)
add_fortran_test(test_toroidal_torque)

# Accepted-state radial electric field and NTV profile coupling
add_executable(test_ntv_feedback.x test_ntv_feedback.f90)
target_link_libraries(test_ntv_feedback.x PUBLIC ql-balance_lib)
add_fortran_test(test_ntv_feedback)

# KIM adapter module tests (unit tests for adapter helpers and linkage)
add_executable(test_kim_adapter.x test_kim_adapter.f90)
target_link_libraries(test_kim_adapter.x PUBLIC ql-balance_lib)
add_fortran_test(test_kim_adapter)

# Velocity-space integral port (vel_integral_m) regression test
add_executable(test_vel_integral.x test_vel_integral.f90)
target_link_libraries(test_vel_integral.x PUBLIC ql-balance_lib)
add_fortran_test(test_vel_integral)
configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/KIM_config_em_small.nml
               ${CMAKE_BINARY_DIR}/tests/KIM_config_em_small.nml COPYONLY)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/profiles)
foreach(profile_file n.dat q.dat Te.dat Ti.dat)
    configure_file(${CMAKE_SOURCE_DIR}/KIM/tests/test_data/${profile_file}
                   ${CMAKE_BINARY_DIR}/tests/profiles/${profile_file} COPYONLY)
endforeach()
set_tests_properties(test_kim_adapter PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

add_executable(test_hdf5_diagnostic_writers.x test_hdf5_diagnostic_writers.f90)
target_link_libraries(test_hdf5_diagnostic_writers.x PUBLIC ql-balance_lib)
add_fortran_test(test_hdf5_diagnostic_writers)
set_tests_properties(test_hdf5_diagnostic_writers PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

add_executable(test_periodic_embedding.x test_periodic_embedding.f90)
target_link_libraries(test_periodic_embedding.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_embedding)

add_executable(test_periodic_response_contract.x test_periodic_response_contract.f90)
target_link_libraries(test_periodic_response_contract.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_response_contract)

add_executable(test_periodic_kim_selector.x test_periodic_kim_selector.f90)
target_link_libraries(test_periodic_kim_selector.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_kim_selector)

add_executable(test_periodic_kim_coupling.x test_periodic_kim_coupling.f90)
target_link_libraries(test_periodic_kim_coupling.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_kim_coupling)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_kim_coupling)
set_tests_properties(test_periodic_kim_coupling PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_kim_coupling)

add_executable(test_transport_smoothing.x test_transport_smoothing.f90)
target_link_libraries(test_transport_smoothing.x PUBLIC ql-balance_lib)
add_fortran_test(test_transport_smoothing)
add_executable(test_periodic_current_normalization.x test_periodic_current_normalization.f90)
target_link_libraries(test_periodic_current_normalization.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_current_normalization)

add_executable(test_periodic_target_current.x test_periodic_target_current.f90)
target_link_libraries(test_periodic_target_current.x PUBLIC ql-balance_lib)
target_include_directories(test_periodic_target_current.x PRIVATE ${CMAKE_BINARY_DIR}/OBJS/kim)
add_fortran_test(test_periodic_target_current)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_target_current)
set_tests_properties(test_periodic_target_current PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_target_current)

add_executable(test_periodic_response_normalization.x test_periodic_response_normalization.f90)
target_link_libraries(test_periodic_response_normalization.x PUBLIC ql-balance_lib)
target_include_directories(test_periodic_response_normalization.x PRIVATE
    ${CMAKE_BINARY_DIR}/OBJS/kim)
add_fortran_test(test_periodic_response_normalization)

foreach(invalid_case geometry_duplicate geometry_descending geometry_nan geometry_inf
        bounds_nan bounds_inf bounds_reversed bounds_empty bounds_outside grid_shape)
    add_test(NAME test_periodic_current_${invalid_case}
        COMMAND test_periodic_current_normalization.x ${invalid_case})
    set_tests_properties(test_periodic_current_${invalid_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

add_executable(test_periodic_multimode_feedback.x test_periodic_multimode_feedback.f90)
target_link_libraries(test_periodic_multimode_feedback.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_multimode_feedback)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_multimode_feedback)
set_tests_properties(test_periodic_multimode_feedback PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_multimode_feedback)

foreach(case wrong_sign outside)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_reject_${case})
    add_test(NAME test_periodic_reject_${case}
        COMMAND ${CMAKE_COMMAND}
            -DTEST_EXECUTABLE=$<TARGET_FILE:test_periodic_multimode_feedback.x>
            -DCASE=${case}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/check_periodic_nonresonant.cmake)
    set_tests_properties(test_periodic_reject_${case} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_reject_${case})
endforeach()

add_executable(test_periodic_multimode_assembly.x test_periodic_multimode_assembly.f90)
target_link_libraries(test_periodic_multimode_assembly.x PUBLIC ql-balance_lib)
target_include_directories(test_periodic_multimode_assembly.x PRIVATE ${CMAKE_BINARY_DIR}/OBJS/kim)
add_fortran_test(test_periodic_multimode_assembly)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_multimode_assembly)
set_tests_properties(test_periodic_multimode_assembly PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_multimode_assembly)
add_executable(test_periodic_amplitude_state.x test_periodic_amplitude_state.f90)
target_link_libraries(test_periodic_amplitude_state.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_amplitude_state)

add_executable(test_periodic_restart_adapter.x test_periodic_restart_adapter.f90)
target_link_libraries(test_periodic_restart_adapter.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_restart_adapter)
set_tests_properties(test_periodic_restart_adapter PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/)

foreach(case mode_count metadata_shape guard_no_stat invalid_initialize reject_fresh)
    add_test(NAME test_periodic_amplitude_${case}
        COMMAND test_periodic_amplitude_state.x ${case})
    set_tests_properties(test_periodic_amplitude_${case} PROPERTIES WILL_FAIL TRUE)
endforeach()
foreach(test periodic_checkpoint periodic_checkpoint_reader periodic_time_evolution time_history_restart
        timestep_controller)
    add_executable(test_${test}.x test_${test}.f90)
    target_link_libraries(test_${test}.x PUBLIC ql-balance_lib)
    add_fortran_test(test_${test})
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/${test})
    set_tests_properties(test_${test} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/${test})
endforeach()

foreach(case version phase partial modes status duplicate nonfinite)
    add_test(NAME test_periodic_checkpoint_${case}
        COMMAND test_periodic_checkpoint.x ${case})
    set_tests_properties(test_periodic_checkpoint_${case} PROPERTIES WILL_FAIL TRUE
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/periodic_checkpoint)
endforeach()

add_executable(test_periodic_workflow_validation.x test_periodic_workflow_validation.f90)
target_link_libraries(test_periodic_workflow_validation.x PUBLIC ql-balance_lib)
add_fortran_test(test_periodic_workflow_validation)
foreach(case duplicate shape bparallel)
    add_test(NAME test_periodic_workflow_${case}
        COMMAND test_periodic_workflow_validation.x ${case})
    set_tests_properties(test_periodic_workflow_${case} PROPERTIES WILL_FAIL TRUE)
endforeach()

# Ipar history of TimeEvolution is initialised and records the selected Ipar
add_executable(test_ipar_history.x test_ipar_history.f90)
target_link_libraries(test_ipar_history.x PUBLIC ql-balance_lib)
add_fortran_test(test_ipar_history)
