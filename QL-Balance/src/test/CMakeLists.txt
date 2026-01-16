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
