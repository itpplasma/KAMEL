add_executable(test_zeal test_zeal.f90)
set_target_properties(test_zeal PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
)
target_link_libraries(test_zeal PUBLIC KIM_lib)
target_include_directories(test_zeal PUBLIC "${CMAKE_BINARY_DIR}/src")
add_test(NAME test_zeal COMMAND
    "${CMAKE_BINARY_DIR}/tests/test_zeal"
)

add_executable(test_kernel_rho_phi test_kernel_rho_phi.f90 kim_init_for_test.f90)
set_target_properties(test_kernel_rho_phi PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
)
target_link_libraries(test_kernel_rho_phi PUBLIC
    KIM_lib
    lapack
    cerf
    slatec
    ddeabm
)
target_include_directories(test_kernel_rho_phi PUBLIC "${CMAKE_BINARY_DIR}/src")
add_test(NAME test_kernel_rho_phi COMMAND
    "${CMAKE_BINARY_DIR}/tests/test_kernel_rho_phi.x" "${CMAKE_SOURCE_DIR}/KIM/nmls/KIM_config.nml"
)

add_executable(test_kernel_rho_B test_kernel_rho_B.f90 kim_init_for_test.f90)
set_target_properties(test_kernel_rho_B PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
)
target_link_libraries(test_kernel_rho_B PUBLIC
    KIM_lib
    lapack
    cerf
    slatec
    ddeabm
)
target_include_directories(test_kernel_rho_B PUBLIC "${CMAKE_BINARY_DIR}/src")
add_test(NAME test_kernel_rho_B COMMAND
    "${CMAKE_BINARY_DIR}/tests/test_kernel_rho_B" "${CMAKE_SOURCE_DIR}/KIM/nmls/KIM_config.nml"
)

# info: This is not a unit test, it needs additional files to run!
# add_executable(test_kernel_localizer test_kernel_localizer.f90 kim_init_for_test.f90)
# set_target_properties(test_kernel_localizer PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
# )
# target_link_libraries(test_kernel_localizer PUBLIC
#     KIM_lib
#     lapack
#     cerf
#     slatec
#     ddeabm
# )
# target_include_directories(test_kernel_localizer PUBLIC "${CMAKE_BINARY_DIR}/src")
# add_test(NAME test_kernel_localizer COMMAND
#     "${CMAKE_BINARY_DIR}/tests/test_kernel_localizer" "${CMAKE_SOURCE_DIR}/KIM/nmls/KIM_config.nml"
# )


add_executable(test_gauss_int test_gauss_int.f90)
set_target_properties(test_gauss_int PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
)
target_link_libraries(test_gauss_int PUBLIC
    KIM_lib
    lapack
    cerf
    slatec
    ddeabm
)
target_include_directories(test_gauss_int PUBLIC "${CMAKE_BINARY_DIR}/src")
add_test(NAME test_gauss_int COMMAND
    "${CMAKE_BINARY_DIR}/tests/test_gauss_int" "${CMAKE_SOURCE_DIR}/KIM/nmls/KIM_config.nml"
)
