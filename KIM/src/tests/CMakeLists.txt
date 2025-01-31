add_executable(test_zeal test_zeal.f90)

set(KIM_lib "${CMAKE_BINARY_DIR}/lib/KIM_lib.a")
add_dependencies(test_zeal KIM_lib)

set_target_properties(test_zeal PROPERTIES 
                        OUTPUT_NAME test_zeal.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_zeal PUBLIC KIM_lib)
target_include_directories(test_zeal PUBLIC "${CMAKE_BINARY_DIR}/src")

add_test(NAME test_zeal COMMAND "${CMAKE_BINARY_DIR}/tests/test_zeal.x")