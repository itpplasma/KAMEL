foreach(test_name IN ITEMS test_inout test_zersol_bridge test_solver test_directory)
    add_executable(${test_name} tests/${test_name}.f90)
    target_link_libraries(${test_name} PRIVATE kilca_lib ${EXTERNAL_LIBS})
    target_include_directories(${test_name} PRIVATE ${PROJECT_BINARY_DIR}/OBJS/kilca/)
    set_target_properties(${test_name} PROPERTIES
        OUTPUT_NAME ${test_name}.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/OBJS/${test_name}/)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/${test_name})
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    set_tests_properties(${test_name} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/${test_name})
endforeach()
