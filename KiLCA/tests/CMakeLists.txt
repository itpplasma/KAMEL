configure_file(
    ${PROJECT_SOURCE_DIR}/tests/test_namelist.nml
    ${CMAKE_CURRENT_BINARY_DIR}/test_namelist.nml
    COPYONLY
)
add_executable(test_namelist
    test_namelist.cpp
)
target_include_directories(test_namelist PRIVATE "${PROJECT_SOURCE_DIR}/io")
target_link_libraries(test_namelist PRIVATE kilca_io)
add_test(NAME test_namelist
    COMMAND test_namelist
)
