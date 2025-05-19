# simplified version
configure_file(
    ${PROJECT_SOURCE_DIR}/tests/data/simplified_namelist.nml
    ${CMAKE_CURRENT_BINARY_DIR}/simplified_namelist.nml
    COPYONLY
)
add_executable(test_namelist_simplified
    test_namelist_simplified.cpp
)
target_include_directories(test_namelist_simplified PRIVATE "${PROJECT_SOURCE_DIR}/io")
target_link_libraries(test_namelist_simplified PRIVATE kilca_io)
add_test(NAME test_namelist_simplified COMMAND test_namelist_simplified)

# full version used for real runs
configure_file(
    ${PROJECT_SOURCE_DIR}/tests/data/kilca_config.nml
    ${CMAKE_CURRENT_BINARY_DIR}/kilca_config.nml
    COPYONLY
)
add_executable(test_namelist_full
    test_namelist_full.cpp
)
target_include_directories(test_namelist_full PRIVATE "${PROJECT_SOURCE_DIR}/io")
target_link_libraries(test_namelist_full PRIVATE kilca_io)
add_test(NAME test_namelist_full COMMAND test_namelist_full)
