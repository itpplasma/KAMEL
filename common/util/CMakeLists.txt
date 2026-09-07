# Utility library for KAMEL project
# Provides common utility functions for all KAMEL codes

# Create the utility library
add_library(kamel_util STATIC
    loading_bar.f90
    terminal_width.c
)

# Set properties for the library
set_target_properties(kamel_util PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# Include directories
target_include_directories(kamel_util
    PUBLIC
        ${CMAKE_BINARY_DIR}/modules
)

# Export for use by other KAMEL components
set(KAMEL_UTIL_LIBRARIES kamel_util CACHE INTERNAL "KAMEL utility library")
set(KAMEL_UTIL_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/modules CACHE INTERNAL
    "KAMEL utility include directory")

if(BUILD_TESTING)
    add_executable(test_terminal_width
        test/test_terminal_width.c
        terminal_width.c
    )
    add_test(NAME test_terminal_width COMMAND test_terminal_width)

    add_executable(test_loading_bar test/test_loading_bar.f90)
    target_link_libraries(test_loading_bar PRIVATE kamel_util)
    add_test(NAME test_loading_bar COMMAND test_loading_bar)
    set_tests_properties(test_loading_bar PROPERTIES
        PASS_REGULAR_EXPRESSION "#####      BAR_TEST_PASS"
    )
endif()
