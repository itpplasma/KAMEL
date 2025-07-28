cmake_minimum_required(VERSION 3.18)

# If this is a standalone build, set the compiler
if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "mpif90")
endif()

project(KiLCA_Fortran_Modules LANGUAGES Fortran)

# Inherit flags from parent or set defaults
if(NOT CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "-fopenmp -O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -fcheck=all")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

# Create the main library
add_library(kilca_fortran_modules STATIC
    kilca_types_m.f90
    kilca_core_m.f90
    kilca_memory_m.f90
    kilca_settings_m.f90
)

# Enable testing
enable_testing()

# Add test executables
add_executable(test_kilca_types tests/test_kilca_types.f90)
target_link_libraries(test_kilca_types kilca_fortran_modules)

add_executable(test_kilca_core tests/test_kilca_core.f90)
target_link_libraries(test_kilca_core kilca_fortran_modules)

add_executable(test_memory_management tests/test_memory_management.f90)
target_link_libraries(test_memory_management kilca_fortran_modules)

add_executable(test_kilca_settings tests/test_kilca_settings.f90)
target_link_libraries(test_kilca_settings kilca_fortran_modules)

add_executable(test_kilca_settings_complete tests/test_kilca_settings_complete.f90)
target_link_libraries(test_kilca_settings_complete kilca_fortran_modules)

add_executable(test_core_initialization tests/test_core_initialization.f90)
target_link_libraries(test_core_initialization kilca_fortran_modules)

add_executable(test_core_deep_copy tests/test_core_deep_copy.f90)
target_link_libraries(test_core_deep_copy kilca_fortran_modules)

add_executable(test_core_validation tests/test_core_validation.f90)
target_link_libraries(test_core_validation kilca_fortran_modules)

add_executable(test_settings_file_parsing tests/test_settings_file_parsing.f90)
target_link_libraries(test_settings_file_parsing kilca_fortran_modules)

add_executable(test_settings_validation tests/test_settings_validation.f90)
target_link_libraries(test_settings_validation kilca_fortran_modules)

add_executable(test_settings_output tests/test_settings_output.f90)
target_link_libraries(test_settings_output kilca_fortran_modules)

add_executable(test_settings_copy_compare tests/test_settings_copy_compare.f90)
target_link_libraries(test_settings_copy_compare kilca_fortran_modules)

add_executable(test_settings_initialization tests/test_settings_initialization.f90)
target_link_libraries(test_settings_initialization kilca_fortran_modules)

add_executable(test_settings_error_handling tests/test_settings_error_handling.f90)
target_link_libraries(test_settings_error_handling kilca_fortran_modules)

add_executable(test_settings_documentation tests/test_settings_documentation.f90)
target_link_libraries(test_settings_documentation kilca_fortran_modules)

# Add the tests
add_test(NAME test_kilca_types COMMAND test_kilca_types)
add_test(NAME test_kilca_core COMMAND test_kilca_core)
add_test(NAME test_memory_management COMMAND test_memory_management)
add_test(NAME test_kilca_settings COMMAND test_kilca_settings)
add_test(NAME test_kilca_settings_complete COMMAND test_kilca_settings_complete)
add_test(NAME test_core_initialization COMMAND test_core_initialization)
add_test(NAME test_core_deep_copy COMMAND test_core_deep_copy)
add_test(NAME test_core_validation COMMAND test_core_validation)
add_test(NAME test_settings_file_parsing COMMAND test_settings_file_parsing)
add_test(NAME test_settings_validation COMMAND test_settings_validation)
add_test(NAME test_settings_output COMMAND test_settings_output)
add_test(NAME test_settings_copy_compare COMMAND test_settings_copy_compare)
add_test(NAME test_settings_initialization COMMAND test_settings_initialization)
add_test(NAME test_settings_error_handling COMMAND test_settings_error_handling)
add_test(NAME test_settings_documentation COMMAND test_settings_documentation)

# Set module directory
set_target_properties(kilca_fortran_modules PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

# Include module directory for tests
target_include_directories(test_kilca_types PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_kilca_core PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_memory_management PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_kilca_settings PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_kilca_settings_complete PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_core_initialization PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_core_deep_copy PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_core_validation PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_settings_file_parsing PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_settings_validation PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_settings_output PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_settings_copy_compare PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_settings_initialization PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(test_settings_error_handling PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/modules
)