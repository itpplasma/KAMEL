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
)

# Enable testing
enable_testing()

# Add test executables
add_executable(test_kilca_types tests/test_kilca_types.f90)
target_link_libraries(test_kilca_types kilca_fortran_modules)

add_executable(test_kilca_core tests/test_kilca_core.f90)
target_link_libraries(test_kilca_core kilca_fortran_modules)

# Add the tests
add_test(NAME test_kilca_types COMMAND test_kilca_types)
add_test(NAME test_kilca_core COMMAND test_kilca_core)

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