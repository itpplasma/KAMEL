cmake_minimum_required(VERSION 3.18)

# If this is a standalone build, set the compiler
if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "mpif90")
endif()

project(KiLCA_Fortran_Modules LANGUAGES Fortran CXX)

# Set build type if not set
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Inherit flags from parent or set defaults
if(NOT CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "-fopenmp -O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall -std=f2008 -Wno-compare-reals -fall-intrinsics")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -fcheck=all -std=f2008 -Wno-compare-reals -fall-intrinsics")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

# Set module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Set runtime output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Set archive output directory for libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Find required packages only if building standalone
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    find_package(GSL REQUIRED)
    set(FORTRAN_MODULE_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
else()
    # When built as part of KiLCA, use parent's dependencies
    set(FORTRAN_MODULE_LIBS lapack blas gsl gslcblas)
endif()

# Create the main library with all source files in dependency order
add_library(kilca_fortran_modules STATIC
    kilca_types_m.f90
    kilca_constants_m.f90
    kilca_shared_m.f90
    bessel_stubs.f90
    kilca_complex_m.f90
    kilca_settings_m.f90
    kilca_core_m.f90
    kilca_stability_m.f90
    kilca_spline_m.f90
    kilca_background_m.f90
    kilca_mode_m.f90
    kilca_eigensolver_m.f90
    kilca_orthogonalization_m.f90
    kilca_mode_solver_m.f90
    kilca_interp_integ_m.f90
    kilca_zone_m.f90
    kilca_eigtransform_m.f90
    kilca_rhs_func_m.f90
    kilca_memory_m.f90
    kilca_hypergeometric_basic_m.f90
    kilca_hypergeometric_m.f90
    kilca_fourier_m.f90
    kilca_antenna_m.f90
    kilca_bdf_solver_m.f90
    kilca_solver_m.f90
    kilca_conductivity_m.f90
    kilca_maxwell_m.f90
    kilca_physics_m.f90
    kilca_quants_m.f90
    kilca_zerofind_m.f90
    kilca_adaptive_grid_m.f90
    hyper1F1.cpp
)

# Set target properties
set_target_properties(kilca_fortran_modules PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    POSITION_INDEPENDENT_CODE ON
)

# Link libraries
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    target_link_libraries(kilca_fortran_modules PUBLIC 
        ${FORTRAN_MODULE_LIBS}
        GSL::gsl 
        GSL::gslcblas
    )
    target_include_directories(kilca_fortran_modules PUBLIC "${GSL_INCLUDE_DIRS}")
else()
    target_link_libraries(kilca_fortran_modules PUBLIC ${FORTRAN_MODULE_LIBS})
endif()

# Include directories
target_include_directories(kilca_fortran_modules PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include>
)

# Note: Bessel functions are handled by bessel_stubs.f90 in this module
# The ../math/bessel/bessel.f90 file has compatibility issues with modern Fortran standards

# Add main program executable
add_executable(kilca_main kilca_main.f90)
target_link_libraries(kilca_main kilca_fortran_modules)
target_include_directories(kilca_main PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Enable testing
enable_testing()

# Function to add test executable
function(add_kilca_test test_name test_source)
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} kilca_fortran_modules)
    target_include_directories(${test_name} PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# Add test executables
# General tests
file(GLOB GENERAL_TEST_SOURCES tests/test_*.f90)
foreach(test_source ${GENERAL_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_kilca_test(${test_name} ${test_source})
endforeach()

# Complex tests
file(GLOB COMPLEX_TEST_SOURCES tests/complex/test_*.f90)
foreach(test_source ${COMPLEX_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_kilca_test(${test_name}_complex ${test_source})
endforeach()

# Custom targets for test groups
add_custom_target(test-complex
    COMMAND ${CMAKE_CTEST_COMMAND} -R "_complex$"
    DEPENDS kilca_fortran_modules
    COMMENT "Running complex number tests..."
)

add_custom_target(test-general
    COMMAND ${CMAKE_CTEST_COMMAND} -E "_complex$"
    DEPENDS kilca_fortran_modules
    COMMENT "Running general tests..."
)

# Install targets
install(TARGETS kilca_fortran_modules
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/kilca
    FILES_MATCHING PATTERN "*.mod"
)