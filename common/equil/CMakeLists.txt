# Equilibrium field-line tracing library for KAMEL project
# Provides common equilibrium field-line tracing functionality for KIM, QL-Balance, and PreProc/fourier
#
# This library provides:
# - geqdsk equilibrium reading (via libneo)
# - Field-line tracing with RK4 integration
# - Magnetic axis finding
# - Safety factor q, flux coordinates, and profile computation

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Util)

# Disable libneo tests - they are not built (EXCLUDE_FROM_ALL) but would be registered with CTest
# Must also temporarily disable BUILD_TESTING because libneo overrides LIBNEO_BUILD_TESTING with it
set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "Disable libneo tests" FORCE)
set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "Disable libneo tests" FORCE)
set(LIBNEO_ENABLE_GOLDEN_TESTS OFF CACHE BOOL "Disable libneo golden tests" FORCE)
set(_KAMEL_BUILD_TESTING_SAVED ${BUILD_TESTING})
set(BUILD_TESTING OFF)

# Only fetch libneo if not already available (prevents duplicate targets)
if(NOT TARGET LIBNEO::magfie)
    find_or_fetch(libneo)
endif()

# Restore BUILD_TESTING for KAMEL's own tests
set(BUILD_TESTING ${_KAMEL_BUILD_TESTING_SAVED})

# Add netcdf library paths (needed because libneo's link_directories doesn't propagate)
find_program(NF_CONFIG "nf-config")
if(NF_CONFIG)
    execute_process(COMMAND nf-config --prefix
        OUTPUT_VARIABLE NETCDF_FORTRAN_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND nc-config --prefix
        OUTPUT_VARIABLE NETCDF_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    link_directories(${NETCDF_FORTRAN_PREFIX}/lib ${NETCDF_PREFIX}/lib)
endif()

set(EQUIL_MODULE_DIR ${CMAKE_BINARY_DIR}/modules/equil)
file(MAKE_DIRECTORY ${EQUIL_MODULE_DIR})

# Source files for the equilibrium library
set(EQUIL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/rk4_integrator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mag_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/field_line_rhs.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/equil_profiles.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/profile_preprocessor.f90
)

# Create the equilibrium library
add_library(kamel_equil STATIC ${EQUIL_SOURCES})

# Ensure libneo is built before kamel_equil (for module files)
add_dependencies(kamel_equil LIBNEO::magfie)

# Set properties for the library
set_target_properties(kamel_equil PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${EQUIL_MODULE_DIR}
    LINKER_LANGUAGE Fortran
)

# Include directories - libneo modules are in the build/common/equil/libneo directory
target_include_directories(kamel_equil
    PUBLIC
        ${EQUIL_MODULE_DIR}
        ${KAMEL_LIBNEO_BINARY_DIR}/include
)

# Link libneo magfie for equilibrium field functionality
target_link_libraries(kamel_equil
    PUBLIC
        LIBNEO::magfie
)

# Export for use by other KAMEL components
set(KAMEL_EQUIL_LIBRARIES kamel_equil CACHE INTERNAL "KAMEL equilibrium library")
set(KAMEL_EQUIL_INCLUDE_DIRS ${EQUIL_MODULE_DIR} CACHE INTERNAL "KAMEL equilibrium include directory")

# Tests
if(BUILD_TESTING)
    # Test for RK4 integrator
    # Note: This test links against kamel_equil which has libneo dependencies.
    # The test itself only uses the RK4 integrator which doesn't need libneo at runtime.
    add_executable(test_rk4_integrator
        ${CMAKE_CURRENT_SOURCE_DIR}/test/test_rk4_integrator.f90
    )
    target_link_libraries(test_rk4_integrator PRIVATE kamel_equil)
    target_include_directories(test_rk4_integrator PRIVATE ${EQUIL_MODULE_DIR})
    set_target_properties(test_rk4_integrator PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    )
    add_test(NAME test_rk4_integrator COMMAND test_rk4_integrator)

    # Test for profile preprocessor
    add_executable(test_profile_preprocessor
        ${CMAKE_CURRENT_SOURCE_DIR}/test/test_profile_preprocessor.f90
    )
    target_link_libraries(test_profile_preprocessor PRIVATE kamel_equil)
    target_include_directories(test_profile_preprocessor PRIVATE ${EQUIL_MODULE_DIR})
    set_target_properties(test_profile_preprocessor PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    )
    add_test(NAME test_profile_preprocessor COMMAND test_profile_preprocessor)
endif()
