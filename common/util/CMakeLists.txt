# Utility library for KAMEL project
# Provides common utility functions for all KAMEL codes

# Create the utility library
add_library(kamel_util STATIC
    loading_bar.f90
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
set(KAMEL_UTIL_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/modules CACHE INTERNAL "KAMEL utility include directory")
