# HDF5 tools library for KAMEL project
# Provides common HDF5 wrapper functions for all KAMEL codes

# Find HDF5 with Fortran support
cmake_policy(SET CMP0074 NEW)
find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)

if (NOT ${HDF5_FOUND})
    message(FATAL_ERROR "HDF5 with Fortran support not found. Please set environment variable HDF5_ROOT to the correct path")
endif()

# Create the HDF5 tools library
add_library(kamel_hdf5_tools STATIC
    hdf5_tools.f90
    hdf5_tools_f2003.f90
)

# Set properties for the library
set_target_properties(kamel_hdf5_tools PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# Include directories
target_include_directories(kamel_hdf5_tools
    PUBLIC
        ${HDF5_INCLUDE_DIRS}
        ${CMAKE_BINARY_DIR}/modules
)

# Link HDF5 libraries
target_link_libraries(kamel_hdf5_tools
    PUBLIC
        ${HDF5_Fortran_LIBRARIES}
        ${HDF5_Fortran_HL_LIBRARIES}
)

# Export for use by other KAMEL components
set(KAMEL_HDF5_TOOLS_LIBRARIES kamel_hdf5_tools CACHE INTERNAL "KAMEL HDF5 tools library")
set(KAMEL_HDF5_TOOLS_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/modules CACHE INTERNAL "KAMEL HDF5 tools include directory")
