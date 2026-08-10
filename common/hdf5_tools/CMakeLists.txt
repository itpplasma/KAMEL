add_library(kamel_hdf5_tools STATIC kamel_hdf5_tools.f90)
target_link_libraries(kamel_hdf5_tools PUBLIC fortio::fortio)
set_target_properties(kamel_hdf5_tools PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)
target_include_directories(kamel_hdf5_tools PUBLIC ${CMAKE_BINARY_DIR}/modules)

set(KAMEL_HDF5_TOOLS_LIBRARIES kamel_hdf5_tools CACHE INTERNAL
    "KAMEL HDF5 tools library")
set(KAMEL_HDF5_TOOLS_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/modules CACHE INTERNAL
    "KAMEL HDF5 tools include directory")
