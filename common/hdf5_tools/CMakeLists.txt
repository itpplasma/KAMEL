add_library(kamel_hdf5_tools STATIC kamel_hdf5_tools.f90)
# The Fortran module is needed to build/link this library, but Fortio's
# zlib dependency must not export Apple's SDK C include directory to unrelated
# downstream C++ targets (libc++'s <cmath>/<cstdio> then fails on macOS).
target_link_libraries(kamel_hdf5_tools PUBLIC $<LINK_ONLY:fortio::fortio>)
target_include_directories(kamel_hdf5_tools PRIVATE
    "${fortio_BINARY_DIR}/modules")
set_target_properties(kamel_hdf5_tools PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)
target_include_directories(kamel_hdf5_tools PUBLIC ${CMAKE_BINARY_DIR}/modules)

set(KAMEL_HDF5_TOOLS_LIBRARIES kamel_hdf5_tools CACHE INTERNAL
    "KAMEL HDF5 tools library")
set(KAMEL_HDF5_TOOLS_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/modules CACHE INTERNAL
    "KAMEL HDF5 tools include directory")
