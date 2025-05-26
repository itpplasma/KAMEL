# Dependencies for the unified KAMEL project
# Fetch or find external libraries once, for all subprojects

# Set up module path for fetch modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Core dependencies
include(FetchGSL)
include(FetchLapack)
include(FetchNetcdf)
include(FetchSuiteSparse)
include(FetchSundials)
include(FetchZeal)

# QL-Balance sparse module (shared dependency)
set(QLBALANCE_BASE "${CMAKE_SOURCE_DIR}/QL-Balance/src/base")
file(GLOB SPARSE_SOURCES
    "${QLBALANCE_BASE}/sparse_mod.f90"
)
add_library(sparse STATIC ${SPARSE_SOURCES})
set_target_properties(sparse PROPERTIES LINKER_LANGUAGE Fortran
                        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/install/lib/"
                        Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/OBJS/sparse/")
                        
# Link sparse to SuiteSparse UMFPACK
target_link_libraries(sparse PUBLIC SuiteSparse::umfpack_wrappers SuiteSparse::umfpack)
# Make module directory available
target_include_directories(sparse PUBLIC "${QLBALANCE_BASE}")

# Utility macros
include(Util)
