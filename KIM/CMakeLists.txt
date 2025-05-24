set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS}"
    "-O3"
    "-Wall"
    "-fbacktrace"
    "-fcheck=all"
    "-ffree-line-length-none"
    "-fopenmp"
    "-g"
    "-mtune=generic"
)
string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

enable_language(Fortran)

# Place Fortran modules for all KIM build files in one directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/KIM_mod")
# Make sure our module directory is on the include path
include_directories("${CMAKE_Fortran_MODULE_DIRECTORY}")

enable_testing()

add_subdirectory(src)
# slatec is provided by the top-level cmake/Dependencies.cmake; do not add its subdir here
add_subdirectory(src/math/libcerf-main)
add_subdirectory(src/math/ddeabm)
add_subdirectory(src/tests)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
