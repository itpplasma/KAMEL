# Common math routines library
# This builds all math-related Fortran sources shared by subprojects
file(GLOB_RECURSE COMMON_MATH_SOURCES
    "${CMAKE_SOURCE_DIR}/common/math/*.f90"
    "${CMAKE_SOURCE_DIR}/common/math/*.f"
)
add_library(slatec STATIC ${COMMON_MATH_SOURCES})
set_target_properties(slatec PROPERTIES LINKER_LANGUAGE Fortran)
# Expose module directory of these sources
target_include_directories(slatec PUBLIC "${CMAKE_SOURCE_DIR}/common/math")
