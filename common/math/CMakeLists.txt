set(SLATEC_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/slatec" CACHE PATH "Include directory for slatec library")
file(GLOB_RECURSE SLATEC_FILES "${SLATEC_INCLUDE_DIR}/*.f")
# fortnum_amos_compat supplies zbesj/zbesi/zbesk on the historical AMOS ABI,
# backed by the fortnum complex-Bessel core; it replaces the deleted bundled
# AMOS sources while the error/machine-constant routines (d1mach, xermsg, ...)
# stay for QUADPACK and ddeabm.
add_library(slatec STATIC ${SLATEC_FILES}
    "${CMAKE_CURRENT_SOURCE_DIR}/fortnum_amos_compat.f90")
set_target_properties(slatec PROPERTIES LINKER_LANGUAGE Fortran)
target_include_directories(slatec PUBLIC "${SLATEC_INCLUDE_DIR}")
target_link_libraries(slatec PUBLIC fortnum)
