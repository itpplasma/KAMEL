# fortnum_amos_compat supplies zbesj/zbesi/zbesk on the historical AMOS ABI,
# backed by the fortnum complex-Bessel core. The bundled SLATEC error and
# machine-constant routines that QUADPACK and ddeabm needed are gone, so this
# is the only remaining math source here.
add_library(fortnum_amos_compat STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/fortnum_amos_compat.f90")
set_target_properties(fortnum_amos_compat PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(fortnum_amos_compat PUBLIC fortnum)
