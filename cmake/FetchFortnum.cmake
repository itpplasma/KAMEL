# fortnum is the single numerical core for KAMEL (special functions,
# quadrature, ODE, root finding). The Fortran target `fortnum` backs the
# KIM/KiLCA Fortran sources; ${fortnum_SOURCE_DIR}/include holds fortnum.h for
# the C/C++ sources. Added before libneo, which asks for fortnum as well.
itp_add_dependency(fortnum REF ${KAMEL_FORTNUM_REF}
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    TARGET fortnum)
FetchContent_GetProperties(fortnum)
