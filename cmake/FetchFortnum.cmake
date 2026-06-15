include(FetchContent)

# fortnum is the single numerical core for KAMEL (special functions,
# quadrature, ODE, root finding). The Fortran target `fortnum` backs the
# KIM/KiLCA Fortran sources; ${fortnum_SOURCE_DIR}/include holds fortnum.h for
# the C/C++ sources.
#
# Guard against double declaration: libneo may pull fortnum transitively, so
# only declare the target if nobody else has provided it yet.
if(NOT TARGET fortnum)
    if(DEFINED ENV{CODE} AND EXISTS $ENV{CODE}/fortnum/CMakeLists.txt)
        message(STATUS "Using fortnum in $ENV{CODE}/fortnum")
        FetchContent_Declare(fortnum SOURCE_DIR $ENV{CODE}/fortnum)
    else()
        FetchContent_Declare(fortnum
            GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
            GIT_TAG 86c287877a7524586b542e3ec7f044b755e338b5
        )
    endif()
    FetchContent_MakeAvailable(fortnum)
endif()
