include(FetchContent)

FetchContent_Declare(
    sundials
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG v7.4.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)

set(SUNDIALS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SUNDIALS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_ARKODE OFF CACHE BOOL "" FORCE)
set(BUILD_CVODES OFF CACHE BOOL "" FORCE)
set(BUILD_IDA OFF CACHE BOOL "" FORCE)
set(BUILD_IDAS OFF CACHE BOOL "" FORCE)
set(BUILD_KINSOL OFF CACHE BOOL "" FORCE)
set(EXAMPLES_ENABLE_C OFF CACHE BOOL "" FORCE)
set(EXAMPLES_INSTALL OFF CACHE BOOL "" FORCE)
# Fortran 2003 module interface (fcvode_mod, fnvector_serial_mod, ...) so the
# KiLCA solver can drive CVODE from Fortran during the C++ -> Fortran port.
set(BUILD_FORTRAN_MODULE_INTERFACE ON CACHE BOOL "" FORCE)

# Show git's clone progress for this (large) fetch only.
set(_kamel_fc_quiet_old ${FETCHCONTENT_QUIET})
set(FETCHCONTENT_QUIET OFF)
FetchContent_MakeAvailable(sundials)
set(FETCHCONTENT_QUIET ${_kamel_fc_quiet_old})
