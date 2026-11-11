include(FetchContent)

# KAMEL uses only native pivoted real/complex LU. BLAS and optional ordering
# libraries are allowed; no external sparse factorizer or helper is provisioned.
if(NOT TARGET fortsparse)
    set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_ENABLE_SUPERLU OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_ENABLE_UMFPACK_IPC OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_ENABLE_MUMPS OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_ENABLE_NATIVE_CUDA OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(fortsparse
        GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
        GIT_TAG d5ab98a8b5d12550eb9aacf795eba97d782f1735
    )
    FetchContent_MakeAvailable(fortsparse)
endif()
