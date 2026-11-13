# Dependencies for the unified KAMEL project
# Fetch or find external libraries once, for all subprojects

# Set up module path for fetch modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Core dependencies
include(FetchLapack)
include(FetchSuiteSparse)
include(FetchSUNDIALS)
include(FetchFortnum)

include(FetchContent)
if(NOT TARGET fortio)
    if(APPLE)
        execute_process(COMMAND xcrun --show-sdk-path OUTPUT_VARIABLE _kamel_sdk
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(ZLIB_INCLUDE_DIR "${_kamel_sdk}/usr/include" CACHE PATH "" FORCE)
    endif()
    set(_KAMEL_BUILD_TESTING_SAVED ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG 11afd0bd1af0c99ea4e9a1c0df683dbdcc299b69
    )
    FetchContent_MakeAvailable(fortio)
    if(APPLE AND TARGET ZLIB::ZLIB)
        get_target_property(_zlib_real ZLIB::ZLIB ALIASED_TARGET)
        if(NOT _zlib_real)
            set(_zlib_real ZLIB::ZLIB)
        endif()
        if(_zlib_real)
            get_target_property(_zlib_includes ${_zlib_real} INTERFACE_INCLUDE_DIRECTORIES)
            if(_zlib_includes)
                list(FILTER _zlib_includes EXCLUDE REGEX "/MacOSX\\.sdk/usr/include/?$")
                set_property(TARGET ${_zlib_real} PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${_zlib_includes}")
            endif()
        endif()
    endif()
    set(BUILD_TESTING ${_KAMEL_BUILD_TESTING_SAVED})
endif()

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

# NEO-RT's standalone mode adds libneo unconditionally.  Configure it before
# KAMEL's other libneo consumers so the whole build shares one libneo target
# namespace instead of trying to add a second copy later.
#
# That shared libneo must be KAMEL's validated pin, not NEO-RT's own: NEO-RT
# resolves libneo through ItpDeps, which honours ITP_LIBNEO_REF. common/equil
# reuses the same LIBNEO_REF when it fetches libneo itself.
set(LIBNEO_REF "d03789a1fef15e188752614328dd50ea9854919f" CACHE STRING "libneo revision")
set(ITP_LIBNEO_REF "${LIBNEO_REF}" CACHE STRING "libneo revision for ItpDeps consumers")
# Pin NEO-RT like libneo, so an unchanged KAMEL commit always builds the same
# NEO-RT (find_or_fetch reads <NAME>_REF).
set(NEO-RT_REF "da2b037c4e662339f8caa794e5aaf48c45c4cbd6" CACHE STRING "NEO-RT revision")
set(NEORT_BUILD_TESTING OFF CACHE BOOL "Build tests for NEO-RT" FORCE)
set(_KAMEL_BUILD_TESTING_SAVED ${BUILD_TESTING})
set(BUILD_TESTING OFF)
if(DEFINED LIBNEO_PATH AND EXISTS "${LIBNEO_PATH}/CMakeLists.txt")
    get_filename_component(_KAMEL_LIBNEO_PATH "${LIBNEO_PATH}" ABSOLUTE
        BASE_DIR "${CMAKE_SOURCE_DIR}")
    set(FETCHCONTENT_SOURCE_DIR_LIBNEO "${_KAMEL_LIBNEO_PATH}")
endif()
find_or_fetch(neo-rt)
set(BUILD_TESTING ${_KAMEL_BUILD_TESTING_SAVED})
