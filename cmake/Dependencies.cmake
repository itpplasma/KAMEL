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
        GIT_TAG fc16c50b53cfc8b92b5ba21f26b1b3566915cfd7
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
            get_target_property(_zlib_system_includes ${_zlib_real} INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
            if(_zlib_system_includes)
                list(FILTER _zlib_system_includes EXCLUDE REGEX "/MacOSX\\.sdk/usr/include/?$")
                set_property(TARGET ${_zlib_real} PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_zlib_system_includes}")
            endif()
        endif()
        if(TARGET fortio)
            get_target_property(_fortio_links fortio INTERFACE_LINK_LIBRARIES)
            if(_fortio_links)
                list(TRANSFORM _fortio_links REPLACE "^ZLIB::ZLIB$" "$<LINK_ONLY:ZLIB::ZLIB>")
                set_property(TARGET fortio PROPERTY INTERFACE_LINK_LIBRARIES "${_fortio_links}")
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
