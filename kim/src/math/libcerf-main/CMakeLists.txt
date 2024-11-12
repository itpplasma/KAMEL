cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
if(NOT DEFINED LIB_MAN)
    option(LIB_MAN "Build libcerf man pages" ON)
endif()
if(NOT DEFINED LIB_RUN)
    option(LIB_RUN "Build executables for command-line computation" ON)
endif()
if(NOT DEFINED CERF_CPP)
    option(CERF_CPP "Build libcerf as native C++" OFF)
endif()

if(CERF_CPP)
    project(cerfcpp CXX)
    message("Build C++ library libcerfcpp")
    # add_compile_definitions(CERF_AS_CPP) # requires CMake 3.13, therefore more explicitly:
    if(MSVC)
        add_compile_options(/DCERF_AS_CPP)
    else()
        add_compile_options(-DCERF_AS_CPP)
    endif()
    set(CMAKE_CXX_STANDARD 11)
else()
    message("Build C library libcerf")
    project(cerf C)
    set(CMAKE_C_STANDARD 99)
endif()

set(CERF_SOVERSION                 2) # API version
set(CERF_VERSION ${CERF_SOVERSION}.4) # minor version

include(GNUInstallDirs)

SET(CERF_COMPILE_OPTIONS "" CACHE
    STRING "User-supplied compiler options, will be appended to hard-coded ones. Not for Windows.")

if(MSVC)
    if (NOT CERF_CPP)
        message(FATAL_ERROR "Under MSVC, only CERF_CPP=ON is supported")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel compilation
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS}")
    set(MS_NOWARN "/wd4018 /wd4068 /wd4101 /wd4244 /wd4267 /wd4305 /wd4715 /wd4996")
    set(MS_D "-D_CRT_SECURE_NO_WARNINGS -D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${MS_NOWARN} ${MS_D}")
    set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
    if(BUILD_SHARED_LIBS)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    else()
        # required for post-build validation under vcpkg:
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd") # multithreaded, debug
    endif()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

elseif(CERF_COMPILE_OPTIONS)
    add_compile_options(${CERF_COMPILE_OPTIONS})

else()
    option(PEDANTIC "Compile with pedantic warnings" ON)
    option(WERROR "Treat warnings as errors" OFF)
    if(PEDANTIC)
        add_compile_options(-pedantic -Wall)
    endif()
    if(WERROR)
        add_compile_options(-Werror)
    endif()

    add_compile_options(-Wno-sign-compare -fno-omit-frame-pointer)
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-g)
    else()
        add_compile_options(-O3)
    endif()
endif()

include(CTest)

configure_file("libcerf.pc.in" "libcerf.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libcerf.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

add_subdirectory(lib)
add_subdirectory(test)
if(LIB_RUN)
    add_subdirectory(run)
endif()
if(LIB_MAN)
    add_subdirectory(man)
endif()
