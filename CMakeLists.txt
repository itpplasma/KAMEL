cmake_minimum_required(VERSION 3.24)

# Set compilers BEFORE project() - only if not already set by user
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    if(NOT DEFINED CMAKE_C_COMPILER)
        # 1. Homebrew LLVM first (Apple Clang lacks OpenMP support)
        find_program(_CLANG_C clang
            PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin
            NO_DEFAULT_PATH)
        if(NOT _CLANG_C)
            # 2. Fall back to PATH (Nix, manual installs, etc.)
            find_program(_CLANG_C_PATH clang)
            set(_CLANG_C "${_CLANG_C_PATH}")
        endif()
        if(_CLANG_C)
            message(STATUS "C compiler (Darwin): ${_CLANG_C}")
            set(CMAKE_C_COMPILER "${_CLANG_C}")
        endif()
    endif()
    if(NOT DEFINED CMAKE_CXX_COMPILER)
        # 1. Homebrew LLVM first (Apple Clang lacks OpenMP support)
        find_program(_CLANG_CXX clang++
            PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin
            NO_DEFAULT_PATH)
        if(NOT _CLANG_CXX)
            # 2. Fall back to PATH (Nix, manual installs, etc.)
            find_program(_CLANG_CXX_PATH clang++)
            set(_CLANG_CXX "${_CLANG_CXX_PATH}")
        endif()
        if(_CLANG_CXX)
            message(STATUS "CXX compiler (Darwin): ${_CLANG_CXX}")
            set(CMAKE_CXX_COMPILER "${_CLANG_CXX}")
        endif()
    endif()
    # Append Homebrew prefix paths only if brew is available (don't overwrite)
    find_program(_BREW_EXE brew HINTS /opt/homebrew/bin /usr/local/bin)
    if(_BREW_EXE)
        list(APPEND CMAKE_PREFIX_PATH /opt/homebrew /usr/local)
    endif()
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    if(NOT DEFINED CMAKE_C_COMPILER)
        set(CMAKE_C_COMPILER gcc)
    endif()
    if(NOT DEFINED CMAKE_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER g++)
    endif()
endif()

if(APPLE AND NOT DEFINED HDF5_DIR)
    # 1. Check if HDF5 is findable via standard paths / CMAKE_PREFIX_PATH
    find_package(HDF5 COMPONENTS C Fortran HL QUIET CONFIG)
    if(HDF5_FOUND)
        message(STATUS "HDF5 found via standard paths: ${HDF5_VERSION}")
    else()
        # 2. Homebrew fallback
        find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
        if(BREW_EXECUTABLE)
            execute_process(
                COMMAND "${BREW_EXECUTABLE}" --prefix hdf5
                OUTPUT_VARIABLE _hdf5_prefix
                OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE _hdf5_res
            )
            if(_hdf5_res EQUAL 0 AND EXISTS "${_hdf5_prefix}/lib/cmake/hdf5")
                set(HDF5_DIR "${_hdf5_prefix}/lib/cmake/hdf5" CACHE PATH "Path to HDF5Config.cmake")
                message(STATUS "HDF5_DIR set from Homebrew: ${HDF5_DIR}")
            endif()
        endif()
        # 3. If still not found, FetchHDF5.cmake will build from source
    endif()
endif()

project(KAMEL LANGUAGES C CXX Fortran)

# Set C++ standard globally
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set architecture flags based on platform
if(APPLE)
    find_package(OpenMP REQUIRED C Fortran)

    # Check if find_package fully resolved OpenMP (compiler knows where it is)
    if(NOT OpenMP_C_LIBRARIES OR NOT OpenMP_Fortran_LIBRARIES)
        # Homebrew libomp fallback — only if find_package didn't populate libraries
        find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
        if(BREW_EXECUTABLE)
            execute_process(
                COMMAND "${BREW_EXECUTABLE}" --prefix libomp
                OUTPUT_VARIABLE LIBOMP_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if(EXISTS "${LIBOMP_PREFIX}")
                message(STATUS "Using Homebrew libomp at ${LIBOMP_PREFIX}")
                include_directories(SYSTEM "${LIBOMP_PREFIX}/include")
                link_directories("${LIBOMP_PREFIX}/lib")
                set(LIBOMP_LINK_FLAGS "-L${LIBOMP_PREFIX}/lib -lomp")
            endif()
        endif()
    else()
        message(STATUS "OpenMP fully resolved by find_package (no Homebrew fallback needed)")
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion -Wno-external-argument-mismatch")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    if(LIBOMP_LINK_FLAGS)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBOMP_LINK_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LIBOMP_LINK_FLAGS}")
    endif()
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion -Wno-external-argument-mismatch")
endif()

enable_testing()

include(${CMAKE_SOURCE_DIR}/cmake/Dependencies.cmake)

set(CMAKE_BUILD_TYPE Release)

# Option to install KIM with system-wide symbolic link
option(INSTALL_KIM_SYMLINK "Create symbolic link 'kim' to KIM.x in /usr/local/bin" OFF)

# Set default install prefix if not specified
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Installation prefix" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Unified build: add subprojects for KiLCA, KIM, and QL-Balance
add_subdirectory(common/math)
add_subdirectory(common/logger)
add_subdirectory(common/hdf5_tools)
add_subdirectory(common/equil)
add_subdirectory(KiLCA)
add_subdirectory(KIM)
add_subdirectory(QL-Balance)

# Preprocessing tools
add_subdirectory(PreProc/fourier)

# Integration and regression tests (Python-based)
add_subdirectory(test)

# Install symbolic link for KIM if requested
if(INSTALL_KIM_SYMLINK)
    # Create a custom target to install the symlink
    add_custom_target(install_kim_symlink ALL
        COMMAND ${CMAKE_COMMAND} -E echo "Creating symbolic link for KIM..."
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_INSTALL_PREFIX}/bin/KIM.x
            ${CMAKE_BINARY_DIR}/kim_symlink_temp
        DEPENDS KIM.x
    )

    install(CODE "
        set(KIM_EXECUTABLE \"${CMAKE_INSTALL_PREFIX}/bin/KIM.x\")
        if(NOT EXISTS \${KIM_EXECUTABLE})
            message(WARNING \"KIM executable not found at \${KIM_EXECUTABLE}\")
        else()
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E create_symlink
                    \${KIM_EXECUTABLE}
                    /usr/local/bin/kim
                RESULT_VARIABLE symlink_result
                ERROR_VARIABLE symlink_error
            )
            if(NOT symlink_result EQUAL 0)
                message(WARNING \"Could not create symbolic link /usr/local/bin/kim. You may need sudo privileges. Error: \${symlink_error}\")
                message(STATUS \"Try running: sudo ln -sf \${KIM_EXECUTABLE} /usr/local/bin/kim\")
            else()
                message(STATUS \"Created symbolic link: /usr/local/bin/kim -> \${KIM_EXECUTABLE}\")
            endif()
        endif()
    ")
endif()
