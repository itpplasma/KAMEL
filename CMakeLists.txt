cmake_minimum_required(VERSION 3.24)

# Set compilers BEFORE project() - only if not already set by user
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    if(NOT DEFINED CMAKE_C_COMPILER)
        # Search for LLVM clang (needed for OpenMP support on macOS)
        find_program(LLVM_CLANG clang
            PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin
            NO_DEFAULT_PATH)
        if(LLVM_CLANG)
            set(CMAKE_C_COMPILER "${LLVM_CLANG}")
        endif()
    endif()
    if(NOT DEFINED CMAKE_CXX_COMPILER)
        find_program(LLVM_CLANGXX clang++
            PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin
            NO_DEFAULT_PATH)
        if(LLVM_CLANGXX)
            set(CMAKE_CXX_COMPILER "${LLVM_CLANGXX}")
        endif()
    endif()
    set(CMAKE_PREFIX_PATH /opt/homebrew /usr/local)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    if(NOT DEFINED CMAKE_C_COMPILER)
        set(CMAKE_C_COMPILER gcc)
    endif()
    if(NOT DEFINED CMAKE_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER g++)
    endif()
endif()

if(APPLE AND NOT DEFINED HDF5_DIR)
    # locate brew
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
            message(STATUS "HDF5_DIR set to ${HDF5_DIR}")
        endif()
    endif()
endif()

project(KAMEL LANGUAGES C CXX Fortran)

# Set C++ standard globally
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED C Fortran)

# Set architecture flags based on platform
if(APPLE)
    # On macOS, use libomp from Homebrew for OpenMP support
    # find_package(OpenMP REQUIRED)
    find_package(OpenMP REQUIRED C Fortran)

    # Locate libomp from Homebrew
    find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
    if(BREW_EXECUTABLE)
        execute_process(
            COMMAND "${BREW_EXECUTABLE}" --prefix libomp
            OUTPUT_VARIABLE LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(EXISTS "${LIBOMP_PREFIX}")
            include_directories(SYSTEM "${LIBOMP_PREFIX}/include")
            link_directories("${LIBOMP_PREFIX}/lib")
            message(STATUS "Found libomp at ${LIBOMP_PREFIX}")
        endif()
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion -Wno-external-argument-mismatch")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LIBOMP_PREFIX}/lib -lomp")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -L${LIBOMP_PREFIX}/lib -lomp")
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
add_subdirectory(common/hdf5_tools)
add_subdirectory(common/equil)
add_subdirectory(KiLCA)
add_subdirectory(KIM)
add_subdirectory(QL-Balance)

# Preprocessing tools
add_subdirectory(PreProc/fourier)

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
