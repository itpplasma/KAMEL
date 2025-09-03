cmake_minimum_required(VERSION 3.24)

set(CMAKE_Fortran_COMPILER "mpif90")

project(KAMEL LANGUAGES C CXX Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion")

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
add_subdirectory(KiLCA)
add_subdirectory(KIM)
add_subdirectory(QL-Balance)

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
