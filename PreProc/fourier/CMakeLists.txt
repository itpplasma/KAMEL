# When built standalone, this sets up the project
# When built as subdirectory of KAMEL, these are already set
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.24)
    project(fouriermodes LANGUAGES Fortran)
    set(FOURIER_STANDALONE TRUE)
else()
    set(FOURIER_STANDALONE FALSE)
endif()

set(src ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(objs ${CMAKE_CURRENT_BINARY_DIR}/OBJS/)
set(lib ${CMAKE_CURRENT_BINARY_DIR}/lib/)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake")
include(Util)

# Only fetch libneo if not already available (prevents duplicate targets)
if(NOT TARGET LIBNEO::magfie)
    find_or_fetch(libneo)
    set(_FOURIER_LIBNEO_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/libneo")
else()
    set(_FOURIER_LIBNEO_BINARY_DIR "${KAMEL_LIBNEO_BINARY_DIR}")
endif()
include_directories(${_FOURIER_LIBNEO_BINARY_DIR}/include)

# Add netcdf library paths (needed because libneo's link_directories doesn't propagate)
find_program(NF_CONFIG "nf-config")
if(NF_CONFIG)
    execute_process(COMMAND nf-config --prefix
        OUTPUT_VARIABLE NETCDF_FORTRAN_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND nc-config --prefix
        OUTPUT_VARIABLE NETCDF_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    link_directories(${NETCDF_FORTRAN_PREFIX}/lib ${NETCDF_PREFIX}/lib)
endif()

# Source files - equilibrium tracing (mag, rhs1, RK4D) now comes from kamel_equil
set(source_files
    ${src}/rhs_flt.f90
    ${src}/change_q_sign.f90
    ${src}/checksizeff.f90
    ${src}/reload_amns.f90
    ${src}/rk4dc.f90
)

add_library(fouriermodes STATIC ${source_files})
set_target_properties(fouriermodes PROPERTIES
    LINKER_LANGUAGE Fortran
    OUTPUT_NAME fouriermodes
    ARCHIVE_OUTPUT_DIRECTORY ${lib}
    Fortran_MODULE_DIRECTORY ${objs}
)

# Link dependencies - kamel_equil provides mag, rhs1, RK4D for equilibrium tracing
if(TARGET kamel_equil)
    target_link_libraries(fouriermodes PUBLIC
        kamel_equil
        LIBNEO::magfie
    )
    target_include_directories(fouriermodes PUBLIC
        ${_FOURIER_LIBNEO_BINARY_DIR}/include
        ${KAMEL_EQUIL_INCLUDE_DIRS}
    )
else()
    # Standalone build - need to build kamel_equil separately or link directly to libneo
    message(FATAL_ERROR "kamel_equil target not found. Build fouriermodes as part of KAMEL or ensure kamel_equil is available.")
endif()

add_executable(fouriermodes.x ${src}/fouriermodes.f90)
target_link_libraries(fouriermodes.x PUBLIC fouriermodes)
set_target_properties(fouriermodes.x PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/install/bin
)
