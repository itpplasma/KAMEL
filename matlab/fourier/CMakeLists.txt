cmake_minimum_required(VERSION 3.24)

project(fouriermodes)
enable_language(Fortran)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    # Linux specific configurations or commands
    message("Building on Linux")
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_CXX_COMPILER "g++")
    set(CMAKE_Fortran_COMPILER "mpif90")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system might not be implemented.")
endif()

set(src ${PROJECT_SOURCE_DIR}/src)
set(objs ${PROJECT_BINARY_DIR}/OBJS/)
set(lib ${PROJECT_BINARY_DIR}/lib/)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

set(source_files 
    ${src}/bdivfree.f90
    ${src}/change_q_sign.f90
    ${src}/checksizeff.f90
    ${src}/field_divB0.f90
    ${src}/mag.f90
    ${src}reload_amns.f90
    ${src}/rhs_flt.f90
    ${src}/rk4dc.f90
    ${src}spline5_RZ.f90)

add_library(fouriermodes STATIC ${source_files})
set_target_properties(fouriermodes PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME fouriermodes
    ARCHIVE_OUTPUT_DIRECTORY ${lib}
    Fortran_MODULE_DIRECTORY ${objs})

add_executable(fouriermodes.x ${src}/fouriermodes.f90)