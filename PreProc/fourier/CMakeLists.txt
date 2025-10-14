cmake_minimum_required(VERSION 3.24)

project(fouriermodes LANGUAGES Fortran)

set(src ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(objs ${CMAKE_CURRENT_BINARY_DIR}/OBJS/)
set(lib ${CMAKE_CURRENT_BINARY_DIR}/lib/)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake")
include(Util)
find_or_fetch(libneo)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/libneo/include)

set(source_files
    ${src}/rhs_flt.f90
    ${src}/change_q_sign.f90
    ${src}/checksizeff.f90
    ${src}/mag.f90
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

target_link_libraries(fouriermodes PUBLIC
    LIBNEO::magfie
)

add_executable(fouriermodes.x ${src}/fouriermodes.f90)
target_link_libraries(fouriermodes.x PUBLIC fouriermodes)
set_target_properties(fouriermodes.x PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

target_include_directories(fouriermodes PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/libneo/include
)
