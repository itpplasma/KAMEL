# Built as part of KAMEL only: kamel_equil and libneo come from the top level.

set(src ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(objs ${CMAKE_CURRENT_BINARY_DIR}/OBJS/)
set(lib ${CMAKE_CURRENT_BINARY_DIR}/lib/)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

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

# kamel_equil provides mag, rhs1 and RK4D for equilibrium tracing.
target_link_libraries(fouriermodes PUBLIC
    kamel_equil
    LIBNEO::magfie
)
target_include_directories(fouriermodes PUBLIC
    ${KAMEL_EQUIL_INCLUDE_DIRS}
)

add_executable(fouriermodes.x ${src}/fouriermodes.f90)
target_link_libraries(fouriermodes.x PUBLIC fouriermodes)
set_target_properties(fouriermodes.x PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/install/bin
)
