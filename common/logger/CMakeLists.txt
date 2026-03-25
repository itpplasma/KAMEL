set(LOGGER_MODULE_DIR ${CMAKE_BINARY_DIR}/modules/logger)
file(MAKE_DIRECTORY ${LOGGER_MODULE_DIR})

add_library(kamel_logger STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/logger_m.f90
)

set_target_properties(kamel_logger PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${LOGGER_MODULE_DIR}
    LINKER_LANGUAGE Fortran
)

target_include_directories(kamel_logger
    PUBLIC ${LOGGER_MODULE_DIR}
)
