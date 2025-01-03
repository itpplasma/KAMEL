# Use ExternalProject instead of FetchContent to allow for custom build options
include(ExternalProject)

ExternalProject_Add(
    SuiteSparse
    DOWNLOAD_EXTRACT_TIMESTAMP
    URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.8.3.tar.gz
    URL_HASH MD5=242e38ecfc8a3e3aa6b7d8d44849c5cf
    CMAKE_ARGS -GNinja -DUMFPACK_USE_CHOLMOD=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target UMFPACK/install
                && ${CMAKE_COMMAND} --build <BINARY_DIR> --target AMD/install
                && ${CMAKE_COMMAND} --build <BINARY_DIR> --target SuiteSparse_config/install
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/install/lib/libumfpack.a
        ${CMAKE_BINARY_DIR}/install/lib/libsuitesparseconfig.a
        ${CMAKE_BINARY_DIR}/install/lib/libamd.a
        ${CMAKE_BINARY_DIR}/install/include/suitesparse/umfpack.h
        <SOURCE_DIR>UMFPACK/Demo/umf4_f77wrapper.c
        <SOURCE_DIR>UMFPACK/Demo/umf4_f77zwrapper.c
)

set(UMFPACK_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libumfpack.a)
set(SUITESPARSE_CONFIG_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libsuitesparseconfig.a)
set(AMD_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libamd.a)
set(UMFPACK_INCLUDE_DIR ${CMAKE_BINARY_DIR}/install/include/suitesparse)

ExternalProject_Get_Property(SuiteSparse SOURCE_DIR)
set(UMFPACK_WRAPPERS
    ${SOURCE_DIR}/UMFPACK/Demo/umf4_f77wrapper.c
    ${SOURCE_DIR}/UMFPACK/Demo/umf4_f77zwrapper.c
)
set_source_files_properties(${SOURCE_DIR}/UMFPACK/Demo/umf4_f77wrapper.c PROPERTIES COMPILE_FLAGS "-I${UMFPACK_INCLUDE_DIR} -DDLONG" GENERATED TRUE)
set_source_files_properties(${SOURCE_DIR}/UMFPACK/Demo/umf4_f77zwrapper.c PROPERTIES COMPILE_FLAGS "-I${UMFPACK_INCLUDE_DIR} -DZLONG" GENERATED TRUE)

add_library(SuiteSparse::suitesparse_config SHARED IMPORTED)
set_target_properties(SuiteSparse::suitesparse_config PROPERTIES
    IMPORTED_LOCATION ${SUITESPARSE_CONFIG_LIBRARY_PATH}
)

add_library(SuiteSparse::amd SHARED IMPORTED)
set_target_properties(SuiteSparse::amd PROPERTIES
    IMPORTED_LOCATION ${AMD_LIBRARY_PATH}
)

add_library(SuiteSparse::umfpack SHARED IMPORTED)
set_target_properties(SuiteSparse::umfpack PROPERTIES
    IMPORTED_LOCATION ${UMFPACK_LIBRARY_PATH}
)
