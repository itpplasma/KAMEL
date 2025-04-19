# Use ExternalProject instead of FetchContent to allow for custom build options
include(ExternalProject)

find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS C)

set(UMFPACK_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libumfpack${CMAKE_STATIC_LIBRARY_SUFFIX})
set(SUITESPARSE_CONFIG_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX})
set(AMD_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX})
set(SUITESPARSE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/install/include/suitesparse)
set(UMFPACK_INCLUDE_DIR ${CMAKE_BINARY_DIR}/install/include/suitesparse)

execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SUITESPARSE_INCLUDE_DIR}
)

ExternalProject_Add(
    SuiteSparse
    DOWNLOAD_EXTRACT_TIMESTAMP
    URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.8.3.tar.gz
    URL_HASH MD5=242e38ecfc8a3e3aa6b7d8d44849c5cf
    CMAKE_ARGS
        -GNinja
	    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        #-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
        #-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        #-DBUILD_SHARED_LIBS=OFF
        -DUMFPACK_USE_CHOLMOD=OFF
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target UMFPACK/install
                && ${CMAKE_COMMAND} --build <BINARY_DIR> --target AMD/install
                && ${CMAKE_COMMAND} --build <BINARY_DIR> --target SuiteSparse_config/install
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS
        ${UMFPACK_LIBRARY_PATH}
        ${SUITESPARSE_CONFIG_LIBRARY_PATH}
        ${AMD_LIBRARY_PATH}
        ${SUITESPARSE_INCLUDE_DIR}
        <SOURCE_DIR>UMFPACK/Demo/umf4_f77wrapper.c
        <SOURCE_DIR>UMFPACK/Demo/umf4_f77zwrapper.c
)

set(UMFPACK_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libumfpack${CMAKE_STATIC_LIBRARY_SUFFIX})
set(SUITESPARSE_CONFIG_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX})
set(AMD_LIBRARY_PATH ${CMAKE_BINARY_DIR}/install/lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX})
set(UMFPACK_INCLUDE_DIR ${CMAKE_BINARY_DIR}/install/include/suitesparse)

ExternalProject_Get_Property(SuiteSparse SOURCE_DIR)

add_library(SuiteSparse::suitesparse_config STATIC IMPORTED)
set_target_properties(SuiteSparse::suitesparse_config PROPERTIES
    IMPORTED_LOCATION ${SUITESPARSE_CONFIG_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES "BLAS::BLAS;OpenMP::OpenMP_C"
)

add_library(SuiteSparse::amd STATIC IMPORTED)
set_target_properties(SuiteSparse::amd PROPERTIES
    IMPORTED_LOCATION ${AMD_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES SuiteSparse::suitesparse_config
)

add_library(SuiteSparse::umfpack STATIC IMPORTED)
set_target_properties(SuiteSparse::umfpack PROPERTIES
    IMPORTED_LOCATION ${UMFPACK_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES SuiteSparse::amd
)

set(UMFPACK_WRAPPERS
    ${SOURCE_DIR}/UMFPACK/Demo/umf4_f77wrapper.c
    ${SOURCE_DIR}/UMFPACK/Demo/umf4_f77zwrapper.c
)
set_source_files_properties(${UMFPACK_WRAPPERS} PROPERTIES GENERATED TRUE)

add_library(umfpack_wrappers STATIC ${UMFPACK_WRAPPERS})
add_library(SuiteSparse::umfpack_wrappers ALIAS umfpack_wrappers)
add_dependencies(umfpack_wrappers SuiteSparse)
target_link_libraries(umfpack_wrappers PUBLIC BLAS::BLAS OpenMP::OpenMP_C)
target_include_directories(umfpack_wrappers PRIVATE ${SUITESPARSE_INCLUDE_DIR})
target_compile_options(umfpack_wrappers PRIVATE -DDLONG -DZLONG)
