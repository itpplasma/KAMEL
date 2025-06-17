include(ExternalProject)

if(NOT DEFINED GSL_INSTALL_DIR)
    set(GSL_INSTALL_DIR
        ${CMAKE_BINARY_DIR}/external-install/gsl
        CACHE PATH "Where to install GSL headers & libs")
endif()

ExternalProject_Add(EXTERNAL_gsl
    PREFIX ${CMAKE_BINARY_DIR}/download
    URL https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
    DOWNLOAD_NAME gsl-2.4.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    BUILD_IN_SOURCE TRUE

    CONFIGURE_COMMAND
        <SOURCE_DIR>/configure --prefix=${GSL_INSTALL_DIR}

    BUILD_COMMAND
        make -C <SOURCE_DIR> -j

    INSTALL_COMMAND
        make -C <SOURCE_DIR> install

    BUILD_BYPRODUCTS
        ${GSL_INSTALL_DIR}/lib/libgsl.a
        ${GSL_INSTALL_DIR}/lib/libgslcblas.a
)

add_library(gsl STATIC IMPORTED GLOBAL)
set_target_properties(gsl PROPERTIES
    IMPORTED_LOCATION
    "${GSL_INSTALL_DIR}/lib/libgsl.a"
)
add_dependencies(gsl EXTERNAL_gsl)

add_library(gslcblas STATIC IMPORTED GLOBAL)
set_target_properties(gslcblas PROPERTIES
    IMPORTED_LOCATION
    "${GSL_INSTALL_DIR}/lib/libgslcblas.a"
)
add_dependencies(gslcblas EXTERNAL_gsl)
