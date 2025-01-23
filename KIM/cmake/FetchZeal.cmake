include(ExternalProject)

ExternalProject_Add(
    Zeal
    DOWNLOAD_EXTRACT_TIMESTAMP
    URL https://elsevier.digitalcommonsdata.com/public-files/datasets/yc9vv7rwyj/files/be55b15b-6d9c-4f0e-b5f7-09b300cfb806/file_downloaded
    #CMAKE_ARGS -GNinja
    INSTALL_COMMAND ""
    BUILD_COMMAND make F90=gfortran LIBS=-llapack
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE TRUE
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/install/lib/libzeal${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/install/lib/libzeal${CMAKE_SHARED_LIBRARY_SUFFIX}
)

add_library(Zeal::zeal SHARED IMPORTED)
set_target_properties(Zeal::zeal PROPERTIES
    IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/install/lib/libzeal${CMAKE_SHARED_LIBRARY_SUFFIX}
)