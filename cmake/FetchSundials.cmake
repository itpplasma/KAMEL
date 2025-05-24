include(ExternalProject)

if(NOT DEFINED SUNDIALS_INSTALL_DIR)
    set(SUNDIALS_INSTALL_DIR
        ${CMAKE_BINARY_DIR}/external-install/sundials
        CACHE PATH "Where Sundials will be installed")
endif()

ExternalProject_Add(sundials # needs to be named like this due to libneo depending on it and finding it via a path
    PREFIX ${CMAKE_BINARY_DIR}/download
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG v5.7.0
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    BUILD_IN_SOURCE TRUE

    CONFIGURE_COMMAND
        ${CMAKE_COMMAND} -S <SOURCE_DIR> -B <SOURCE_DIR>/build
        -DCMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_DIR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}

    BUILD_COMMAND
        ${CMAKE_COMMAND} --build <SOURCE_DIR>/build -- -j

    INSTALL_COMMAND
        ${CMAKE_COMMAND} --build <SOURCE_DIR>/build --target install

    BUILD_BYPRODUCTS
        ${SUNDIALS_INSTALL_DIR}/lib/libsundials_cvode.a
        ${SUNDIALS_INSTALL_DIR}/lib/libsundials_nvecserial.a
)

add_library(cvode STATIC IMPORTED GLOBAL)
set_target_properties(cvode PROPERTIES
    IMPORTED_LOCATION
    "${SUNDIALS_INSTALL_DIR}/lib/libsundials_cvode.a"
)
add_dependencies(cvode sundials)

add_library(nvecserial STATIC IMPORTED GLOBAL)
set_target_properties(nvecserial PROPERTIES
    IMPORTED_LOCATION
    "${SUNDIALS_INSTALL_DIR}/lib/libsundials_nvecserial.a"
)
add_dependencies(nvecserial sundials)
