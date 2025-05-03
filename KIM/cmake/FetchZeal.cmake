include(FetchContent)

find_program(GFORTRAN_EXECUTABLE gfortran REQUIRED
    DOC "Plain gfortran for building Zeal"
)

FetchContent_Declare(zeal
    URL https://elsevier.digitalcommonsdata.com/public-files/datasets/yc9vv7rwyj/files/be55b15b-6d9c-4f0e-b5f7-09b300cfb806/file_downloaded
    URL_HASH MD5=5450326083c6114b026fd765a9733928
)

FetchContent_Populate(zeal)

set(ZEAL_BUILD_STAMP ${CMAKE_BINARY_DIR}/zeal-build.stamp)

add_custom_command(
    OUTPUT ${ZEAL_BUILD_STAMP}
    COMMENT "Building Zeal (lapack + main)..."
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${GFORTRAN_EXECUTABLE} ${zeal_SOURCE_DIR}/f90
    COMMAND ${CMAKE_COMMAND} -E env PATH=<SOURCE_DIR>:$ENV{PATH} make lapack zeal -j
    COMMAND ${CMAKE_COMMAND} -E touch ${ZEAL_BUILD_STAMP}
    WORKING_DIRECTORY ${zeal_SOURCE_DIR}

    # DEPENDS ${GFORTRAN_EXECUTABLE}
)
