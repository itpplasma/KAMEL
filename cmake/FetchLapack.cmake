# First try to find system LAPACK/BLAS installation.
# On macOS with flang, the top-level CMakeLists.txt sets BLA_VENDOR=OpenBLAS
# and adds OpenBLAS to CMAKE_LIBRARY_PATH so FindBLAS can find it.
find_package(LAPACK QUIET)
find_package(BLAS QUIET)

if(LAPACK_FOUND AND BLAS_FOUND)
    message(STATUS "Using system LAPACK/BLAS libraries")
    message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
    message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")

    # Create interface libraries that link to system libraries
    add_library(lapack INTERFACE)
    target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})

    add_library(blas INTERFACE)
    target_link_libraries(blas INTERFACE ${BLAS_LIBRARIES})

    # Create a dummy tmg library for compatibility
    add_library(tmg INTERFACE)

else()
    message(STATUS "System LAPACK/BLAS not found, building from source")

    include(ExternalProject)

    if(NOT DEFINED LAPACK_INSTALL_DIR)
        set(LAPACK_INSTALL_DIR ${CMAKE_BINARY_DIR}/external-install/lapack CACHE PATH
            "Where to install LAPACK .a files")
    endif()

    ExternalProject_Add(EXTERNAL_lapack
        PREFIX ${CMAKE_BINARY_DIR}/download
        URL http://www.netlib.org/lapack/lapack-3.2.1.tgz
        URL_HASH MD5=a3202a4f9e2f15ffd05d15dab4ac7857
        DOWNLOAD_NAME lapack-3.2.1.tgz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        BUILD_IN_SOURCE TRUE

        CONFIGURE_COMMAND "" # no configure step

        BUILD_COMMAND
            cp INSTALL/make.inc.gfortran make.inc &&
            make -j lib blaslib tmglib

        BUILD_BYPRODUCTS
            <SOURCE_DIR>/lapack_LINUX.a
            <SOURCE_DIR>/blas_LINUX.a
            <SOURCE_DIR>/tmglib_LINUX.a

        INSTALL_COMMAND
            ${CMAKE_COMMAND} -E make_directory ${LAPACK_INSTALL_DIR}/lib &&
            ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/lapack_LINUX.a ${LAPACK_INSTALL_DIR}/lib/ &&
            ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/blas_LINUX.a ${LAPACK_INSTALL_DIR}/lib/ &&
            ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/tmglib_LINUX.a ${LAPACK_INSTALL_DIR}/lib/
    )

    ExternalProject_Get_Property(EXTERNAL_lapack source_dir)

    add_library(lapack STATIC IMPORTED)
    set_target_properties(lapack PROPERTIES
        IMPORTED_LOCATION "${source_dir}/lapack_LINUX.a"
    )
    add_dependencies(lapack EXTERNAL_lapack)

    add_library(blas STATIC IMPORTED)
    set_target_properties(blas PROPERTIES
        IMPORTED_LOCATION "${source_dir}/blas_LINUX.a"
    )
    add_dependencies(blas EXTERNAL_lapack)

    add_library(tmg STATIC IMPORTED)
    set_target_properties(tmg PROPERTIES
        IMPORTED_LOCATION "${source_dir}/tmglib_LINUX.a"
    )
    add_dependencies(tmg EXTERNAL_lapack)

endif()
