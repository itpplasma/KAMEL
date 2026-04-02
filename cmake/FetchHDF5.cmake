# FetchHDF5.cmake — Find or auto-fetch HDF5 with Fortran and HL support
#
# If HDF5 is found on the system (without MPI/parallel requirement), use it.
# Otherwise, download and build HDF5 from source (serial, with Fortran and HL bindings).
#
# Users can pass -DHDF5_DIR=<path> to point to a specific HDF5 installation.

# If user explicitly provided HDF5_DIR, trust it and try to use it.
# Otherwise, check if the system HDF5 config requires parallel/MPI (which we don't use).
# If it does, skip it and build from source to avoid MPI dependency.
if(NOT DEFINED HDF5_DIR)
    # Look for hdf5-config.cmake on CMAKE_PREFIX_PATH / system paths
    find_file(_hdf5_config_file hdf5-config.cmake
        PATHS ${CMAKE_PREFIX_PATH}
        PATH_SUFFIXES lib/cmake/hdf5 share/cmake/hdf5
        NO_DEFAULT_PATH)
    if(NOT _hdf5_config_file)
        find_file(_hdf5_config_file hdf5-config.cmake
            PATH_SUFFIXES lib/cmake/hdf5 share/cmake/hdf5)
    endif()
    if(_hdf5_config_file)
        file(READ "${_hdf5_config_file}" _hdf5_config_content)
        string(FIND "${_hdf5_config_content}" "ENABLE_PARALLEL       ON" _parallel_pos)
        if(_parallel_pos GREATER -1)
            message(STATUS "System HDF5 requires MPI (parallel) — skipping, will build from source")
            set(_skip_system_hdf5 TRUE)
        endif()
    endif()
endif()

if(NOT _skip_system_hdf5)
    find_package(HDF5 COMPONENTS C Fortran HL QUIET)
endif()

if(HDF5_FOUND)
    message(STATUS "HDF5 found on system: ${HDF5_VERSION}")
else()
    message(STATUS "HDF5 not found on system — downloading and building from source...")

    include(ExternalProject)

    set(HDF5_VERSION "1.14.6")
    set(HDF5_URL "https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz")
    set(HDF5_INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/hdf5-install")

    ExternalProject_Add(EXTERNAL_hdf5
        URL ${HDF5_URL}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${HDF5_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            -DHDF5_BUILD_FORTRAN=ON
            -DHDF5_BUILD_HL_LIB=ON
            -DHDF5_ENABLE_PARALLEL=OFF
            -DHDF5_BUILD_TOOLS=OFF
            -DHDF5_BUILD_EXAMPLES=OFF
            -DBUILD_TESTING=OFF
            -DBUILD_SHARED_LIBS=OFF
        BUILD_BYPRODUCTS
            ${HDF5_INSTALL_DIR}/lib/libhdf5.a
            ${HDF5_INSTALL_DIR}/lib/libhdf5_fortran.a
            ${HDF5_INSTALL_DIR}/lib/libhdf5_hl.a
            ${HDF5_INSTALL_DIR}/lib/libhdf5_hl_fortran.a
    )

    # Create imported targets matching find_package(HDF5) interface
    file(MAKE_DIRECTORY ${HDF5_INSTALL_DIR}/include)

    add_library(hdf5::hdf5 STATIC IMPORTED GLOBAL)
    set_target_properties(hdf5::hdf5 PROPERTIES
        IMPORTED_LOCATION ${HDF5_INSTALL_DIR}/lib/libhdf5.a
        INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INSTALL_DIR}/include
    )
    add_dependencies(hdf5::hdf5 EXTERNAL_hdf5)

    add_library(hdf5::hdf5_fortran STATIC IMPORTED GLOBAL)
    set_target_properties(hdf5::hdf5_fortran PROPERTIES
        IMPORTED_LOCATION ${HDF5_INSTALL_DIR}/lib/libhdf5_fortran.a
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INSTALL_DIR}/include;${HDF5_INSTALL_DIR}/include/static"
    )
    add_dependencies(hdf5::hdf5_fortran EXTERNAL_hdf5)

    add_library(hdf5::hdf5_hl STATIC IMPORTED GLOBAL)
    set_target_properties(hdf5::hdf5_hl PROPERTIES
        IMPORTED_LOCATION ${HDF5_INSTALL_DIR}/lib/libhdf5_hl.a
        INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INSTALL_DIR}/include
    )
    add_dependencies(hdf5::hdf5_hl EXTERNAL_hdf5)

    add_library(hdf5::hdf5_hl_fortran STATIC IMPORTED GLOBAL)
    set_target_properties(hdf5::hdf5_hl_fortran PROPERTIES
        IMPORTED_LOCATION ${HDF5_INSTALL_DIR}/lib/libhdf5_hl_fortran.a
        INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INSTALL_DIR}/include;${HDF5_INSTALL_DIR}/include/static"
    )
    add_dependencies(hdf5::hdf5_hl_fortran EXTERNAL_hdf5)

    # Fortran module directory
    target_include_directories(hdf5::hdf5_fortran INTERFACE
        ${HDF5_INSTALL_DIR}/include/static
    )
    target_include_directories(hdf5::hdf5_hl_fortran INTERFACE
        ${HDF5_INSTALL_DIR}/include/static
    )

    # Mark as found for downstream
    set(HDF5_FOUND TRUE)
    set(HDF5_IS_FETCHED TRUE)
    message(STATUS "HDF5 ${HDF5_VERSION} will be built from source")
endif()
