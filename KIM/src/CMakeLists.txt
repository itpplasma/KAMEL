include(CMakeSources.in)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


add_library(KIM_lib STATIC "${KIM_setup}"
                        "${KIM_grid}"
                        "${KIM_background_equilibrium}"
                        "${KIM_general}"
                        "${KIM_poisson}"
                        "${KIM_util}"
                        #"${KIM_kernels}"
                        "${KIM_math}"
                        "${KIM_fokkerplanck}"
                        "${KIM_dispersion}"
                        "${KIM_WKB_dispersion}"
                        "${KIM_asymptotics}"
                        "${KIM_diagnostics}"
                        "${KIM_electromagnetic}"
)

set_target_properties(KIM_lib PROPERTIES
                        OUTPUT_NAME KIM
                        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/install/lib/")

# Use the shared sparse target instead of bundling a local copy
add_dependencies(KIM_lib sparse)
target_link_libraries(KIM_lib PUBLIC sparse)
target_include_directories(KIM_lib PUBLIC "${CMAKE_BINARY_DIR}/OBJS/sparse/")

# fortnum provides the special functions (bessel_in) used by FLR2_asymptotics.
target_link_libraries(KIM_lib PUBLIC fortnum)

find_package(LAPACK REQUIRED)
if("${LAPACK_FOUND}")
    message(" ====== LAPACK was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${LAPACK_LIBRARIES}")
    message("LAPACK libraries: ${LAPACK_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${LAPACK_INCLUDE_DIRS}")
    message("${LAPACK_INCLUDE_DIRS}")
endif()

find_package(SuiteSparse REQUIRED)
if("${SuiteSparse_FOUND}")
    message(" ====== SuiteSparse was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${SuiteSparse_LIBRARIES}")
    message("SuiteSparse libraries: ${SUITESPARSE_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${SuiteSparse_INCLUDE_DIRS}")
    message("${SUITESPARSE_INCLUDE_DIRS}")
endif()

find_package(Umfpack REQUIRED)
if("${UMFPACK_FOUND}")
    message(" ====== Umfpack was found! ======")
    target_link_libraries(KIM_lib PUBLIC "${UMFPACK_LIBRARIES}")
    message("UMFPACK libraries: ${UMFPACK_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${UMFPACK_INCLUDES}")
    message("${UMFPACK_INCLUDES}")
endif()

find_library(UMFPACK_LIBRARY umfpack)

find_package(SuperLU REQUIRED)
if("${SUPERLU_FOUND}")
    message(" ====== SuperLU was found! =====")
    target_link_libraries(KIM_lib PUBLIC "${SUPERLU_LIBRARIES}")
    message("SuperLu libraries: ${SUPERLU_LIBRARIES}")
    target_include_directories(KIM_lib PUBLIC "${SUPERLU_INCLUDES}")
    message("${SUPERLU_INCLUDES}")
endif()

set(KILCA_PATH ${CMAKE_SOURCE_DIR}/KiLCA)
# KiLCA is NOT linked to KIM_lib: KIM does not use any KiLCA symbols.
# Previously PUBLIC, this caused duplicate library linking when QL-Balance
# linked both KIM_lib and kilca_lib, corrupting symbol resolution order
# and breaking KiLCA's wave solver (besselj failures, vacuum stitching errors).
# KIM_exe links KiLCA directly via fortnum_amos_compat (complex Bessel ABI).
target_include_directories(KIM_lib PUBLIC ${KILCA_PATH}/math/)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(" ====== HDF5 was found! ======")

target_link_libraries(KIM_lib PUBLIC
                                hdf5::hdf5
                                hdf5::hdf5_fortran
                                hdf5::hdf5_hl
                                hdf5::hdf5_hl_fortran
                                kamel_hdf5_tools
                                kamel_logger
                                kamel_equil
                                cerf
)

add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES
                        OUTPUT_NAME KIM.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/install/bin/")
target_link_libraries(KIM_exe
                        KIM_lib
                        kilca_lib
                        lapack
                        fortnum_amos_compat
)
