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
                        "${KIM_flr2}"
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

# SuiteSparse/UMFPACK come SOLELY from the fetched-static SuiteSparse v7.8.3 via the
# `sparse` target linked above (Dependencies.cmake -> SuiteSparse::umfpack{,_wrappers},
# built with -DDLONG -DZLONG). Previously a find_package(SuiteSparse)+find_package(Umfpack)
# here ALSO linked the system (Homebrew) UMFPACK, so the binary loaded two conflicting
# UMFPACK builds at once (static 7.8.3 + dylib 6.3.7). The duplicate SuiteSparse_config
# global state corrupted every complex sparse solve (umf4zsym -> UMFPACK_ERROR_invalid_matrix,
# solutions Inf/0), which silently made the global electrostatic Poisson solve return Phi==0.
# The umf4_f77*.c wrappers' umfpack.h and the UMFPACK symbols both propagate transitively from
# SuiteSparse::umfpack via the `sparse` target, so no separate find_package is needed.

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

target_link_libraries(KIM_lib PUBLIC
                                kamel_hdf5_tools
                                kamel_logger
                                kamel_equil
                                cerf
                                kamel_util
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
