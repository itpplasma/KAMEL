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

# QUADPACK (vendored from Netlib, public domain). Previously downloaded at
# configure time from netlib.org with no integrity check and an "if not exists"
# guard, so a single bad/partial download poisoned the CI build cache (the
# cached archive was missing dqagi_, breaking the link and never rebuilding).
# Vendoring the sources removes the network/cache fragility.
set(QUADPACK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/math/quadpack")

add_library(quadpack_netlib STATIC
    ${QUADPACK_DIR}/dqag.f
    ${QUADPACK_DIR}/dqage.f
    ${QUADPACK_DIR}/dqags.f
    ${QUADPACK_DIR}/dqagse.f
    ${QUADPACK_DIR}/dqagi.f
    ${QUADPACK_DIR}/dqagie.f
    ${QUADPACK_DIR}/dqelg.f
    ${QUADPACK_DIR}/dqpsrt.f
    ${QUADPACK_DIR}/dqk15.f
    ${QUADPACK_DIR}/dqk21.f
    ${QUADPACK_DIR}/dqk31.f
    ${QUADPACK_DIR}/dqk41.f
    ${QUADPACK_DIR}/dqk51.f
    ${QUADPACK_DIR}/dqk61.f
    ${QUADPACK_DIR}/dqk15i.f
)
set_target_properties(quadpack_netlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(KIM_lib PUBLIC quadpack_netlib)

# Optionally provide minimal XERROR compatibility if SLATEC doesn't provide it
option(KIM_XERROR_FALLBACK "Provide fallback XERROR for QUADPACK if SLATEC lacks it" ON)
if(KIM_XERROR_FALLBACK)
  target_sources(quadpack_netlib PRIVATE ${CMAKE_SOURCE_DIR}/KIM/src/math/quadpack_compat_xerror.f90)
endif()

# Provide a portable D1MACH for QUADPACK (avoids platform-specific Netlib d1mach variants)
target_sources(quadpack_netlib PRIVATE ${CMAKE_SOURCE_DIR}/KIM/src/math/quadpack_d1mach.f90)

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
# KIM_exe links KiLCA directly via slatec (which provides the math routines).
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
                        slatec
)
