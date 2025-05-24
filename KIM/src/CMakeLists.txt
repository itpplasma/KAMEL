include(CMakeSources.in)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_library(KIM_lib STATIC
    ${KIM_setup}
    ${KIM_general}
    ${KIM_poisson}
    ${KIM_debye}
    ${KIM_grid}
    ${KIM_util}
    ${KIM_kernels}
    ${KIM_background_equilibrium}
    ${KIM_math}
    ${KIM_dispersion}
    ${KIM_WKB_dispersion}
    ${external_zeal}
)
set_target_properties(KIM_lib PROPERTIES 
                        OUTPUT_NAME KIM_lib.a
                        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/")

# TODO: unify this: either use the system library everywhere, or use the cloned version
# KiLCA uses the cloned libraries for now, QL-Balance the system libraries
find_package(GSL REQUIRED)
target_link_libraries(KIM_lib PUBLIC GSL::gsl GSL::gslcblas)
target_include_directories(KIM_lib PUBLIC "${GSL_INCLUDE_DIRS}")

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


add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/")
target_link_libraries(KIM_exe
                        KIM_lib 
                        lapack 
                        cerf 
                        slatec 
                        ddeabm
)

