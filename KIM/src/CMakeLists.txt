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
                        "${external_zeal}"
)

set_target_properties(KIM_lib PROPERTIES 
                        OUTPUT_NAME KIM
                        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/install/lib/")

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

## Fetch QUADPACK directly from Netlib (avoid Zeal's copies)
set(QUADPACK_DIR "${CMAKE_BINARY_DIR}/external/quadpack")
file(MAKE_DIRECTORY "${QUADPACK_DIR}")

set(QUADPACK_FILES
    dqag.f
    dqage.f
    dqags.f
    dqagse.f
    dqagi.f
    dqagie.f
    dqelg.f
    dqpsrt.f
    dqk15.f
    dqk21.f
    dqk31.f
    dqk41.f
    dqk51.f
    dqk61.f
    dqk15i.f
)

foreach(f ${QUADPACK_FILES})
    if(NOT EXISTS "${QUADPACK_DIR}/${f}")
        message(STATUS "Downloading QUADPACK source ${f} from netlib.org")
        file(DOWNLOAD "https://netlib.org/quadpack/${f}" "${QUADPACK_DIR}/${f}" SHOW_PROGRESS)
    endif()
endforeach()

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
set(KiLCA_lib_path ${CMAKE_SOURCE_DIR}/build/install/lib/libKiLCA_Lib_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit.a)
target_link_libraries(KIM_lib PUBLIC "${KiLCA_lib_path}")
target_include_directories(KIM_lib PUBLIC ${CMAKE_SOURCE_DIR}/KiLCA/math/)
target_include_directories(KIM_lib PUBLIC ${KILCA_PATH}/math/)

add_executable(KIM_exe kim_main.f90)
set_target_properties(KIM_exe PROPERTIES 
                        OUTPUT_NAME KIM.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/install/bin/")
target_link_libraries(KIM_exe
                        KIM_lib 
                        lapack 
                        cerf 
                        slatec 
                        ddeabm
)

# Optional: build QUADPACK/RKF45 integration comparison test
add_executable(test_integration_methods ${CMAKE_SOURCE_DIR}/KIM/tests/test_integration_methods.f90)
set_target_properties(test_integration_methods PROPERTIES 
                        OUTPUT_NAME test_integration_methods.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_integration_methods KIM_lib lapack cerf slatec ddeabm)
add_test(NAME test_integration_methods
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_integration_methods.x)

# QUADPACK DQAGI sanity test (infinite integral)
add_executable(test_quadpack_qagi ${CMAKE_SOURCE_DIR}/KIM/tests/test_quadpack_qagi.f90)
set_target_properties(test_quadpack_qagi PROPERTIES 
                        OUTPUT_NAME test_quadpack_qagi.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_quadpack_qagi KIM_lib lapack cerf slatec ddeabm)
add_test(NAME test_quadpack_qagi
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_quadpack_qagi.x)
