cmake_minimum_required (VERSION 3.12)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    #set(CMAKE_C_COMPILER gcc-14)
    set(CMAKE_C_COMPILER clang)
    #set(CMAKE_CXX_COMPILER g++-14)
    set(CMAKE_CXX_COMPILER clang++)
    # flag does not work for macos for some reason
    set(CMAKE_CXX_OPT_FLAGS "-O3")
    set(CMAKE_F_OPT_FLAGS "-O3")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux" AND
        CMAKE_SYSTEM_PROCESSOR MATCHES "(x86|i386|i686|amd64|x86_64)")
    # Linux specific configurations or commands
    message("Building on Linux")
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    set(CMAKE_CXX_OPT_FLAGS "-O3 -march=x86-64-v2 -mtune=generic")
    set(CMAKE_F_OPT_FLAGS "-O3 -march=x86-64-v2 -mtune=generic")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system. Might not be implemented.")
endif()

set(CMAKE_Fortran_COMPILER gfortran)

project(KiLCA)
enable_language(CXX Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/mod)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=generic ${CMAKE_CXX_OPT_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2 -mtune=generic ${CMAKE_F_OPT_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wshadow -mtune=generic")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O2 -Wall -Wextra -Wshadow -Wtabs -mtune=generic")

# --- Options ---
set(drift) # empty, "_drift"
set(coll "_gen") # empty, "_wc_n", "_wc_ne", "_wc_nep", "_nc", "_gen"
set(interp "_pol") # empty, "_pol"
set(Jac "_J0") # empty, "_J0"
add_compile_definitions(
    DEBUG_FLAG=0
    CALC_EIGEN_DECOMPOSITION=0
    SORT_DISPERSION_PROFILES=0
    USE_JACOBIAN_IN_ODE_SOLVER=1
    USE_SPLINES_IN_RHS_EVALUATION=0
)

# --- Dependencies ---
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(FetchGSL)
include(FetchLapack)
include(FetchSundials)

# include paths:
set(libs_root_path ${CMAKE_BINARY_DIR}/external-install)
set(sundials_include_path ${libs_root_path}/sundials/include)
set(gsl_include_path ${libs_root_path}/gsl/include)

file(GLOB libslatec_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/math/bessel/slatec/*.f")
add_library(slatec STATIC "${libslatec_files}")

set(EXTERNAL_LIBS slatec lapack blas cvode nvecserial gsl gslcblas)

include_directories(
    "${gsl_include_path}"
    "${sundials_include_path}"
    "${PROJECT_SOURCE_DIR}/antenna"
    "${PROJECT_SOURCE_DIR}/background"
    "${PROJECT_SOURCE_DIR}/core"
    "${PROJECT_SOURCE_DIR}/eigmode"
    "${PROJECT_SOURCE_DIR}/flre"
    "${PROJECT_SOURCE_DIR}/flre/conductivity"
    "${PROJECT_SOURCE_DIR}/flre/dispersion"
    "${PROJECT_SOURCE_DIR}/flre/maxwell_eqs"
    "${PROJECT_SOURCE_DIR}/flre/quants"
    "${PROJECT_SOURCE_DIR}/hom_medium"
    "${PROJECT_SOURCE_DIR}/imhd"
    "${PROJECT_SOURCE_DIR}/interface"
    "${PROJECT_SOURCE_DIR}/interp"
    "${PROJECT_SOURCE_DIR}/io"
    "${PROJECT_SOURCE_DIR}/math"
    "${PROJECT_SOURCE_DIR}/math/adapt_grid"
    "${PROJECT_SOURCE_DIR}/math/zersol-0.0.0/src"
    "${PROJECT_SOURCE_DIR}/mode"
    "${PROJECT_SOURCE_DIR}/solver"
    "${PROJECT_SOURCE_DIR}/spline"
)

set(cpp_sources
    antenna/antenna.cpp
    background/back_sett.cpp
    background/background.cpp
    background/calc_back.cpp
    background/eval_back.cpp
    core/core.cpp
    core/output_sett.cpp
    core/settings.cpp
    core/shared.cpp
    eigmode/calc_eigmode.cpp
    eigmode/eigmode_sett.cpp
    eigmode/find_eigmodes.cpp
    flre/conductivity/calc_cond.cpp
    flre/conductivity/cond_profs${interp}.cpp
    flre/conductivity/eval_cond.cpp
    flre/dispersion/disp_profs.cpp
    flre/flre_zone.cpp
    flre/maxwell_eqs/eval_sysmat.cpp
    flre/maxwell_eqs/maxwell_eqs_data.cpp
    flre/maxwell_eqs/sysmat_profs.cpp
    flre/quants/calc_flre_quants.cpp
    flre/quants/flre_quants.cpp
    hom_medium/hmedium_zone.cpp
    imhd/compressible_flow.cpp
    imhd/imhd_zone.cpp
    imhd/incompressible.cpp
    interface/wave_code_interface.cpp
    interp/interp.cpp
    math/adapt_grid/adaptive_grid.cpp
    math/adapt_grid/adaptive_grid_pol.cpp
    math/hyper/hyper1F1.cpp
    mode/calc_mode.cpp
    mode/mode.cpp
    mode/transforms.cpp
    mode/zone.cpp
    solver/rhs_func.cpp
    solver/solver.cpp
    spline/spline.cpp
)

set(fortran_sources
    antenna/antenna_m.f90
    antenna/antenna_spectrum.f90
    background/background_m.f90
    background/density_p${Jac}.f90
    background/f0moments${Jac}.f90
    core/constants_m.f90
    core/core_m.f90
    flre/conductivity/bisection.f90
    flre/conductivity/calc_D_array${drift}.f90
    flre/conductivity/calc_D_array_1_fac${drift}.f90
    flre/conductivity/calc_D_array_3_fac${drift}.f90
    flre/conductivity/calc_D_array_5_fac${drift}.f90
    flre/conductivity/calc_I_array${drift}.f90
    flre/conductivity/calc_W2_array${coll}${drift}.f90
    flre/conductivity/conduct_arrays${drift}.f90
    flre/conductivity/conduct_parameters.f90
    flre/conductivity/conductivity.f90
    flre/conductivity/delC.f90
    flre/conductivity/gaussfilter.f90
    flre/conductivity/gcorr.f90
    flre/conductivity/kmatrices${drift}.f90
    flre/conductivity/newton.f90
    flre/dispersion/dispersion.f90
    flre/dispersion/dpsort.f
    flre/maxwell_eqs/diff_sys.f90
    flre/maxwell_eqs/eqs_matrix_9.f90
    flre/maxwell_eqs/flre.f90
    flre/maxwell_eqs/maxwell_eqs_m.f90
    flre/maxwell_eqs/start_vals.f90
    flre/flre_sett_m.f90
    hom_medium/hmedium.f90
    imhd/imhd.f90
    interface/dummy.f90
    math/W/wfunc.f90
    math/W/wm_asympt.f90
    math/W/wmfunc.f90
    math/W/wofz.f
    math/bessel/bessel.f90
    math/det/calc_det.f90
    math/det/solve_hom_syst.f90
    math/eigsys/eigsys.f90
    math/odeint_allroutines.f
    math/solsubsys.f90
    math/transform.f90
    math/zgeev_tuned.f
    mode/mode_m.f90
    mode/stitching.f90
)

add_subdirectory(io)

add_library(${PROJECT_NAME} STATIC
    ${fortran_sources}
    ${cpp_sources}
    $<TARGET_OBJECTS:kilca_io>
)
target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)
set_target_properties(${PROJECT_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY lib
    Fortran_MODULE_DIRECTORY mod)
target_link_libraries(${PROJECT_NAME} PUBLIC ${EXTERNAL_LIBS})

add_executable(kilca_normal_exe "progs/main_linear.cpp")
set_target_properties(kilca_normal_exe PROPERTIES
    OUTPUT_NAME "kilca_normal"
    RUNTIME_OUTPUT_DIRECTORY exe)
target_link_libraries(kilca_normal_exe ${PROJECT_NAME} ${EXTERNAL_LIBS})

add_executable(kilca_eigparam_exe "progs/main_eig_param.cpp")
set_target_properties(kilca_eigparam_exe PROPERTIES
    OUTPUT_NAME "kilca_eig_param"
    RUNTIME_OUTPUT_DIRECTORY exe)
target_link_libraries(kilca_eigparam_exe ${PROJECT_NAME} ${EXTERNAL_LIBS})

add_executable(post_proc post_processing/main_post_proc.cpp post_processing/save_fortran.f90)
set_target_properties(post_proc PROPERTIES
    OUTPUT_NAME "post_processor"
    RUNTIME_OUTPUT_DIRECTORY exe)
target_link_libraries(post_proc ${PROJECT_NAME} ${EXTERNAL_LIBS})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
else()
  message(FATAL_ERROR "No doxygen tool is found in the system!")
endif()

enable_testing()
add_subdirectory(tests)
