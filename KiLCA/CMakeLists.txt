if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    #set(CMAKE_C_COMPILER gcc-14)
    set(CMAKE_C_COMPILER clang)
    #set(CMAKE_CXX_COMPILER g++-14)
    set(CMAKE_CXX_COMPILER clang++)
    # flag does not work for macos for some reason
    set (CMAKE_CXX_OPT_FLAGS "-O3")
    set (CMAKE_F_OPT_FLAGS "-O3")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux" AND
        CMAKE_SYSTEM_PROCESSOR MATCHES "(x86|i386|i686|amd64|x86_64)")
    # Linux specific configurations or commands
    message("Building on Linux")
    set (CMAKE_C_COMPILER   gcc)
    set (CMAKE_CXX_COMPILER g++)
    set (CMAKE_EXE_LINKER_FLAGS "")
    set (CMAKE_CXX_OPT_FLAGS "-O3 -march=x86-64-v2 -mtune=generic")
    set (CMAKE_F_OPT_FLAGS "-O3 -march=x86-64-v2 -mtune=generic")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system. Might not be implemented.")
endif()

set (KiLCA_ARCHITECTURE 64bit) #possible values: 32bit, 64bit; must conform to compilers and libs versions in use

set (KiLCA_VERSION_MAJOR 2)
set (KiLCA_VERSION_MINOR 4.2)

set (MAGNETIC_DRIFT MDNO) #possible values: MDNO, MDYES

set (COLLISION_MODEL FPGEN) #possible values: FPN, FPNE, FPNEP, FPGEN - with collisions, NC - no collisions

set (JACOBIAN SIMPLIFIED) #possible values: SIMPLIFIED - simplifyed jacobian (preffered), FULL

set (ADAPTIVE_GRID_GENERATOR POLYNOMIAL) #possible values: POLYNOMIAL - polynomial, USUAL -old
#set (ADAPTIVE_GRID_GENERATOR USUAL) #possible values: POLYNOMIAL - polynomial, USUAL

set (ADDITIONAL_COMMENT _) #possible values: anything

set (DEBUG_FLAG 0) #possible values: 0 - no debugging, 1 - with debugging

set (CALC_EIGEN_DECOMPOSITION 0) #possible values: 0 - if no, 1 - if yes

set (SORT_DISPERSION_PROFILES 0) #possible values: 0 - if no, 1 - if yes

set (USE_JACOBIAN_IN_ODE_SOLVER 1) #possible values: 0 - if no, 1 - if yes

set (USE_SPLINES_IN_RHS_EVALUATION 0) #possible values: 0 - if no, 1 - if yes

set (kilca_sign V_${KiLCA_VERSION_MAJOR}.${KiLCA_VERSION_MINOR}_${MAGNETIC_DRIFT}_${COLLISION_MODEL}_${ADAPTIVE_GRID_GENERATOR}_${CMAKE_BUILD_TYPE}${ADDITIONAL_COMMENT}${KiLCA_ARCHITECTURE})

set (kilca_lib KiLCA_Lib_${kilca_sign})
set (kilca_normal_exe KiLCA_Normal_${kilca_sign})
set (kilca_eigparam_exe KiLCA_EigParam_${kilca_sign})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
include_directories("${CMAKE_Fortran_MODULE_DIRECTORY}")

#Due to gcc compiler bug the code may produce wrong results if it is compiled with optimizations -O2 flags!!!

if (CMAKE_BUILD_TYPE STREQUAL Release)
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=generic ${CMAKE_CXX_OPT_FLAGS}")
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2 -mtune=generic ${CMAKE_F_OPT_FLAGS}")
elseif (CMAKE_BUILD_TYPE STREQUAL Debug)
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wshadow -mtune=generic")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O2 -Wall -Wextra -Wshadow -Wtabs -mtune=generic")
endif ()

if (KiLCA_ARCHITECTURE STREQUAL 32bit)
    set (arch _32bit)
elseif (KiLCA_ARCHITECTURE STREQUAL 64bit)
    set (arch _64bit)
else (KiLCA_ARCHITECTURE STREQUAL 32bit)
    message (FATAL_ERROR "Unknown KiLCA architecture: ${KiLCA_ARCHITECTURE}")
endif ()

if (MAGNETIC_DRIFT STREQUAL MDNO)
    set (drift )
else ()
    set (drift _drift)
endif ()

if (COLLISION_MODEL STREQUAL FPN)
    set (coll _wc_n)
elseif (COLLISION_MODEL STREQUAL FPNE)
    set (coll _wc_ne)
elseif (COLLISION_MODEL STREQUAL FPNEP)
    set (coll _wc_nep)
elseif (COLLISION_MODEL STREQUAL NC)
    set (coll _nc)
elseif (COLLISION_MODEL STREQUAL FPGEN)
    set (coll _gen)
else ()
    message (FATAL_ERROR "Unknown collision model: ${COLLISION_MODEL}")
endif ()

if (ADAPTIVE_GRID_GENERATOR STREQUAL POLYNOMIAL)
    set (interp _pol)
else ()
    set (interp )
endif ()

if (JACOBIAN STREQUAL SIMPLIFIED)
    set (Jac _J0)
else ()
    set (Jac )
endif ()

set (kilca_normal_main progs/main_linear.cpp)
set (kilca_eigparam_main progs/main_eig_param.cpp)

# --- Dependencies ---
# Dependencies are provided by the top-level project via cmake/Dependencies.cmake

# include paths:
set (libs_root_path ${CMAKE_BINARY_DIR}/external-install)
set (sundials_include_path ${libs_root_path}/sundials/include)
set (gsl_include_path ${libs_root_path}/gsl/include)

set(EXTERNAL_LIBS slatec lapack blas cvode nvecserial gsl gslcblas)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/code_settings.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/code_settings.h"
)
# Make generated header available
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

include_directories ("${sundials_include_path}")
include_directories ("${gsl_include_path}")
include_directories ("${PROJECT_BINARY_DIR}")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/antenna")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/background")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/core")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/eigmode")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/flre")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/flre/conductivity")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/flre/dispersion")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/flre/maxwell_eqs")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/flre/quants")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/hom_medium")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/imhd")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/interface")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/interp")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/io")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/math")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/math/adapt_grid")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/math/zersol-0.0.0/src")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/mod")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/mode")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/progs")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/solver")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/spline")

#cpp source files for kilca lib:
set(kilca_lib_cpp_sources
    core/shared.cpp
    core/core.cpp
    core/settings.cpp
    core/output_sett.cpp
    background/back_sett.cpp
    background/background.cpp
    background/calc_back.cpp
    background/eval_back.cpp
    eigmode/eigmode_sett.cpp
    eigmode/calc_eigmode.cpp
    eigmode/find_eigmodes.cpp
    antenna/antenna.cpp
    flre/flre_zone.cpp
    flre/conductivity/calc_cond.cpp
    flre/conductivity/cond_profs${interp}.cpp
    flre/conductivity/eval_cond.cpp
    flre/dispersion/disp_profs.cpp
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
    io/inout.cpp
    mode/mode.cpp
    mode/calc_mode.cpp
    mode/zone.cpp
    mode/transforms.cpp
    solver/solver.cpp
    solver/rhs_func.cpp
    spline/spline.cpp
    math/hyper/hyper1F1.cpp
    math/adapt_grid/adaptive_grid.cpp
    math/adapt_grid/adaptive_grid_pol.cpp
)

#fortran source files for kilca lib:
set(kilca_lib_fortran_sources
    core/constants_m${arch}.f90
    core/core_m.f90
    antenna/antenna_m.f90
    background/background_m.f90
    flre/flre_sett_m.f90
    mode/mode_m.f90
    flre/maxwell_eqs/maxwell_eqs_m.f90
    background/f0moments${Jac}.f90
    background/density_p${Jac}.f90
    antenna/antenna_spectrum.f90
    flre/conductivity/calc_D_array_1_fac${drift}.f90
    flre/conductivity/calc_D_array_3_fac${drift}.f90
    flre/conductivity/calc_D_array_5_fac${drift}.f90
    flre/conductivity/calc_D_array${drift}.f90
    #flre/conductivity/calc_I_array${drift}_serg.f90
    flre/conductivity/calc_I_array${drift}.f90
    flre/conductivity/calc_W2_array${coll}${drift}.f90
    flre/conductivity/conduct_arrays${drift}.f90
    flre/conductivity/conduct_parameters.f90
    flre/conductivity/conductivity.f90
    flre/conductivity/delC.f90
    flre/conductivity/gcorr.f90
    flre/conductivity/kmatrices${drift}.f90
    flre/conductivity/newton.f90
    flre/conductivity/bisection.f90
    flre/conductivity/gaussfilter.f90
    flre/dispersion/dispersion.f90
    flre/dispersion/dpsort.f
    flre/maxwell_eqs/eqs_matrix_9.f90
    flre/maxwell_eqs/diff_sys.f90
    flre/maxwell_eqs/start_vals.f90
    flre/maxwell_eqs/flre.f90
    hom_medium/hmedium.f90
    imhd/imhd.f90
    interface/dummy.f90
    mode/stitching.f90
    math/eigsys/eigsys.f90
    math/bessel/bessel.f90
    math/solsubsys.f90
    math/transform.f90
    math/zgeev_tuned.f
    math/W/wofz.f
    math/W/wfunc.f90
    math/W/wmfunc.f90
    math/W/wm_asympt.f90
    math/det/calc_det.f90
    math/det/solve_hom_syst.f90
)

add_library (kilca_lib STATIC ${kilca_lib_fortran_sources} ${kilca_lib_cpp_sources})
#add_library (kilca_lib SHARED ${kilca_lib_fortran_sources} ${kilca_lib_cpp_sources})
target_compile_options(kilca_lib PRIVATE -fPIC)

set_target_properties(kilca_lib PROPERTIES
    OUTPUT_NAME ${kilca_lib}
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/install/lib/
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS/kilca/)
add_dependencies(kilca_lib ${EXTERNAL_LIBS})
target_link_libraries(kilca_lib PUBLIC ${EXTERNAL_LIBS})

# Add include directories on `kilca_lib` target for all source and module paths
target_include_directories(kilca_lib PUBLIC
    # Fortran modules generated in build
    ${CMAKE_CURRENT_BINARY_DIR}/mod
    # Generated headers
    ${CMAKE_CURRENT_BINARY_DIR}
    # KiLCA source subdirectories
    ${CMAKE_CURRENT_SOURCE_DIR}/antenna
    ${CMAKE_CURRENT_SOURCE_DIR}/background
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/eigmode
    ${CMAKE_CURRENT_SOURCE_DIR}/flre
    ${CMAKE_CURRENT_SOURCE_DIR}/flre/conductivity
    ${CMAKE_CURRENT_SOURCE_DIR}/flre/dispersion
    ${CMAKE_CURRENT_SOURCE_DIR}/flre/maxwell_eqs
    ${CMAKE_CURRENT_SOURCE_DIR}/flre/quants
    ${CMAKE_CURRENT_SOURCE_DIR}/hom_medium
    ${CMAKE_CURRENT_SOURCE_DIR}/imhd
    ${CMAKE_CURRENT_SOURCE_DIR}/interface
    ${CMAKE_CURRENT_SOURCE_DIR}/interp
    ${CMAKE_CURRENT_SOURCE_DIR}/io
    ${CMAKE_CURRENT_SOURCE_DIR}/math
    ${CMAKE_CURRENT_SOURCE_DIR}/math/adapt_grid
    ${CMAKE_CURRENT_SOURCE_DIR}/math/zersol-0.0.0/src
    ${CMAKE_CURRENT_SOURCE_DIR}/mod
    ${CMAKE_CURRENT_SOURCE_DIR}/mode
    ${CMAKE_CURRENT_SOURCE_DIR}/progs
    ${CMAKE_CURRENT_SOURCE_DIR}/solver
    ${CMAKE_CURRENT_SOURCE_DIR}/spline
    # External dependencies
    ${sundials_include_path}
    ${gsl_include_path}
)

add_executable (kilca_normal_exe ${kilca_normal_main})
set_target_properties(kilca_normal_exe PROPERTIES
    OUTPUT_NAME ${kilca_normal_exe}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/install/bin/)
add_dependencies(kilca_normal_exe ${EXTERNAL_LIBS})
target_link_libraries (kilca_normal_exe kilca_lib ${EXTERNAL_LIBS})

add_executable (kilca_eigparam_exe ${kilca_eigparam_main})
set_target_properties(kilca_eigparam_exe PROPERTIES
    OUTPUT_NAME ${kilca_eigparam_exe}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/install/bin/)
add_dependencies(kilca_eigparam_exe ${EXTERNAL_LIBS})
target_link_libraries (kilca_eigparam_exe kilca_lib ${EXTERNAL_LIBS})

add_executable (post_proc post_processing/main_post_proc.cpp post_processing/save_fortran.f90)
set_target_properties(post_proc PROPERTIES
    OUTPUT_NAME post_processor_${KiLCA_ARCHITECTURE}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/install/bin/)
add_dependencies(post_proc ${EXTERNAL_LIBS})
target_link_libraries (post_proc kilca_lib ${EXTERNAL_LIBS})

# add a target to generate API documentation with Doxygen
find_package (Doxygen)

if (DOXYGEN_FOUND)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile @ONLY)
  add_custom_target (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
else ()
  message (FATAL_ERROR "No doxygen tool is found in the system!")
endif ()
