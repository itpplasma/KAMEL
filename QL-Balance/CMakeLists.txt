set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    find_package(OpenMP REQUIRED)

    # Locate libomp from Homebrew
    find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
    if(BREW_EXECUTABLE)
        execute_process(
            COMMAND "${BREW_EXECUTABLE}" --prefix libomp
            OUTPUT_VARIABLE LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(EXISTS "${LIBOMP_PREFIX}")
            include_directories(SYSTEM "${LIBOMP_PREFIX}/include")
            link_directories("${LIBOMP_PREFIX}/lib")
        endif()
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -DDOUBLE_APPEND_FORTRAN")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -g -Og -Wall -Wno-unused-function -fbacktrace")
    set(CMAKE_EXE_LINKER_FLAGS "-L/opt/homebrew/lib -L${LIBOMP_PREFIX}/lib -lomp")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    # Linux configurations: use superproject's Fortran compiler
    message("Building on Linux; using default Fortran compiler: ${CMAKE_Fortran_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system might not be implemented.")
endif()
enable_language(C Fortran)
set(objs ${PROJECT_BINARY_DIR}/OBJS/ql-balance/)
set(lib ${PROJECT_BINARY_DIR}/install/lib/)
set(ql-balance_lib "${lib}libql-balance.a")

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

# Inlined CMakeSources.in contents:
set(src ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(base_balance_src ${src}/base/)
set(balance_lib_base_fortran_sources
    ${base_balance_src}/amn_of_r.f90
    ${base_balance_src}/balance_eqs_source_terms.f90
    ${base_balance_src}/balanceBase.f90
    ${base_balance_src}/baseparam_mod.f90
    ${base_balance_src}/calc_current_densities.f90
    ${base_balance_src}/control_mod.f90
    ${base_balance_src}/diagnostic_mod.f90
    ${base_balance_src}/diff_coeffs.f90
    ${base_balance_src}/evolvestep.f90
    ${base_balance_src}/gengrid.f90
    ${base_balance_src}/get_dql.f90
    ${base_balance_src}/getIfunc.f90
    ${base_balance_src}/grid_mod.f90
    ${base_balance_src}/h5mod.f90
    ${base_balance_src}/integrate_parallel_current.f90
    ${base_balance_src}/magnetic_island_width.f90
    ${base_balance_src}/matrix_mod.f90
    ${base_balance_src}/paramscan.f90
    ${base_balance_src}/plag_coeff.f90
    ${base_balance_src}/plasma_parameters.f90
    ${base_balance_src}/ql-balance_mod.f90
    ${base_balance_src}/QLBalance_kinds.f90
    ${base_balance_src}/ramp_coil.f90
    ${base_balance_src}/read_config.f90
    ${base_balance_src}/recstep_mod.f90
    ${base_balance_src}/resonances_mod.f90
    ${base_balance_src}/restart.f90
    ${base_balance_src}/rhs_balance_m.f90
    ${base_balance_src}/sample_matrix.f90
    ${base_balance_src}/SingleStep.f90
    ${base_balance_src}/smooth_array_gauss.f90
    ${base_balance_src}/time_evolution.f90
    ${base_balance_src}/toroidal_torque.f90
    ${base_balance_src}/transp_coeffs_mod.f90
    ${base_balance_src}/W2_arr.f90
    ${base_balance_src}/wave_code_data_64bit.f90
    ${base_balance_src}/wave_code_data_mod.f90
    ${base_balance_src}/Wmn.f90
    ${base_balance_src}/writeData.f90
)
set(balance_lib_stellarator_fortran_sources
    ${src}/stellarator/D_one_over_nu.f90
    ${src}/stellarator/rhs_balance_stell.f90
    ${src}/stellarator/balance_eqs_source_terms_stell.f90
    ${src}/stellarator/time_evol_stell.f90
    ${src}/stellarator/evolvestep_stell.f90
)
set(balance_lib_cpp_sources
    ${base_balance_src}/vel_integral.cpp
    ${base_balance_src}/cvodeint.cpp
)

add_library(ql-balance_lib STATIC ${balance_lib_base_fortran_sources}
                                    ${balance_lib_stellarator_fortran_sources}
                                    ${balance_lib_cpp_sources})

set_target_properties(ql-balance_lib PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME ql-balance
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")

add_dependencies(ql-balance_lib kilca_lib kamel_hdf5_tools)
target_include_directories(ql-balance_lib PUBLIC ${CMAKE_BINARY_DIR}/KiLCA/mod ${CMAKE_BINARY_DIR}/modules)

# Include central slatec math sources via cache variable defined by common/math
target_include_directories(ql-balance_lib PUBLIC ${SLATEC_INCLUDE_DIR})

find_package(HDF5 REQUIRED COMPONENTS C HL Fortran Fortran_HL)
message(" ====== HDF5 was found! ======")

set(BUILD_TESTING OFF)

set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "Disable libneo tests" FORCE)
set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "Disable libneo tests" FORCE)


# Only fetch libneo if not already available (prevents duplicate targets)
if(NOT TARGET LIBNEO::neo)
    find_or_fetch(libneo)
endif()
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../libneo/include)

add_dependencies(ql-balance_lib sparse)

# Link common dependencies defined in top-level Dependencies.cmake
target_link_libraries(ql-balance_lib PUBLIC
    LIBNEO::collision_freqs
    LIBNEO::neo
    LIBNEO::species
    LIBNEO::transport
    SUNDIALS::cvode
    SUNDIALS::nvecserial
    blas
    gsl
    gslcblas
    hdf5::hdf5
    hdf5::hdf5_fortran
    hdf5::hdf5_hl
    hdf5::hdf5_hl_fortran
    kamel_hdf5_tools
    kilca_lib
    lapack
    netcdf
    netcdff
    slatec
    sparse
)

target_include_directories(sparse PUBLIC "${CMAKE_BINARY_DIR}/OBJS/sparse/")

add_executable(ql-balance.x ${base_balance_src}/ql-balance.f90)
target_link_libraries(ql-balance.x PUBLIC ql-balance_lib)
set_target_properties(ql-balance.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/install/bin)

# Add KiLCA module directory for Fortran modules
include_directories(${CMAKE_BINARY_DIR}/KiLCA/mod)

add_subdirectory("${src}/test")
