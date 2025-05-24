set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    #set(CMAKE_C_COMPILER "gcc-13")
    #set(CMAKE_CXX_COMPILER "g++-13")
    #set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    #set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp -O3 -DDOUBLE_APPEND_FORTRAN")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -g -Og -Wall \
            -Wno-unused-function -fbacktrace")
    set(CMAKE_EXE_LINKER_FLAGS "-L/opt/homebrew/lib")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    # Linux configurations: use superproject’s Fortran compiler
    message("Building on Linux; using default Fortran compiler: ${CMAKE_Fortran_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system might not be implemented.")
endif()
enable_language(C Fortran)
find_package(MPI REQUIRED)

set(hdf5_tools ${CMAKE_CURRENT_SOURCE_DIR}/src/hdf5_tools/)
set(objs ${PROJECT_BINARY_DIR}/OBJS/)
set(lib ${PROJECT_BINARY_DIR}/lib/)
set(ql-balance_lib "${lib}libql-balance.a")

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

# Inlined CMakeSources.in contents:
set(src ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(base_balance_src ${src}/base/)
set(balance_lib_base_fortran_sources
    ${base_balance_src}/QLBalance_kinds.f90
    ${base_balance_src}/ql-balance_mod.f90
    ${base_balance_src}/diagnostic_mod.f90
    ${base_balance_src}/balanceBase.f90
    ${base_balance_src}/grid_mod.f90
    ${base_balance_src}/smooth_array_gauss.f90
    ${base_balance_src}/getIfunc.f90
    ${base_balance_src}/get_dql.f90
    ${base_balance_src}/diff_coeffs.f90
    ${base_balance_src}/control_mod.f90
    ${base_balance_src}/resonances_mod.f90
    ${base_balance_src}/matrix_mod.f90
    ${base_balance_src}/recstep_mod.f90
    ${base_balance_src}/baseparam_mod.f90
    ${base_balance_src}/h5mod.f90
    ${base_balance_src}/rhs_balance.f90
    ${base_balance_src}/calc_current_densities.f90
    ${base_balance_src}/wave_code_data_mod.f90
    ${base_balance_src}/Wmn.f90
    ${base_balance_src}/balance_eqs_source_terms.f90
    ${base_balance_src}/evolvestep.f90
    ${base_balance_src}/transp_coeffs_mod.f90
    ${base_balance_src}/read_config.f90
    ${base_balance_src}/gengrid.f90
    ${base_balance_src}/plag_coeff.f90
    ${base_balance_src}/amn_of_r.f90
    ${src}/hdf5_tools/hdf5_tools.f90
    ${base_balance_src}/paramscan.f90
    ${base_balance_src}/parallelTools.f90
    ${base_balance_src}/SingleStep.f90
    ${base_balance_src}/writeData.f90
    ${base_balance_src}/ramp_coil.f90
    ${base_balance_src}/time_evolution.f90
    ${base_balance_src}/restart.f90
    ${base_balance_src}/plasma_parameters.f90
    ${base_balance_src}/W2_arr.f90
    ${base_balance_src}/magnetic_island_width.f90
    ${base_balance_src}/integrate_parallel_current.f90
    ${base_balance_src}/wave_code_data_64bit.f90
    ${base_balance_src}/sample_matrix.f90
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

add_dependencies(ql-balance_lib kilca_lib)
target_link_libraries(ql-balance_lib PUBLIC kilca_lib MPI::MPI_Fortran)
target_include_directories(ql-balance_lib PUBLIC ${CMAKE_BINARY_DIR}/KiLCA/mod)

# Include central slatec math sources via cache variable defined by common/math
target_include_directories(ql-balance_lib PUBLIC ${SLATEC_INCLUDE_DIR})

# Link common dependencies defined in top-level Dependencies.cmake
target_link_libraries(ql-balance_lib PUBLIC slatec cvode nvecserial gsl gslcblas lapack blas)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(" ====== HDF5 was found! ======")
target_link_libraries(ql-balance_lib PUBLIC hdf5::hdf5
                                            hdf5::hdf5_hl
                                            hdf5::hdf5_fortran
                                            hdf5::hdf5_hl_fortran)

find_or_fetch(libneo)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/libneo/include)

# Dependencies (e.g., NetCDF) provided by top-level cmake/Dependencies.cmake
target_link_libraries(ql-balance_lib PUBLIC slatec 
                                            cvode 
                                            nvecserial 
                                            sparse
                                            netcdf
                                            netcdff
                                            LIBNEO::neo
                                            LIBNEO::species
                                            LIBNEO::collision_freqs
                                            LIBNEO::transport
                                            )

add_executable(ql-balance.x ${base_balance_src}/ql-balance.f90)
target_link_libraries(ql-balance.x PUBLIC ql-balance_lib)

# Add KiLCA module directory for Fortran modules
include_directories(${CMAKE_BINARY_DIR}/KiLCA/mod)
