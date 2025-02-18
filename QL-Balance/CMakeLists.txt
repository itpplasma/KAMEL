cmake_minimum_required(VERSION 3.24)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(ql-balance)
enable_language(C Fortran)
enable_testing()

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    #set(CMAKE_C_COMPILER /opt/homebrew/bin/gfortran-12)
    #set(CMAKE_CXX_COMPILER /opt/homebrew/bin/gfortran-12)
    #set(CMAKE_CXX_COMPILER /Library/Developer/CommandLineTools/usr/bin/c++)
    #set(CMAKE_C_COMPILER /Library/Developer/CommandLineTools/usr/bin/cc)

    set(CMAKE_C_COMPILER /opt/homebrew/opt/llvm/bin/clang)
    set(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp -O3 -DDOUBLE_APPEND_FORTRAN")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    # Linux specific configurations or commands
    message("Building on Linux")
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_CXX_COMPILER "g++")
    set(CMAKE_Fortran_COMPILER "mpif90")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system might not be implemented.")
endif()

include(FetchSuiteSparse)
include(Util)

set (libs_root_path ${CMAKE_SOURCE_DIR}/../external)

set(src ${PROJECT_SOURCE_DIR}/src)
set(full_balance_src ${src}/full/)

set(hdf5_tools ${PROJECT_SOURCE_DIR}/src/hdf5_tools/)
set(objs ${PROJECT_BINARY_DIR}/OBJS/)
set(lib ${PROJECT_BINARY_DIR}/lib/)
set(ql-balance_lib "${lib}libql-balance.a")
set(KiLCA_lib_path ${CMAKE_SOURCE_DIR}/../KiLCA/build/lib/libKiLCA_Lib_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit.a)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

add_subdirectory(${src})
add_subdirectory(${src}/tests)

set (sundials_include_path ${libs_root_path}/sundials/include)
set (sundials_include_path_build ${libs_root_path}/sundials/build/include)
include_directories ("${sundials_include_path}")
include_directories ("${sundials_include_path_build}")

set(balance_lib_fortran_sources
    ${full_balance_src}/ql-balance_mod.f90
    ${full_balance_src}/diagnostic_mod.f90
    ${full_balance_src}/balanceBase.f90
    ${full_balance_src}/grid_mod.f90
    ${full_balance_src}/control_mod.f90
    ${full_balance_src}/resonances_mod.f90
    ${full_balance_src}/matrix_mod.f90
    ${full_balance_src}/recstep_mod.f90
    ${full_balance_src}/baseparam_mod.f90
    ${full_balance_src}/h5mod.f90
    ${full_balance_src}/rhs_balance.f90
    ${full_balance_src}/wave_code_data_mod.f90
    ${full_balance_src}/Wmn.f90
    ${full_balance_src}/evolvestep.f90 
    ${full_balance_src}/transp_coeffs_mod.f90 
    ${full_balance_src}/read_config.f90 
    ${full_balance_src}/gengrid.f90
    ${full_balance_src}/plag_coeff.f90
    ${full_balance_src}/amn_of_r.f90
    ${src}/hdf5_tools/hdf5_tools.f90
    ${full_balance_src}/paramscan.f90
    ${full_balance_src}/parallelTools.f90
    ${full_balance_src}/SingleStep.f90
    ${full_balance_src}/writeData.f90
    ${full_balance_src}/time_evolution.f90
    ${full_balance_src}/restart.f90
    ${full_balance_src}/plasma_parameters.f90
    ${full_balance_src}/W2_arr.f90
    ${full_balance_src}/magnetic_island_width.f90
    ${full_balance_src}/wave_code_data_64bit.f90
    ${full_balance_src}/sample_matrix.f90)

set(balance_lib_cpp_sources
    ${full_balance_src}/vel_integral.cpp
    ${full_balance_src}/cvodeint.cpp
    )

add_library(sparse STATIC
    ${full_balance_src}/sparse_mod.f90
	${UMFPACK_WRAPPERS}
)
set_target_properties(sparse PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME sparse
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")

add_library(ql-balance_lib STATIC ${balance_lib_fortran_sources} ${balance_lib_cpp_sources})

set_target_properties(ql-balance_lib PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME ql-balance
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")


target_link_libraries(ql-balance_lib PUBLIC "${KiLCA_lib_path}")
target_include_directories(ql-balance_lib PUBLIC "${sundials_include_path}")
target_include_directories(ql-balance_lib PUBLIC "${sundials_include_path_build}")
target_include_directories(ql-balance_lib PUBLIC ${CMAKE_SOURCE_DIR}/../KiLCA/math/)

file(GLOB libslatec_files "${CMAKE_SOURCE_DIR}/../KiLCA/math/bessel/slatec/*.f")
add_library (slatec STATIC "${libslatec_files}")
set_target_properties(slatec PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME slatec
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")

set (cvode_lib_fullpath ${libs_root_path}/sundials/build/src/cvode/libsundials_cvode.a)
add_library (cvode STATIC IMPORTED)
set_target_properties (cvode PROPERTIES IMPORTED_LOCATION "${cvode_lib_fullpath}")

set (nvecserial_lib_fullpath ${libs_root_path}/sundials/build/src/nvector/serial/libsundials_nvecserial.a)
add_library (nvecserial STATIC IMPORTED)
set_target_properties (nvecserial PROPERTIES IMPORTED_LOCATION "${nvecserial_lib_fullpath}")

find_package(GSL REQUIRED)
message(" ====== GSL was found! ======")
target_link_libraries(ql-balance_lib PUBLIC GSL::gsl)

find_package(MPI REQUIRED)
message(" ====== MPI was found! ======")
target_link_libraries(ql-balance_lib PUBLIC MPI::MPI_Fortran)

add_dependencies(sparse SuiteSparse)
target_link_libraries(sparse PUBLIC SuiteSparse::suitesparse_config)
target_link_libraries(sparse PUBLIC SuiteSparse::amd)
target_link_libraries(sparse PUBLIC SuiteSparse::umfpack)
target_include_directories(sparse PUBLIC "${UMFPACK_INCLUDE_DIR}")

find_package(BLAS REQUIRED)
message(" ====== BLAS was found! ======")
target_link_libraries(ql-balance_lib PUBLIC BLAS::BLAS)
target_link_libraries(sparse PUBLIC BLAS::BLAS)

find_package(LAPACK REQUIRED)
message(" ====== LAPACK was found! ======")
target_link_libraries(ql-balance_lib PUBLIC LAPACK::LAPACK)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(" ====== HDF5 was found! ======")
target_link_libraries(ql-balance_lib PUBLIC hdf5::hdf5
                                            hdf5::hdf5_hl
                                            hdf5::hdf5_fortran
                                            hdf5::hdf5_hl_fortran)


#find_or_fetch(libneo)
#include_directories(${CMAKE_CURRENT_BINARY_DIR}/libneo/include)

include(FetchNetcdf)
#target_link_libraries(ql-balance_lib PUBLIC ${NETCDFLIB_DIR}/libnetcdf.a)
#target_link_libraries(ql-balance_lib PUBLIC ${NETCDF_FLIBS})
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_FLIBS}")
#target_include_directories(ql-balance_lib PUBLIC ${NETCDFINCLUDE_DIR})

include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})

# Extract the path (anything that starts with -L)
string(REGEX MATCHALL "-L[^ ]+" NETCDF_F_PATH_ONLY "${NETCDF_FLIBS}")

# Remove the "-L" prefix to get only the raw path
string(REPLACE "-L" "" NETCDF_F_PATH "${NETCDF_F_PATH_ONLY}")

add_link_options(${NETCDF_F_PATH_ONLY})

target_link_libraries(ql-balance_lib PUBLIC slatec 
                                            cvode 
                                            nvecserial 
                                            sparse)
                                            #LIBNEO::neo
                                            #LIBNEO::collision_freqs)

add_executable(ql-balance.x ${full_balance_src}/ql-balance.f90)
target_link_libraries(ql-balance.x PUBLIC ql-balance_lib)


