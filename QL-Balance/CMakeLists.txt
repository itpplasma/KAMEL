cmake_minimum_required(VERSION 3.24)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(ql-balance)
enable_testing()

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
    # Linux specific configurations or commands
    message("Building on Linux")
    set(CMAKE_Fortran_COMPILER /usr/bin/mpifort)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system might not be implemented.")
endif()
enable_language(C Fortran)

include(FetchSuiteSparse)
include(Util)

set (libs_root_path ${CMAKE_SOURCE_DIR}/../external)

set(hdf5_tools ${PROJECT_SOURCE_DIR}/src/hdf5_tools/)
set(objs ${PROJECT_BINARY_DIR}/OBJS/)
set(lib ${PROJECT_BINARY_DIR}/lib/)
set(ql-balance_lib "${lib}libql-balance.a")
set(KiLCA_lib_path ${CMAKE_SOURCE_DIR}/../KiLCA/build/lib/libKiLCA_Lib_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit.a)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

set (sundials_include_path ${libs_root_path}/sundials/include)
set (sundials_include_path_build ${libs_root_path}/sundials/build/include)
include_directories ("${sundials_include_path}")
include_directories ("${sundials_include_path_build}")

include(CMakeSources.in)

add_library(sparse STATIC
    ${base_balance_src}/sparse_mod.f90
)
set_target_properties(sparse PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME sparse
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")
add_dependencies(sparse SuiteSparse)
target_link_libraries(sparse PUBLIC SuiteSparse::umfpack_wrappers)
target_link_libraries(sparse PUBLIC SuiteSparse::umfpack)
target_include_directories(sparse PUBLIC "${UMFPACK_INCLUDE_DIR}")

add_library(ql-balance_lib STATIC ${balance_lib_base_fortran_sources} 
                                    ${balance_lib_stellarator_fortran_sources}
                                    ${balance_lib_cpp_sources})

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
target_link_libraries(ql-balance_lib PUBLIC MPI::MPI_C)
target_link_libraries(ql-balance_lib PUBLIC MPI::MPI_CXX)

find_package(BLAS REQUIRED)
message(" ====== BLAS was found! ======")
target_link_libraries(ql-balance_lib PUBLIC BLAS::BLAS)

find_package(LAPACK REQUIRED)
message(" ====== LAPACK was found! ======")
target_link_libraries(ql-balance_lib PUBLIC LAPACK::LAPACK)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
message(" ====== HDF5 was found! ======")
target_link_libraries(ql-balance_lib PUBLIC hdf5::hdf5
                                            hdf5::hdf5_hl
                                            hdf5::hdf5_fortran
                                            hdf5::hdf5_hl_fortran)

#set (CMAKE_Fortran_FLAGS "-cpp -g -Og -fPIC -fcheck=all -fbacktrace \
#        -fno-realloc-lhs -fopenmp -Wall -Wextra -Werror \
#        -Wno-unused-function -fmax-errors=1")

find_or_fetch(libneo)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/libneo/include)

include(FetchNetcdf)

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

#set_property(TARGET ql-balance_lib PROPERTY LINKER_LANGUAGE Fortran)

add_executable(ql-balance.x ${base_balance_src}/ql-balance.f90)
target_link_libraries(ql-balance.x PUBLIC ql-balance_lib)
