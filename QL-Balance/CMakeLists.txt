set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

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

# KI LCA path relative to subproject
set(KILCA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../KiLCA)

set(libs_root_path ${KILCA_PATH}/build/download/src) # libneo depends on this, cannot easily be changed

set(hdf5_tools ${CMAKE_CURRENT_SOURCE_DIR}/src/hdf5_tools/)
set(objs ${PROJECT_BINARY_DIR}/OBJS/)
set(lib ${PROJECT_BINARY_DIR}/lib/)
set(ql-balance_lib "${lib}libql-balance.a")

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

set (sundials_include_path ${libs_root_path}/sundials/include)
set (sundials_include_path_build ${libs_root_path}/sundials/build/include)
include_directories ("${sundials_include_path}")
include_directories ("${sundials_include_path_build}")

include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeSources.in")

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

target_include_directories(ql-balance_lib PUBLIC "${sundials_include_path}")
target_include_directories(ql-balance_lib PUBLIC "${sundials_include_path_build}")
# Include central slatec math sources
target_include_directories(ql-balance_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../common/math/slatec")

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
