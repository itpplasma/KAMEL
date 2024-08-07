cmake_minimum_required(VERSION 3.7...3.21)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file. Please make a build subdirectory. Feel free to remove CMakeCache and CMakeFiles.")
endif()

project(ql-balance)
enable_language(Fortran)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    set(CMAKE_C_COMPILER gcc-14)
    set(CMAKE_CXX_COMPILER g++-14)
    # flag does not work for macos for some reason
    set (CMAKE_CXX_OPT_FLAGS "")
    set (CMAKE_F_OPT_FLAGS "")
    set (CMAKE_CXX_OPT_FLAGS "")
    set (CMAKE_F_OPT_FLAGS "")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    # Linux specific configurations or commands
    message("Building on Linux")
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_CXX_COMPILER "g++")
    #set(CMAKE_Fortran_COMPILER "gfortran")
    set(CMAKE_Fortran_COMPILER "mpif90")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/lib64 -lstdc++")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -L/usr/lib64 -lstdc++ -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp -O3")
else()
    # Other operating system specific configurations or commands
    message("Building on another operating system. Might not be implemented.")
endif()


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set (libs_root_path ${CMAKE_SOURCE_DIR}/../../libs)

set(src ${PROJECT_SOURCE_DIR}/src)
set(full_balance_src ${src}/full/)
set(light_balance_src ${src}/light/)

set(hdf5_tools ${PROJECT_SOURCE_DIR}/hdf5_tools/)
set(objs ${PROJECT_SOURCE_DIR}/OBJS/)
set(lib ${PROJECT_SOURCE_DIR}/lib/)
set(ql-balance_lib "${lib}libql-balance.a")
set(KiLCA_lib_path ${CMAKE_SOURCE_DIR}/../KiLCA/build/lib/libKiLCA_Lib_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit.a)

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

add_subdirectory(${src})

set (sundials_include_path ${libs_root_path}/sundials/include)
set (sundials_include_path_build ${libs_root_path}/sundials/build/include)
include_directories ("${sundials_include_path}")
include_directories ("${sundials_include_path_build}")


set(balance_lib_fortran_sources
    ${full_balance_src}/ql-balance_mod.f90
    ${full_balance_src}/grid_mod.f90
    ${full_balance_src}/h5mod.f90
    ${full_balance_src}/rhs_balance.f90
    ${full_balance_src}/wave_code_data_mod.f90
    ${full_balance_src}/Wmn.f90
    ${full_balance_src}/evolvestep.f90 
    ${full_balance_src}/read_config.f90 
    ${full_balance_src}/gengrid.f90
    ${full_balance_src}/plag_coeff.f90
    ${full_balance_src}/amn_of_r.f90
    ${src}/hdf5_tools/hdf5_tools.f90
    ${full_balance_src}/sparse_mod.f90
    ${full_balance_src}/paramscan.f90
    ${full_balance_src}/parallelTools.f90
    ${full_balance_src}/linear_run.f90
    ${full_balance_src}/writeData.f90
    ${full_balance_src}/time_evolution.f90
    ${full_balance_src}/restart.f90
    ${full_balance_src}/plasma_parameters.f90
    ${full_balance_src}/W2_arr.f90
    ${full_balance_src}/magnetic_island_width.f90
    ${full_balance_src}/diagnostic_mod.f90
    ${full_balance_src}/wave_code_data_64bit.f90
    ${full_balance_src}/sample_matrix.f90)

set(balance_lib_cpp_sources
    ${full_balance_src}/vel_integral.cpp
    ${full_balance_src}/cvodeint.cpp
    ${full_balance_src}/umf4_f77wrapper.c
    ${full_balance_src}/umf4_f77zwrapper.c
    )

add_library(ql-balance_lib STATIC ${balance_lib_fortran_sources} ${balance_lib_cpp_sources})

set_target_properties(ql-balance_lib PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME ql-balance-lib
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")


target_link_libraries(ql-balance_lib PUBLIC "${KiLCA_lib_path}")

target_include_directories(ql-balance_lib PUBLIC "${sundials_include_path}")
target_include_directories(ql-balance_lib PUBLIC "${sundials_include_path_build}")
target_include_directories(ql-balance_lib PUBLIC ${CMAKE_SOURCE_DIR}/../KiLCA/math/)


file(GLOB libslatec_files "${CMAKE_SOURCE_DIR}/../KiLCA/math/bessel/slatec/*.f")
add_library (slatec STATIC "${libslatec_files}")

set(superlu_path ${libs_root_path}/SuperLU_4.1/lib/libsuperlu_4.1.a)
add_library(superlu STATIC IMPORTED)
set_target_properties(superlu PROPERTIES IMPORTED_LOCATION "${superlu_path}")

set (cvode_lib_fullpath ${libs_root_path}/sundials/build/src/cvode/libsundials_cvode.a)
add_library (cvode STATIC IMPORTED)
set_target_properties (cvode PROPERTIES IMPORTED_LOCATION "${cvode_lib_fullpath}")

set (nvecserial_lib_fullpath ${libs_root_path}/sundials/build/src/nvector/serial/libsundials_nvecserial.a)
add_library (nvecserial STATIC IMPORTED)
set_target_properties (nvecserial PROPERTIES IMPORTED_LOCATION "${nvecserial_lib_fullpath}")

find_package(GSL REQUIRED)
if("${GSL_FOUND}")
    message(" ====== GSL was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(ql-balance_lib PUBLIC "${GSL_CBLAS_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${GSL_INCLUDE_DIRS}")
endif()

find_package(MPI REQUIRED)
if("${MPI_FOUND}")
    message(" ====== MPI was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${MPI_Fortran_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${MPI_Fortran_INCLUDE_PATH}")
endif()

find_package(SuiteSparse REQUIRED)
if("${SUITESPARSE_FOUND}")
    message(" ====== SuiteSparse was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${SUITESPARSE_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${SUITESPARSE_INCLUDE_DIRS}")
endif()

find_package(BLAS REQUIRED)
if("${BLAS_FOUND}")
    message(" ====== BLAS was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${BLAS_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${BLAS_INCLUDE_DIRS}")
endif()

find_package(LAPACK REQUIRED)
if("${LAPACK_FOUND}")
    message(" ====== LAPACK was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${LAPACK_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${LAPACK_INCLUDE_DIRS}")
endif()

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
if("${HDF5_FOUND}")
    message(" ====== HDF5 was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${HDF5_LIBRARIES}")
    target_link_libraries(ql-balance_lib PUBLIC "${HDF5_C_LIBRARIES}")
    target_link_libraries(ql-balance_lib PUBLIC "${HDF5_HL_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${HDF5_INCLUDE_DIRS}")
endif()

target_sources(ql-balance_lib PUBLIC
    ${libs_root_path}/SuperLU_4.1/FORTRAN/c_fortran_zgssv.c
    ${libs_root_path}/SuperLU_4.1/FORTRAN/c_fortran_dgssv.c)
target_include_directories(ql-balance_lib PUBLIC ${libs_root_path}/SuperLU_4.1/SRC/)

target_link_libraries(ql-balance_lib PUBLIC slatec superlu cvode nvecserial gfortran)

add_executable(ql-balance ${full_balance_src}/ql-balance.f90)
target_link_libraries(ql-balance PUBLIC ql-balance_lib)

