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


set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_Fortran_COMPILER "mpif90.openmpi")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/lib64 -lstdc++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -L/usr/lib64 -lstdc++ -O3 -D DOUBLE_APPEND_FORTRAN")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set (libs_root_path ${CMAKE_SOURCE_DIR}/../../libs)

set(src ${PROJECT_SOURCE_DIR}/src)
set(hdf5_tools ${PROJECT_SOURCE_DIR}/hdf5_tools/)
set(objs ${PROJECT_SOURCE_DIR}/OBJS/)
set(lib ${PROJECT_SOURCE_DIR}/lib/)
set(ql-balance_lib "${lib}libql-balance.a")

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

add_subdirectory(${src})

#add_subdirectory(hdf5_tools)

set(umfpack_path ${libs_root_path}/SuiteSparse/UMFPACK/Lib/libumfpack.a)
add_library(umfpack STATIC IMPORTED)
set_target_properties(umfpack PROPERTIES IMPORTED_LOCATION "${umfpack_path}")
include_directories("${libs_root_path}/SuiteSparse/UMFPACK/Include/")
include_directories("${libs_root_path}/SuiteSparse/UFconfig/")
include_directories("${libs_root_path}/SuiteSparse/AMD/Include/")

set (sundials_include_path ${libs_root_path}/sundials/include)
set (sundials_include_path_build ${libs_root_path}/sundials/build/include)
include_directories ("${sundials_include_path}")
include_directories ("${sundials_include_path_build}")

set (cvode_lib_fullpath ${libs_root_path}/sundials/build/src/cvode/libsundials_cvode.a)
add_library (cvode STATIC IMPORTED)
set_target_properties (cvode PROPERTIES IMPORTED_LOCATION "${cvode_lib_fullpath}")

set (nvecserial_lib_fullpath ${libs_root_path}/sundials/build/src/nvector/serial/libsundials_nvecserial.a)
add_library (nvecserial STATIC IMPORTED)
set_target_properties (nvecserial PROPERTIES IMPORTED_LOCATION "${nvecserial_lib_fullpath}")

set(superlu_path ${libs_root_path}/SuperLU_4.1/lib/libsuperlu_4.1.a)
add_library(superlu STATIC IMPORTED)
set_target_properties(superlu PROPERTIES IMPORTED_LOCATION "${superlu_path}")

set(EXTERNAL_LIBS ${EXTERNAL_LIBS} umfpack superlu cvode nvecserial)

set(balance_lib_fortran_sources
    #${src}/ql-balance.f90
    ${src}/ql-balance_mod.f90
    ${src}/rhs_balance.f90
    ${src}/wave_code_data_64bit.f90
    ${src}/Wmn.f90
    ${src}/evolvestep.f90 
    ${src}/gengrid.f90
    ${src}/localizer.f90
    ${src}/plag_coeff.f90
    ${src}/amn_of_r.f90
    ${src}/hdf5_tools/hdf5_tools.f90
    ${src}/sparse_mod.f90
    ${src}/W2_arr.f90
    ${src}/sample_matrix.f90)

set(balance_lib_cpp_sources
    ${src}/vel_integral.cpp
    ${src}/cvodeint.cpp
    ${src}/umf4_f77wrapper.c
    ${src}/umf4_f77zwrapper.c
    )

add_library(ql-balance_lib STATIC ${balance_lib_fortran_sources} ${balance_lib_cpp_sources})

set_target_properties(ql-balance_lib PROPERTIES LINKER_LANGUAGE Fortran
    OUTPUT_NAME ql-balance-lib
    ARCHIVE_OUTPUT_DIRECTORY "${lib}"
    Fortran_MODULE_DIRECTORY "${objs}")

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
if("${HDF5_FOUND}")
    message(" ====== HDF5 was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${HDF5_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${HDF5_INCLUDE_DIRS}")
endif()


add_executable(ql-balance ${src}/ql-balance.f90)
target_link_libraries(ql-balance ql-balance_lib ${EXTERNAL_LIBS})
