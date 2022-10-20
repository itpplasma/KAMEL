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

set(src ${PROJECT_SOURCE_DIR}/src)
set(hdf5_tools ${PROJECT_SOURCE_DIR}/hdf5_tools/)
set(objs ${PROJECT_SOURCE_DIR}/OBJS/)
set(lib ${PROJECT_SOURCE_DIR}/lib/)
set(soft ${PROJECT_SOURCE_DIR}/soft/)
set(ql-balance_lib "${lib}libql-balance.a")

set(CMAKE_Fortran_MODULE_DIRECTORY "${objs}")

add_library(ql-balance_lib "")
add_subdirectory(${src})
#add_subdirectory(hdf5_tools)

find_package(SuiteSparse REQUIRED)
if("${SuiteSparse_FOUND}")
    message(" ====== SuiteSparse was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${SUITESPARSE_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${SUITESPARSE_INCLUDE_DIRS}")
endif()

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
if("${HDF5_FOUND}")
    message(" ====== HDF5 was found! ======")
    target_link_libraries(ql-balance_lib PUBLIC "${HDF5_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${HDF5_INCLUDE_DIRS}")
endif()


find_package(Sundials REQUIRED)
if("${SUNDIALS_FOUND}")
    message(" ====== Sundials was found! =====")
    target_link_libraries(ql-balance_lib PUBLIC "${SUNDIALS_LIBRARIES}")
    target_include_directories(ql-balance_lib PUBLIC "${SUNDIALS_INCLUDE_DIRS}")
endif()


add_executable(ql-balance "${src}/ql-balance.f90")

target_link_libraries(ql-balance ql-balance_lib)