cmake_minimum_required(VERSION 3.7...3.21)



project(ql-light)
enable_language(Fortran)


if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    # macOS specific configurations or commands
    message("Building on macOS")
    set(CMAKE_C_COMPILER gcc-13)
    set(CMAKE_CXX_COMPILER g++-13)
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


set(src ${PROJECT_SOURCE_DIR}/src)

set(hdf5_tools ${PROJECT_SOURCE_DIR}/hdf5_tools/)

set(light_lib_fortran_sources
    ${PROJECT_SOURCE_DIR}/../ql-balance/src/wave_code_data_64bit.f90
    ${PROJECT_SOURCE_DIR}/../ql-balance/src/hdf5_tools/hdf5_tools.f90
    ${PROJECT_SOURCE_DIR}/src/h5mod.f90
    )

add_library(ql-light-lib STATIC ${light_lib_fortran_sources})

add_subdirectory(${src})


find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
if("${HDF5_FOUND}")
    message(" ====== HDF5 was found! ======")
    target_link_libraries(ql-light-lib PUBLIC "${HDF5_LIBRARIES}")
    target_link_libraries(ql-light-lib PUBLIC "${HDF5_C_LIBRARIES}")
    target_link_libraries(ql-light-lib PUBLIC "${HDF5_HL_LIBRARIES}")
    target_include_directories(ql-light-lib PUBLIC "${HDF5_INCLUDE_DIRS}")
endif()

add_executable(ql-light ${src}/ql-light.f90)
target_link_libraries(ql-light PUBLIC ql-light-lib)