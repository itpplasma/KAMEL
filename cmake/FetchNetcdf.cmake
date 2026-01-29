find_program(NF_CONFIG "nf-config")

if (NF_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
execute_process(COMMAND nc-config --libdir
                OUTPUT_VARIABLE NETCDFLIB_DIR)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS)
else()
message(SEND_ERROR "nf-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NETCDFINCLUDE_DIR} NETCDFINCLUDE_DIR)
string(STRIP ${NETCDFLIB_DIR} NETCDFLIB_DIR)
string(STRIP ${NETCDF_FLIBS} NETCDF_FLIBS)

message(STATUS "QL-Balance: NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "QL-Balance: NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "QL-Balance: NetCDF Fortran libs: " ${NETCDF_FLIBS})


include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})
string(REGEX MATCHALL "-L[^ ]+" NETCDF_F_PATH_ONLY "${NETCDF_FLIBS}")
string(REPLACE "-L" "" NETCDF_F_PATH "${NETCDF_F_PATH_ONLY}")
add_link_options(${NETCDF_F_PATH_ONLY})
