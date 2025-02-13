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

# Replace space by semicolon in the Fortran libs
#string(REPLACE " " ";" NETCDF_FLIBS ${NETCDF_FLIBS})

message(STATUS "QL-Balance: NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "QL-Balance: NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "QL-Balance: NetCDF Fortran libs: " ${NETCDF_FLIBS})
