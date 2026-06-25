#!/usr/bin/env bash
curr_loc=$PWD
libs_path=$curr_loc/../../build/install/lib
lib_kilca=$curr_loc/../../build/install/lib
f2py -c src/getIfunc.f90 src/W2_arr.f90 $curr_loc/../../KiLCA/flre/conductivity/calc_Imn_array.f90 -m susc_funcs \
    -llapack \
    -lblas \
    -L$lib_kilca -lKiLCA_Lib_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit \
    -L$libs_path/sundials/build/src/cvode/ -lsundials_cvode \
    -L$libs_path/sundials/build/src/nvector/serial/ -lsundials_nvecserial \
    -L$libs_path/slatec/ -lslatec \
    -L$libs_path/bessel/lib/ -lbessel \
    -lgsl \
    -lgslcblas \
    -lm \
    -lc \
    --verbose
