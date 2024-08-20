curr_loc=$PWD
libs_path=$curr_loc/../../../libs
lib_kilca=$curr_loc/../../KiLCA/lib/
f2py -c src/getIfunc.f90 src/W2_arr.f90 $curr_loc/../../KiLCA/flre/conductivity/calc_Imn_array.f90 -m susc_funcs \
	-llapack \
	-lblas \
	-L$lib_kilca -lkilca  \
	-L$libs_path/sundials/build/src/cvode/ -lsundials_cvode \
	-L$libs_path/sundials/build/src/nvector/serial/ -lsundials_nvecserial \
	-L$libs_path/slatec/ -lslatec \
	-L$libs_path/bessel/lib/ -lbessel \
	-lgsl \
	-lgslcblas \
	-lm \
	-lc \
	--verbose
