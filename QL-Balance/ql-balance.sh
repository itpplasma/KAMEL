#!/bin/bash

# SuiteSparse
echo "Building SuiteSparse..."
cd external
git clone git@github.com:DrTimothyAldenDavis/SuiteSparse.git
cd SuiteSparse
make -j 4
echo "Finished building SuiteSparse..."
echo ""
cd ..

echo "Building SuperLU..."
wget https://portal.nersc.gov/project/sparse/superlu/superlu_4.1.tar.gz
tar -xzvf superlu_4.1.tar.gz
rm superlu_4.1.tar.gz
cd SuperLU_4.1
cp MAKE_INC/make.linux make.inc
# Check the operating system
if [[ "$(uname)" == "Darwin" ]]; then
    # macOS specific compiler flags
    CFLAGS="-O2"
    CFLAGS+=" -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include"
    sed -E "s|SuperLUroot\t=.*|SuperLUroot\t= $(pwd)|g" make.inc >> make.inc.tmp
    mv make.inc.tmp make.inc
    sed -E "s|BLASLIB[[:space:]]*=.*|BLASLIB = #-lblas|g" make.inc >> make.inc.tmp
    mv make.inc.tmp make.inc
    sed -E "s|RANLIB       = ranlib|RANLIB       = gcc-ranlib-12|g" make.inc >> make.inc.tmp
    mv make.inc.tmp make.inc
    sed -E "s|CC           = gcc|CC           = gcc-12|g" make.inc >> make.inc.tmp
    mv make.inc.tmp make.inc
    sed -E "s|ARCH         = ar|ARCH         = gcc-ar-12|g" make.inc >> make.inc.tmp
    mv make.inc.tmp make.inc
    sed -E "s|CFLAGS       = -O3|CFLAGS       = -O3 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -D__FLT_EVAL_METHOD__=0|g" make.inc >> make.inc.tmp
    mv make.inc.tmp make.inc
elif [[ "$(uname)" == "Linux" ]]; then
    # Linux specific compiler flags
    sed -i "s|SuperLUroot\t= \$(HOME)/Codes/SuperLU_4.1|SuperLUroot = $(pwd)|" make.inc
    sed -i 's/g77/gfortran/g' make.inc
fi
make -j
echo "Finished building SuperLU..."
echo ""
cd ..
cd ..

mkdir build
cd build
cmake -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 ..
make
