#!/bin/bash

mkdir external
cd external
echo "Building Slatec..."
mkdir slatec
cd slatec
mkdir lib
curl -L https://www.netlib.org/slatec/slatec_src.tgz -o - | tar xz

CFLAGS="-O2"
gfortran -c -Wall -Wtabs -mtune=generic $CFLAGS src/*.f
ar rcs libslatec.a *.o
rm *.o
mv libslatec.a lib/
echo "Finished building Slatec..."

cd ../ # external
echo ""


echo ""
echo "Building Sundials..."
git clone git@github.com:LLNL/sundials.git
cd sundials
git checkout v5.7.0
mkdir build
cd build
cmake ..
make
cd ../../
echo "Finished building Sundials..."

echo "Building lapack..."
curl -L http://www.netlib.org/lapack/lapack-3.2.1.tgz -o - | tar xz
cd lapack-3.2.1
cp INSTALL/make.inc.gfortran make.inc
make -j lib
make -j blaslib
make clean
cd ..
echo "Finished building lapack..."
echo ""

echo "Building gsl-2.4"
curl -L https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz -o - | tar xz
cd gsl-2.4
./configure
make -j
cd ..
echo "Finished building gsl-2.4..."
echo ""

cd ../KiLCA
mkdir build
cd build
cmake ..
make -j
