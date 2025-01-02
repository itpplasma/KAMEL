#!/bin/bash
export KAMELTOPDIR=$(pwd)

echo "Building UMFPACK in SuiteSparse..."
cd external
git clone git@github.com:DrTimothyAldenDavis/SuiteSparse.git

cd SuiteSparse/AMD
make local
make install
cd ..
cd SuiteSparse_config
make local
make install
cd ..

cd SuiteSparse/UMFPACK
cd build
cmake $(CMAKE_OPTIONS) -USUITESPARSE_PKGFILEDIR -DSUITESPARSE_LOCAL_INSTALL=1 -DUMFPACK_USE_CHOLMOD=OFF ..
cmake --build . --config Release -j
make install
echo "Finished building UMFPACK in SuiteSparse..."
echo ""
cd ../../..

echo "Building SuperLU..."
git clone git@github.com:xiaoyeli/superlu.git
cd superlu
mkdir build
cd build
cmake ..
make -j
echo "Finished building SuperLU..."
echo ""
cd ..
cd ..

cd $KAMELTOPDIR/QL-Balance
mkdir build
cd build
cmake ..
make
