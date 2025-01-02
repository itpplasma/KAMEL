#!/bin/bash

echo "Building SuiteSparse..."
cd external
git clone git@github.com:DrTimothyAldenDavis/SuiteSparse.git
cd SuiteSparse
make -j 4
echo "Finished building SuiteSparse..."
echo ""
cd ..

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

cd QL-Balance
mkdir build
cd build
cmake ..
make
