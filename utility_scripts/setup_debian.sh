#!/usr/bin/env bash

sudo apt update && apt install -y wget git cmake make gcc g++ gfortran \
    libnetcdf-dev libnetcdff-dev openmpi-bin libopenmpi-dev libfftw3-dev \
    python3 python3-pip python3-numpy python3-scipy ninja-build libgsl-dev \
    libopenblas-dev
