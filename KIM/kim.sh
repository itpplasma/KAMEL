#!/bin/bash

build_zeal() {
    if test ! -d "Zeal"; then
        echo "Downloading Zeal..."
        if ! wget https://elsevier.digitalcommonsdata.com/public-files/datasets/yc9vv7rwyj/files/be55b15b-6d9c-4f0e-b5f7-09b300cfb806/file_downloaded; then
            echo "Download failed. Exiting."
            return 1
        fi
        if ! tar -xzf file_downloaded; then
            echo "Extraction failed. Exiting."
            return 1
        fi
        rm file_downloaded
        echo "Finished downloading Zeal."
    fi

    cd Zeal

    # zeal uses 'f90' => redirect to gfortran
    if test ! -f "f90"; then
        ln -s "$(which gfortran)" f90
    fi

    echo "Building dependencies..."
    if ! PATH="$PWD:$PATH" make lapack -j; then
        echo "Build failed. Exiting."
        return 1
    fi

    echo "Building Zeal..."
    if ! PATH="$PWD:$PATH" make -j; then
        echo "Build failed. Exiting."
        return 1
    fi

    cd ..

    echo "Finished building Zeal."
}

build_kim() {
    mkdir -p build
    cd build
    cmake .. || return $?
    make -j || return $?
}

mkdir -p external
cd external
build_zeal || exit $?
cd ../KIM
build_kim || exit $?
