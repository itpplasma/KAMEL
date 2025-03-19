#!/bin/bash

OLD_DIR=$(pwd)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd $SCRIPT_DIR
./KiLCA/kilca.sh || exit $?
cd QL-Balance
make || exit $?
cd $OLD_DIR
