#!/bin/bash

OLD_DIR=$(pwd)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd $SCRIPT_DIR
./KiLCA/kilca.sh || exit $?
./KIM/kim.sh || exit $?
./QL-Balance/ql-balance.sh || exit $?
cd $OLD_DIR
