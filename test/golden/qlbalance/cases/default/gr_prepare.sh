#!/bin/bash
# Prepare a QL-Balance synthetic runfolder in the current dir ($PWD == out dir).
# gr_run_case.sh exports GR_BUILD_SRC (the worktree root for this build), so the
# ref build drives ref's KiLCA and the cur build drives cur's KiLCA.
set -eu
SRC="${GR_BUILD_SRC:?gr_prepare needs GR_BUILD_SRC}"
export GR_PYTHON_DIR="$SRC/python"
export PYTHONPATH="$SRC/python:${PYTHONPATH:-}"
# Generate balance_conf.nml + profiles + run KiLCA (flre+vacuum) into $PWD.
python3 "$SRC/test/golden/qlbalance/cases/default/setup_runfolder.py" "$PWD" \
        --nml-template "$PWD/balance_conf.nml"
