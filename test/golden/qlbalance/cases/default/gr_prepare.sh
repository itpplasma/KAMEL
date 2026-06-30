#!/bin/bash
# Prepare a QL-Balance synthetic runfolder in the current dir ($PWD == out dir).
# gr_run_case.sh exports GR_BUILD_SRC (the worktree root for this build), so the
# ref build drives ref's KiLCA and the cur build drives cur's KiLCA.
set -eu
SRC="${GR_BUILD_SRC:?gr_prepare needs GR_BUILD_SRC}"
# Drive the generator with THIS build's python tree + binaries (so ref's KAMELpy
# runs ref's KiLCA, cur's runs cur's). Use the setup_runfolder.py that
# gr_run_case copied into the run dir ($PWD), not one under $SRC: the baseline
# build predates this harness, so $SRC has no test/golden/ -- but it always has
# python/ (the KAMELpy the generator imports via GR_PYTHON_DIR).
export GR_PYTHON_DIR="$SRC/python"
export PYTHONPATH="$SRC/python:${PYTHONPATH:-}"
# Generate balance_conf.nml + profiles + run KiLCA (flre+vacuum) into $PWD.
python3 "$PWD/setup_runfolder.py" "$PWD" \
        --nml-template "$PWD/balance_conf.nml"
