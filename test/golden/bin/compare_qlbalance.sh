#!/bin/bash
# Usage: compare_qlbalance.sh <out/ref> <out/cur> [cases...]
# QL-Balance output is one HDF5 file per case (out/<case>/out/*.hdf5), not a dir
# tree, so this wrapper locates each side's HDF5 and defers to compare.py, which
# compares the curated QUANTITIES_TO_COMPARE list at the QL-Balance tolerances.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$1"; CUR="$2"
rc=0
for casedir in "$REF"/*/; do
  c=$(basename "$casedir")
  g=$(ls "$REF/$c/out/"*.hdf5 2>/dev/null | head -1)
  a=$(ls "$CUR/$c/out/"*.hdf5 2>/dev/null | head -1)
  if [ -z "$g" ] || [ -z "$a" ]; then echo "$c: missing hdf5 (ref=$g cur=$a)"; rc=1; continue; fi
  python3 "$ROOT/bin/compare.py" "$g" "$a" || rc=1
done
exit $rc
