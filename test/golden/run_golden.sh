#!/bin/bash
# Top-level golden-record driver: build ref+cur once, run every code's cases
# against both, compare. Exit nonzero if any code diverges.
# Usage: run_golden.sh [REF_GITREF] [CUR_GITREF]
#   defaults: REF=golden-baseline  CUR=HEAD
set -eu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
REF_REF="${1:-golden-baseline}"
CUR_REF="${2:-HEAD}"
CODES="kim qlbalance kilca"

echo "== building ref ($REF_REF) and cur ($CUR_REF) =="
bin/gr_build_all.sh ref "$REF_REF"
bin/gr_build_all.sh cur "$CUR_REF"

overall=0
for code in $CODES; do
  echo "== $code =="
  # point this code's build_ref/build_cur at the shared builds
  ln -sfn "$ROOT/build_ref" "$ROOT/$code/build_ref"
  ln -sfn "$ROOT/build_cur" "$ROOT/$code/build_cur"
  ( cd "$ROOT/$code" && bin/gr_dispatch.sh --backend local --local-jobs 2 ) || overall=1
  ( cd "$ROOT/$code" && bin/gr_compare.sh ) || overall=1
done

echo "== golden-record overall: $([ $overall -eq 0 ] && echo PASS || echo FAIL) =="
exit $overall
