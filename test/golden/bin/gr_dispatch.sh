#!/bin/bash
# Golden-record dispatcher: run every case for one or both builds via Condor,
# local parallel, or a hybrid claim-pool, then compare ref vs cur.
# Usage: gr_dispatch.sh --backend condor|local|hybrid [--builds "ref cur"]
#                       [--local-jobs N] [--cases GLOB]
set -eu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
. "$ROOT/config.sh"
BACKEND=condor; BUILDS="ref cur"; LJOBS=24; CASEGLOB="cases/*"
while [ $# -gt 0 ]; do case "$1" in
  --backend) BACKEND="$2"; shift 2;; --builds) BUILDS="$2"; shift 2;;
  --local-jobs) LJOBS="$2"; shift 2;; --cases) CASEGLOB="$2"; shift 2;; *) shift;; esac; done
mkdir -p logs out claims
CASES=$(for d in $CASEGLOB; do [ -f "$d/${GR_INPUT}" ] && echo "${d#cases/}"; done)
: > caselist_builds.txt
for b in $BUILDS; do for c in $CASES; do echo "$c $b"; done; done >> caselist_builds.txt
echo "cases=$(echo "$CASES"|wc -w) builds=($BUILDS) backend=$BACKEND"

run_one(){ b="$1"; c="$2"; bin/gr_run_case.sh "build_$b" "cases/$c" "out/$b/$c"; }
drain(){ while read -r c b; do mkdir "claims/${b}__${c//\//_}" 2>/dev/null && run_one "$b" "$c"; done < caselist_builds.txt; }

case "$BACKEND" in
  local)  export -f run_one; export GR_EXE GR_INPUT
          xargs -P "$LJOBS" -a caselist_builds.txt -L1 bash -c 'run_one "$1" "$0"' ;;
  condor) sed -e "s#REPLACE_INITIALDIR#$ROOT#" -e "s#REPLACE_BATCH#gr-$(basename "$ROOT")#" \
              bin/condor.sub > logs/job.sub
          condor_submit logs/job.sub ;;
  hybrid) drain & DPID=$!                       # local drainer on this node
          sed -e "s#REPLACE_INITIALDIR#$ROOT#" -e "s#REPLACE_BATCH#gr-$(basename "$ROOT")#" \
              bin/condor_drain.sub > logs/drain.sub
          condor_submit logs/drain.sub          # K condor drainers share claims/
          wait $DPID ;;
esac
