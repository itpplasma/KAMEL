#!/bin/bash
# Run ONE golden-record case with a given build. Backend-agnostic (Condor or local).
# Usage: gr_run_case.sh <build_dir> <case_dir> <out_dir>
# Reads config.sh at repo root for GR_EXE (binary name under build_dir).
set -u
BUILD_DIR="$1"; CASE_DIR="$2"; OUT_DIR="$3"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
. "$ROOT/config.sh"
# Per-case wall-clock guard so a hanging code (e.g. a ported solver that never
# converges) fails the case fast instead of hanging the whole job. Seconds;
# GR_CASE_TIMEOUT=0 disables. Uses GNU/BSD `timeout`/`gtimeout` if present.
GR_CASE_TIMEOUT="${GR_CASE_TIMEOUT:-600}"
_GR_TIMEOUT="$(command -v timeout || command -v gtimeout || true)"
gr_guard() {  # gr_guard <cmd...> — wraps in `timeout` when available
  if [ -n "$_GR_TIMEOUT" ] && [ "$GR_CASE_TIMEOUT" != 0 ]; then
    "$_GR_TIMEOUT" -k 10 "$GR_CASE_TIMEOUT" "$@"
  else
    "$@"
  fi
}
# Resolve relative paths against ROOT so cd to OUT_DIR does not break them
[[ "$BUILD_DIR" != /* ]] && BUILD_DIR="$ROOT/$BUILD_DIR"
[[ "$CASE_DIR"  != /* ]] && CASE_DIR="$ROOT/$CASE_DIR"
[[ "$OUT_DIR"   != /* ]] && OUT_DIR="$ROOT/$OUT_DIR"
EXE="$BUILD_DIR/${GR_EXE}"
[ -x "$EXE" ] || { echo "no exe: $EXE" >&2; exit 3; }
mkdir -p "$OUT_DIR" || exit 2
# self-contained run dir: copy inputs and dereference any equilibrium symlinks
cp -rL "$CASE_DIR"/. "$OUT_DIR"/ 2>/dev/null
cd "$OUT_DIR" || exit 2
# Optional input-API migration before the run (see migrations/, MIGRATIONS.md).
# Explicit GR_MIGRATE=up|down; auto 'down' when running a legacy build dir.
MIG="${GR_MIGRATE:-}"
if [ -z "$MIG" ]; then case "$(basename "$BUILD_DIR")" in legacy*|build_legacy*) MIG=down ;; esac; fi
if [ -n "$MIG" ] && [ -d "$ROOT/migrations" ]; then
  python3 "$ROOT/bin/gr_migrate.py" apply --direction "$MIG" \
    --migrations "$ROOT/migrations" --dir "$OUT_DIR" > migrate.log 2>&1 \
    || echo "gr_migrate $MIG failed (see migrate.log)" >&2
fi
# Optional per-case prepare hook (e.g. QL-Balance synthetic runfolder generation).
# Runs before the executable, with GR_BUILD_SRC pointing at this build's worktree
# root so KAMELpy uses the matching python tree + binaries.
if [ -f "$OUT_DIR/gr_prepare.sh" ]; then
  export GR_BUILD_SRC="$(cd "$BUILD_DIR" && pwd)/src"
  gr_guard bash "$OUT_DIR/gr_prepare.sh" > prepare.log 2>&1 \
    || { echo "gr_prepare failed or timed out (see prepare.log)" >&2; echo 97 > exit_code.txt; exit 97; }
fi
s=$(date +%s.%N); gr_guard "$EXE" > run.log 2>&1; rc=$?; e=$(date +%s.%N)
[ "$rc" -eq 124 ] && echo "gr_run_case: $GR_EXE timed out after ${GR_CASE_TIMEOUT}s" | tee -a run.log >&2
awk "BEGIN{printf \"%.3f\n\", $e-$s}" > runtime_seconds.txt
echo "$rc" > exit_code.txt
exit "$rc"
