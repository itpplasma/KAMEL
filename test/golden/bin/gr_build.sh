#!/bin/bash
# Build a reference (ref) or candidate (cur) binary from a git ref of the code repo.
# Usage: gr_build.sh <ref|cur> <git-ref>
set -eu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
. "$ROOT/config.sh"
TAG="$1"; GITREF="$2"
: "${GR_CODE_SRC:?set GR_CODE_SRC in config.sh}"
WT="$ROOT/build_$TAG"
rm -rf "$WT"; mkdir -p "$WT"
git -C "$GR_CODE_SRC" worktree prune
git -C "$GR_CODE_SRC" worktree add -f "$WT/src" "$GITREF"
( cd "$WT/src" && eval "$GR_BUILD_CMD" )
ln -sf "$WT/src/$GR_BUILD_EXE" "$WT/$GR_EXE"
echo "built $TAG ($GITREF) -> $WT/$GR_EXE"
