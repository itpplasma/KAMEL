#!/bin/bash
# Build the WHOLE KAMEL repo at a git ref once, into a shared build dir with
# stable exe symlinks for all three codes.
# Usage: gr_build_all.sh <ref|cur> <git-ref>
#
# Target names (logical) differ from the installed binary names:
#   KIM        target KIM_exe          -> build/install/bin/KIM.x
#   QL-Balance target ql-balance.x     -> build/install/bin/ql-balance.x
#   KiLCA      target kilca_normal_exe -> build/install/bin/KiLCA_Normal_*_64bit
set -eu
TAG="$1"; GITREF="$2"
GOLDEN_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # test/golden
SRC_REPO="$(cd "$GOLDEN_ROOT/../.." && pwd)"                     # repo root
WT="$GOLDEN_ROOT/build_$TAG"
rm -rf "$WT"; mkdir -p "$WT"
git -C "$SRC_REPO" worktree prune
git -C "$SRC_REPO" worktree add -f "$WT/src" "$GITREF"
( cd "$WT/src" \
  && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  && cmake --build build -j --target KIM_exe ql-balance.x kilca_normal_exe )
# No `cmake --install`: each exe target already sets RUNTIME_OUTPUT_DIRECTORY to
# build/install/bin, and the project libs (libneo, sundials, *_lib) link
# statically, so the built exes are self-contained. A full install would fail
# trying to install sundials components we deliberately didn't build (e.g.
# nvecmanyvector), and we don't need them.
BIN="$WT/src/build/install/bin"
ln -sf "$BIN/KIM.x"            "$WT/KIM.x"
ln -sf "$BIN/ql-balance.x"     "$WT/ql-balance.x"
ln -sf "$(ls "$BIN"/KiLCA_Normal_*_64bit | head -1)" "$WT/KiLCA_Normal.x"
echo "built $TAG ($GITREF) -> $WT  (KIM.x, ql-balance.x, KiLCA_Normal.x)"
