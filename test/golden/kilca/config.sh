# Per-code golden-record config (KiLCA).
GR_CODE=kilca
GR_EXE=KiLCA_Normal.x                # placeholder; resolved to the version-suffixed exe by gr_build_all symlink
GR_INPUT=background.in               # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target kilca_normal_exe'
GR_BUILD_EXE=build/install/bin/KiLCA_Normal.x
GR_COMPARE=gr_numcompare.py
