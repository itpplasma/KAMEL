# Per-code golden-record config (KIM).
GR_CODE=kim
GR_EXE=KIM.x
GR_INPUT=KIM_config.nml             # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target KIM_exe'
GR_BUILD_EXE=build/install/bin/KIM.x
GR_COMPARE=gr_numcompare.py
