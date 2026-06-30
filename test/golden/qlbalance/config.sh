# Per-code golden-record config (QL-Balance).
GR_CODE=qlbalance
GR_EXE=ql-balance.x
GR_INPUT=balance_conf.nml            # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target ql-balance.x && cmake --install build'
GR_BUILD_EXE=build/install/bin/ql-balance.x
GR_COMPARE=compare_qlbalance.sh      # wrapper: HDF5-aware compare of the two output files
