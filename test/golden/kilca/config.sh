# Per-code golden-record config (KiLCA).
GR_CODE=kilca
GR_EXE=KiLCA_Normal.x                # placeholder; resolved to the version-suffixed exe by gr_build_all symlink
GR_INPUT=background.in               # case auto-discovery key
GR_CODE_SRC=${GR_CODE_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
GR_BUILD_CMD='cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build -j --target kilca_normal_exe'
GR_BUILD_EXE=build/install/bin/KiLCA_Normal.x
GR_COMPARE=gr_numcompare.py

# Resonance-chaotic outputs excluded from the strict bar (see gr_numcompare.py
# GR_EXCLUDE, and docs/flre-golden-record-investigation.md). The flre_m6n2
# (m=6, n=2) FLRE zone-0 solution passes through its resonant layer, where the
# conductivity is built from a confluent-hypergeometric 1F1m evaluation with
# b ~ z: the subtraction 1F1 - 1 - z/b cancels to ~1e-11 of its operands, so its
# low bits are compiler-codegen-dependent cancellation garbage. Any independent
# implementation (gfortran port vs the pre-port C++ oracle) lands on adjacent
# floating-point values (~2e-12), and that seed amplifies through the stiff
# resonant FLRE integration to O(1) in the fields near the resonant surface,
# changing the adaptive output grid. CONTROL EXPERIMENT on the unmodified C++
# oracle: multiplying its own 1F1m result by (1 + 1e-12) reproduces the identical
# failure (EB.dat max_rel 1.6e-6, zone_0_poy_test_err O(1)), proving this is a
# pre-existing property of this ill-conditioned test case, not a port regression.
# The bar (rtol 1e-7, atol 1e-12) is unchanged; only these two proven-affected
# FLRE zone-0 outputs are excluded, all other 56 files stay strict. See
# itpplasma/KAMEL#164 and the follow-up issue.
export GR_EXCLUDE='linear-data/m_6_n_2_flab_[1,0]/EB.dat;linear-data/m_6_n_2_flab_[1,0]/zone_0_poy_test_err.dat'
