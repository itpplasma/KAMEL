# Per-Species I-Function Conservation Models Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add KiLCA-compatible per-species particle, energy, and momentum conservation selection to KIM while preserving its default energy-conserving behavior.

**Architecture:** Extend the shared QL-Balance I-function routine with an explicit conservation-model entry point and retain its legacy Boolean wrapper. KIM resolves separate electron and ion model settings, passes the selected model into each susceptibility evaluation, and records the resolved settings in its configuration output.

**Tech Stack:** Fortran, CMake, CTest, HDF5.

---

### Task 1: Verify the susceptibility formulas

**Files:**
- Create: `KIM/tests/test_ifunc_conservation_models.f90`
- Modify: `KIM/tests/CMakeLists.txt`
- Modify: `QL-Balance/src/base/getIfunc.f90`

1. Add a test that independently evaluates the KiLCA FPGEN formulas for models 0–3.
2. Confirm the test initially fails because the explicit model entry point is absent.
3. Implement the explicit entry point and keep `getIfunc` as a compatibility wrapper.
4. Verify model 1 matches the existing energy-conserving result and model 0 matches the bare result.

### Task 2: Add per-species KIM configuration

**Files:**
- Modify: `KIM/src/setup/config_mod.f90`
- Modify: `KIM/src/general/read_config.f90`
- Modify: `KIM/src/general/config_display.f90`
- Modify: `KIM/src/background_equilibrium/species_mod.f90`

1. Add electron and ion integer settings with KiLCA-compatible values: 0=N, 1=N+E, 2=N+P, 3=N+E+P.
2. Use `-1` internally as “inherit the legacy Boolean” so old namelists remain valid.
3. Validate and resolve both values after reading the namelist.
4. Pass the resolved model explicitly for each species.

### Task 3: Record and exercise the settings

**Files:**
- Modify: KIM example/test namelists and user documentation as appropriate.
- Modify: `KIM/src/util/IO_collection.f90` if resolved settings are not already captured.

1. Document the new settings and the deprecated Boolean fallback.
2. Add an integration fixture for electron model 1 and ion model 3.
3. Confirm the default 1/1 configuration reproduces the existing behavior.

### Task 4: Verify

1. Build the isolated worktree.
2. Run the formula and integration tests.
3. Run the full CTest suite and compare failures with the recorded baseline.
4. Inspect the final diff for unintended changes.
