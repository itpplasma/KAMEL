# Collisionless Causal Pole Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the collisionless ion kernel's mixed endpoint/cell-average `k_parallel` treatment with an explicit causal signed pole and a consistently broadened magnitude while retaining the singular limit.

**Architecture:** Add a required positive inverse-length regulator to `KIM_CONFIG`.  Use `k_parallel_cc + i epsilon` for analytical signed poles and `sqrt(k_parallel_cc**2 + epsilon**2)` for analytical magnitude denominators such as `z0`; leave FP and ordinary Krook behavior unchanged.

**Tech Stack:** Fortran 2008, CMake/CTest, Python benchmark runner and `unittest`, GNU Make.

---

### Task 1: Specify the causal pole with failing unit tests

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_kernel.f90`
- Modify: `KIM/tests/CMakeLists.txt`

**Steps:**
1. Extend the focused kernel test to request a cell-centred causal pole and assert `kp_cc + i*epsilon`.
2. Assert that collisionless `z0` is computed from cell-centred `om_E`, `vT`, and the broadened magnitude.
3. Assert both inverse forms grow when epsilon decreases at `kp_cc=0`.
4. Build and run `test_collisionless_ion_kernel`; verify compilation or assertions fail because the API is absent.

### Task 2: Implement the causal pole in every collisionless prefactor

**Files:**
- Modify: `KIM/src/kernels/Krook_kernel_plasma_prefacs.f90`
- Modify: `KIM/tests/test_collisionless_ion_kernel.f90`

**Steps:**
1. Add pure helpers returning `cmplx(kp_cc, epsilon)` and `sqrt(kp_cc**2 + epsilon**2)`.
2. Change collisionless `z0` to use cell-centred background values and the magnitude helper.
3. Route `G1/G2/G3 rho_phi` through the magnitude helper and signed magnetic inverse-k terms through the causal pole.
4. Add a realistic resonance regression proving `plasma_Z(z0)` stays finite.
5. Preserve the existing collisional calculations byte-for-byte in their branches.
6. Build and run the focused test; verify it passes.

### Task 3: Make epsilon explicit and validated configuration

**Files:**
- Modify: `KIM/src/setup/config_mod.f90`
- Modify: `KIM/src/general/read_config.f90`
- Modify: `KIM/src/general/config_display.f90`
- Modify: `KIM/src/util/IO_collection.f90`
- Modify: `KIM/tests/test_collisionless_ion_assembly.f90`

**Steps:**
1. Add `collisionless_kpar_epsilon` with a non-positive default.
2. Add it to `KIM_CONFIG`, display output, and HDF5 provenance.
3. Reject collisionless ion configuration unless it is positive.
4. Add it to the assembly test namelist and assert it was read.
5. Build and run the focused configuration/assembly tests.

### Task 4: Make the benchmark reproducible

**Files:**
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/scripts/run_kim.py`
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/tests/test_run_kim.py`
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/Makefile`
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/README.md`

**Steps:**
1. Add a deterministic default epsilon to generated configurations and provenance metadata.
2. Add Python tests for the generated namelist value.
3. Add Make targets for an explicit epsilon scan without changing FP cases.
4. Document units, causal convention, and resolution requirement.
5. Run the Python unit suite.

### Task 5: Verify numerical behavior

**Files:**
- Regenerate: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/runs/*`
- Regenerate: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/results/*`

**Steps:**
1. Run focused CTest targets and the Python suite.
2. Run the original and quarter-cell-shifted m=7 FP and collisionless cases at the same epsilon.
3. Verify the collisionless operator no longer changes by roughly three orders of magnitude at fixed resolved epsilon.
4. Run at least two decreasing epsilon values and verify growth toward the collisionless singularity.
5. Run `make check`, regenerate the PDF, and inspect the final diff.
