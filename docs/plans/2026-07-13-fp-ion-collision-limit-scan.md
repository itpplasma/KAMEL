# FP Ion-Collision Limit Scan Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reproducibly scan the FP ion collision rate toward zero while preserving the physical FP electron response.

**Architecture:** Add one positive ion-only collision multiplier to KIM and apply it before any collision-dependent ion background is derived. Add a benchmark scan driver and analyzer that preserve every run and compare the sequence with both FP and analytical collisionless endpoints.

**Tech Stack:** Fortran 2008, CMake/CTest, Python 3 `unittest`, GNU Make, Matplotlib.

---

### Task 1: Specify species-selective scaling

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_assembly.f90`
- Modify: `KIM/src/background_equilibrium/species_mod.f90`

1. Add assertions that an ion-frequency scaling helper multiplies ions and leaves electrons unchanged.
2. Build the focused test and verify it fails because the helper is absent.
3. Implement the minimal pure helper and use it when assigning computed collision frequencies.
4. Rebuild and verify the focused test passes.

### Task 2: Configure and record the multiplier

**Files:**
- Modify: `KIM/src/setup/config_mod.f90`
- Modify: `KIM/src/general/read_config.f90`
- Modify: `KIM/src/general/config_display.f90`
- Modify: `KIM/src/util/IO_collection.f90`
- Modify: `KIM/tests/test_collisionless_ion_assembly.f90`

1. Write a failing configuration assertion for an explicit multiplier.
2. Add `ion_fp_collision_scale=1` to `KIM_CONFIG`, require it to be positive, display it, and write it to HDF5 provenance.
3. Run the focused assembly and kernel tests.

### Task 3: Build the reproducible scan harness

**Files:**
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/tests/test_benchmark_scripts.py`
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/scripts/run_kim.py`
- Create: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/scripts/scan_fp_limit.py`
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/Makefile`

1. Write failing Python tests for scale rendering, scan-case naming, and scale parsing.
2. Add the multiplier to rendered configs and metadata.
3. Implement isolated scan directories and a Make target using the fixed scale list.
4. Run the Python suite.

### Task 4: Analyze and plot the limiting sequence

**Files:**
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/tests/test_benchmark_scripts.py`
- Create: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/scripts/analyze_fp_limit.py`
- Modify: `/Users/markusmarkl/1_projects/KIM_collisions_benchmark/README.md`

1. Write failing tests for endpoint-distance metrics.
2. Implement JSON/CSV output and a PDF with one row per mode.
3. Document the scan command, values, and interpretation.
4. Run the Python suite and focused CTest targets.

### Task 5: Execute and verify the scan

1. Build KIM once.
2. Run scale `0.3, 0.1, 0.03, 0.01` for both modes, preserving scale-one and analytical collisionless endpoints.
3. Generate the scan metrics and PDF.
4. Run `make check`, inspect the PDF, and report any low-collision numerical failure explicitly.
