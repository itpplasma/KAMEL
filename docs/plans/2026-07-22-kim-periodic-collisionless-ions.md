# Collisionless Periodic Ion Kernels Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `m_phi=0` collisionless ion charge and parallel-current kernels to the forced-periodicity electrostatic solver while preserving FP electrons and the default FP-ion path.

**Architecture:** A new `collisionless_fourier_kernel_m` implements the four local ion kernel cores from thesis Chapter 13 and exposes configured hybrid wrappers. `periodic_assembly_m` switches to those wrappers, which preserve existing FP electron and FP-ion behavior and substitute collisionless ions only when configured.

**Tech Stack:** Fortran 2008, CMake/Ninja, CTest, LAPACK, libcerf plasma-dispersion function, Wolfram Language/Mathematica, LaTeX.

---

### Task 1: Freeze symbolic formulas and test oracles

**Files:**
- Create: `docs/derivations/collisionless_periodic_kernels.wl`
- Create: `KIM/tests/test_collisionless_fourier_kernel.f90`
- Modify: `KIM/tests/CMakeLists.txt`

**Step 1: Add the Mathematica derivation source**

Encode thesis (13.53)--(13.59), set `m_phi=0`, simplify the four moments, verify the magnetic-current recurrence, and print high-precision numerical values for a nontrivial manufactured parameter set.

**Step 2: Run Mathematica**

Run:

```sh
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
  -run 'Get["docs/derivations/collisionless_periodic_kernels.wl"];Quit[]'
```

Expected: all symbolic checks print `True` or `0`, followed by four finite complex oracle values.

**Step 3: Write failing pure-kernel tests**

Test the intended public API:

```fortran
call collisionless_ion_cores(plasma, 1, kr, krp, j, epsilon, &
    rho_phi, rho_B, j_phi, j_B)
```

Assert the four Mathematica oracle values, homogeneous anchors, the `jB/rhoB` recurrence, finiteness at resonance, and epsilon sensitivity.

**Step 4: Configure and run the test to verify RED**

Run:

```sh
cmake -S . -B build
cmake --build build --target test_collisionless_fourier_kernel.x
ctest --test-dir build -R '^test_collisionless_fourier_kernel$' --output-on-failure
```

Expected: compilation fails because `collisionless_fourier_kernel_m` does not exist.

### Task 2: Implement the four collisionless ion cores

**Files:**
- Create: `KIM/src/asymptotics/collisionless_fourier_kernel.f90`
- Modify: `KIM/src/CMakeSources.in`

**Step 1: Implement shared local response inputs**

Use `flr_arg_pair_sp` and `scaled_bessel_pair`; evaluate `k_abs`, `k_pole`, real `z0`, and `plasma_Z(z0)` from the local periodic background. Reject nonpositive epsilon and invalid thermal/Debye/gyrofrequency inputs.

**Step 2: Implement `rho-Phi` and `j-Phi`**

Evaluate `a0`, `a1`, and `a2` from thesis (13.53)--(13.55) at `m_phi=0`, then the zeroth and first velocity moments from (13.58)--(13.59). Apply their respective thesis prefactors and the established `1/(8*pi^2)` continuous Fourier normalization exactly once.

**Step 3: Implement `rho-B` and `j-B`**

Use the Mathematica-reduced common magnetic bracket from thesis (13.151) and (13.153), the causal pole, and `1/c`. Evaluate both functions independently, retaining the recurrence as a test rather than as the implementation.

**Step 4: Build and verify GREEN**

Run the focused test. Expected: all oracle and identity cases pass.

### Task 3: Add configured hybrid dispatch

**Files:**
- Modify: `KIM/src/asymptotics/collisionless_fourier_kernel.f90`
- Modify: `KIM/src/electrostatic_poisson/periodic_assembly.f90`
- Test: `KIM/tests/test_collisionless_fourier_kernel.f90`
- Test: `KIM/tests/test_periodic_assembly.f90`

**Step 1: Write failing dispatch tests**

Assert:

- FP configuration reproduces the existing `hatG_*` outputs;
- collisionless configuration equals FP electron plus collisionless ion cores;
- ion current is nonzero;
- `turn_off_*` flags are honored;
- changing FP ion collision frequency does not alter collisionless ion cores.

**Step 2: Verify RED**

Run the two focused tests. Expected: the collisionless dispatch assertions fail because periodic assembly still uses all-species FP kernels.

**Step 3: Implement wrappers and assembly selection**

Expose four `configured_hatG_*` wrappers and make `periodic_assembly_m` call them. Validate the collision-model string and preserve the default FP path.

**Step 4: Verify GREEN and FP regression safety**

Run `test_collisionless_fourier_kernel`, `test_periodic_assembly`, `test_flr2_fourier_kernel`, and `test_periodic_solve`.

### Task 4: Exercise the periodic solve and output path

**Files:**
- Modify: `KIM/tests/test_kim_solver_periodic.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90` if output metadata needs clarification

**Step 1: Add a failing low-cost end-to-end case**

Configure `ion_collision_model='collisionless'`, positive epsilon, small `M`, and a resolved `n_rg`. Assert finite/nonzero `Phi` and `jpar`, and verify that the current changes when ions are disabled.

**Step 2: Verify RED, implement missing plumbing, and verify GREEN**

The only expected plumbing is dispatch visibility and output description; no new namelist variables are required.

### Task 5: Document the derivation

**Files:**
- Create: `docs/derivations/collisionless_periodic_kernels.tex`
- Modify: `KIM/nmls/README.md`

**Step 1: Write the standalone LaTeX derivation**

Document source equations, `m_phi=0` reduction, FLR definitions, all four implementable kernels, Fourier conventions, causal prescription, dimensions, identities, and limitations.

**Step 2: Build the document**

Run:

```sh
cd docs/derivations
pdflatex -interaction=nonstopmode -halt-on-error collisionless_periodic_kernels.tex
```

Expected: exit 0. Do not commit `.aux`, `.log`, or generated PDF unless explicitly requested.

**Step 3: Update configuration documentation**

State that `electrostatic_periodic` now supports FP electrons plus collisionless ions for charge and current, requires positive `collisionless_kpar_epsilon`, and retains `m_phi=0`.

### Task 6: Full verification and representative run

**Files:**
- No production files unless verification exposes a defect.

**Step 1: Run focused KIM tests**

Run all periodic and collisionless KIM tests with `ctest -R`.

**Step 2: Run the full build and CTest suite**

Run `make all` and `make test`. Expected: no new failures; the documented baseline `test_rhs_balance` failure may remain.

**Step 3: Run a low-cost collisionless periodic case**

Use the actual AUG profile input already available to the project, with modest `M` and `n_rg`. Record wall time and verify finite `Phi_m.dat` and `jpar.dat`.

**Step 4: Review the diff**

Check formatting, generated-file exclusions, branch status, and preservation of the primary checkout.
