# KiLCA Lorentz Torque Verification Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan
> task-by-task.

**Goal:** Add an independent, source-linked CTest for KiLCA's complex-phasor Lorentz-force
torque density and cylindrical radial integration.

**Architecture:** Extract the scalar/vector mathematics embedded in
`calc_lorentz_torque_density` into small pure C++ functions. Keep the existing production data
plumbing intact, call the pure functions from production, and test them without constructing a
`flre_quants` object.

**Tech Stack:** C++11, CMake/CTest, KiLCA's existing complex-number conventions.

---

### Task 1: Register a focused KiLCA torque test

**Files:**
- Create: `KiLCA/tests/CMakeLists.txt`
- Create: `KiLCA/tests/test_lorentz_torque.cpp`
- Modify: `KiLCA/CMakeLists.txt`

**Step 1: Write the failing test**

Add cases for the time-averaged phasor force
`0.5 Re(q n E* + (j x B*) / c)`, cylindrical moment arms, species summation, and trapezoidal
integration with the cylindrical volume factor.

**Step 2: Run test to verify it fails**

Run:
`cmake -S . -B build -G Ninja && cmake --build build --target test_kilca_lorentz_torque`

Expected: FAIL because the pure torque functions are not declared.

### Task 2: Extract and route the production formula

**Files:**
- Modify: `KiLCA/flre/quants/calc_flre_quants.h`
- Modify: `KiLCA/flre/quants/calc_flre_quants.cpp`

**Step 1: Implement the minimal pure functions**

Add functions that compute the averaged Lorentz force, apply `(0, r, R0)` moment arms, and
integrate a radial density with the existing trapezoidal cylinder rule.

**Step 2: Route production through them**

Replace the inline force and moment-arm arithmetic in `calc_lorentz_torque_density` with calls to
the tested helpers, preserving inputs and output layout.

**Step 3: Run the focused test**

Run:
`cmake --build build --target test_kilca_lorentz_torque`

Then:
`ctest --test-dir build -R kilca_lorentz_torque --output-on-failure`

Expected: PASS.

### Task 3: Document conventions and remaining scope

**Files:**
- Create: `docs/verification/kilca-flr.md`

**Step 1: Record the verified mapping**

Link dissertation equations 6.17-6.18 to the production symbols. State CGS units, the
`exp(-i omega t)` phasor convention, complex conjugation, the factor `1/2`, torque moment arms,
and the straight-cylinder integration assumptions.

**Step 2: Record unverified work explicitly**

List the generated FLR1/FLR3/FLR5 coefficient identities and constant-psi/layer-width reductions
as pending rather than implying that the whole issue is verified.

### Task 4: Verify and commit

**Files:**
- Verify all changed files.

**Step 1: Run formatting and static checks**

Run: `git diff --check`

Expected: no output.

**Step 2: Run focused and full CTest suites**

Run:
`cmake --build build`

Then:
`ctest --test-dir build -R kilca_lorentz_torque --output-on-failure`

Then:
`ctest --test-dir build --output-on-failure`

Expected: all tests pass.

**Step 3: Commit**

Run:
Stage the KiLCA CMake/test changes, the two `calc_flre_quants` files, and the plan and
verification documents.

Run: `git commit -m "test(KiLCA): verify Lorentz torque convention"`
