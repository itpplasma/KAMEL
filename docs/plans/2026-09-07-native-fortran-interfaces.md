# Native Fortran interfaces implementation plan

**Goal:** Remove internal C calling conventions from the PR #161 migration while
preserving its numerical behavior, object ownership, and actual external interfaces.

**Architecture:** Import native procedures from their owning modules instead of
duplicating C interfaces. Group the external settings readers/getters into native
modules and pass Fortran strings directly. Resolve legacy Fortran callers together
with their definitions. Keep C interoperability at OS/library APIs and callbacks.
Preserve existing handle ownership during this interface-only refactor.

**Tech stack:** Fortran modules, C interoperability, CMake/CTest, Python regression
harness, gfortran on macOS ARM64 and Ubuntu.

1. Inventory binding definitions, duplicate declarations, legacy callers, and module
   dependencies. Record the genuine C boundary keep-list. Preserve the validated
   `1f655463` full-output results for an A/B comparison.
2. Convert antenna/background/output/eigenmode settings readers and getters into
   native module procedures. Update `settings_m.f90` and settings tests to use
   ordinary Fortran strings. Compile and run these tests.
3. Replace duplicate internal interfaces throughout `KiLCA/` with imports from
   implementation modules. Convert legacy external callers in conductivity,
   Maxwell equations, and `QL-Balance/` with their definitions. Resolve dependency
   cycles explicitly; do not depend on accidental compiler symbol spelling.
4. Convert internal adaptive-grid and ODE RHS callback plumbing to native procedure
   arguments. Keep the callbacks passed to C/C++ libraries interoperable.
5. Remove redundant internal C-string packing and matching unpacking. Preserve
   parsing semantics and existing numerical algorithms. Keep only documented,
   genuine C boundaries and interoperability types that still describe C objects.
6. Update standalone fakes/tests to exercise native interfaces, build all programs,
   and run all CTests and Python harness tests. Compare full-output KiLCA data with
   the saved validated result, preserving tolerances and all output files.
7. Obtain independent specification and quality review. Record retained boundaries
   and verification results, run all formatting hooks, commit the focused cleanup,
   and push to PR #161. Require fresh Ubuntu/macOS CI and strict golden validation.

Use existing behavioral regressions as the baseline for this refactor. Add focused
coverage when an interface change exposes an untested contract. Any new failure is
investigated before altering production behavior.
