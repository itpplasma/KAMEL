# Owned CI controller handoff — 21 October 2026

Draft [PR 361](https://github.com/itpplasma/KAMEL/pull/361) is sealed at
fc943377b329b4e0ca6066b68c30317a89a5465b on
`ci/owned-capacity-20261021`. This documentation branch preserves that exact
candidate and its current qualification. Numerical source is unchanged.

All PR events skip required jobs before runner allocation. Merge groups,
explicit dispatches and direct default-branch nonqueue pushes validate on
owned capacity. Existing required context names, Python 3.10/3.11 coverage,
native build targets, harness tests, periodic integration test, full CTest
commands and independent golden workflow remain. `full=true` selects all
registered CTests; it does not invoke the separate golden workflow.

Linux uses owned Ubuntu24 runners, two CPUs/16 GiB. Build parallelism is two;
BLAS/OpenMP use one thread and CTest remains sequential. The golden build
helper respects a positive allocated CMAKE_BUILD_PARALLEL_LEVEL (default two)
instead of unbounded `-j`; all three original targets remain. Both golden
jobs retain 30-minute timeouts, original baselines and tolerances; the
pre-fortnum diagnostic remains advisory.

Required ARM jobs use the owned native ARM host (actual macOS26.6.2,
deployment minimum14.0). Legacy macOS14 context names remain stable; they do
not mean macOS14 was executed. The provisioned native helper
`/opt/itp-ci-macos/setup-python.sh` creates job-private Python3.11 environments;
the test job requests the build job's exact Python version. Homebrew
dependencies are checked read-only. Tests rebuild in their own isolated
workspace so CMake/CTest/library paths cannot point into another job.

Current actual qualification runs:

- [CI 39645141209](https://github.com/itpplasma/KAMEL/actions/runs/39645141209),
  explicit `full=true`, has started owned Pre-commit; native builds and full
  CTest remain pending.
- [Golden 39645145813](https://github.com/itpplasma/KAMEL/actions/runs/39645145813)
  is independently queued on owned capacity.
- PR runs 39645151518 and 39645151501 skipped every job with null runner IDs.

Static/local checks passed: actionlint, Bash/ShellCheck, configured staged
hooks, 38 existing comparator/configuration tests, and an actual-CMake probe
using synthetic targets (default/explicit budgets pass; zero/malformed values
are rejected). The six-minor native Python helper pilot is infrastructure
evidence, not KAMEL build, CTest or MPI qualification.

Controller source approval is conditional on actual owned ARM complete build
and full CTest success. Retain raw native failures and the independent golden
gate. Do not mark ready or claim qualification from queue waiting or local
Python tests. No replacement runs are authorized by this handoff. Existing
hygiene [PR 360](https://github.com/itpplasma/KAMEL/pull/360) stays held until
361 lands; only then rebase and pin both shared fields to
104634c01aa669d542e8a26a34a50b6471fc19ce without changing the native gates.
The controller owns promotion and global rulesets.
