# KIM HTCondor Runner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a bounded and resumable Condor scan API inside KAMEL's separately installable `kamel-kim` package.

**Architecture:** Put the KIM Condor client, adapter, and staged worker in `KIM/python/src/kim/`; do not import KAMELpy or the massive-shot analysis package. Stage one job per `SweepSpec` Er scale factor and keep both source backends explicit: a `kamel-kim` Python API run or `KIM.x` with a reviewed base namelist.

**Tech Stack:** Python 3.10+, `dataclasses`, Pydantic v2, `numpy`, `h5py`, HTCondor command-line tools, the existing `kim` package and `pytest`.

**Branch:** Execute this plan on `feature/kim-condor-runner`, created from `main`. Carry the approved spec and this plan onto that branch before implementation; do not base KIM code on the NEO-2 feature branch.

**Spec:** `docs/superpowers/specs/2026-09-30-kamel-condor-integration-design.md`

## Global Constraints

- Keep KIM Condor code inside `KIM/python/`; `kamel-kim` must not import KAMELpy or `massive_shot_kim`.
- Use one job per KIM scan point; preserve scan order and scale only the staged `Er.dat` for the namelist backend.
- Keep both backends explicit; do not silently choose a backend or replace reviewed solver settings.
- Make paths, resources, timeouts, polling, machine exclusions, transfer policy, and shared-filesystem prefixes explicit plan inputs.
- Version manifests and worker records; refuse overwrite and fail closed on ambiguous submit outcomes.
- Preserve complex KIM API current results; do not claim a scalar resonance minimum without an explicit validated metric.
- Never persist credentials or copy environment secrets into worker payloads.
- Keep local KIM runs, sweeps, CLI defaults, and Fortran solver behavior unchanged.

## Review Focus

- An API backend without `kamel-kim` available on the execute node must fail with a recorded launch error, not be mistaken for a solver result (Tasks 2–3 tests).
- The namelist backend must change only `Er.dat` values and declared base-namelist settings; it must not mutate source profiles (Task 2 tests).
- API backend complex currents and missing scalar metrics must be preserved/reported without a silent real/magnitude reduction (Task 5 tests).
- Ambiguous submission and restart adoption must not queue duplicate scale-factor jobs (Task 4 tests).
- Missing profiles, incompatible scan/backend config, and malformed KIM outputs must fail before scalar/resonance interpretation (Tasks 2 and 5 tests).

---

### Task 1: Package-local Condor client primitives

**Files:**
- Create: `KIM/python/src/kim/condor_client.py`
- Create: `KIM/python/tests/unit/test_condor_client.py`

**Interfaces:**
- Produces `CondorError`, frozen `CondorToolConfig(condor_bin_directory: Path = Path("/usr/bin"), command_timeout_s: float = 120.0)`, and frozen `CondorSubmitSpec` for one absolute executable and job directory.
- Produces `run_condor(arguments: Sequence[str], *, config: CondorToolConfig, tool: str) -> subprocess.CompletedProcess[str]`, `parse_condor_submit_output(text: str) -> int`, `query_job_ads(clusters: Sequence[int], *, config: CondorToolConfig, include_history: bool = True) -> tuple[dict[str, object], ...]`, `remove_clusters(clusters: Sequence[int], *, config: CondorToolConfig, reason: str) -> None`, and `verify_shared_filesystem(root, *, allowed_prefixes, should_transfer_files) -> dict[str, object]`.

- [ ] **Step 1: Write failing tests** asserting an 8-GiB request renders `request_memory = 8192`, values/quotes render safely, cluster output parses to its integer ID, queue/history JSON ads normalize, command timeout/nonzero results raise `CondorError`, and `NEVER` transfer rejects an out-of-prefix run root.
- [ ] **Step 2: Run `python -m pytest KIM/python/tests/unit/test_condor_client.py -q`** and confirm the new module is missing.
- [ ] **Step 3: Implement the minimal `kim`-local client** with bounded argument-array subprocess calls and no dependency on `neo2_for_Er`.
- [ ] **Step 4: Re-run the focused client tests** and confirm they pass.
- [ ] **Step 5: Commit** as `feat(KIM): add Condor client primitives`.

### Task 2: Typed scan plan, metric contracts, and immutable job staging

**Files:**
- Create: `KIM/python/src/kim/condor.py`
- Create: `KIM/python/tests/unit/test_condor_staging.py`

**Interfaces:**
- Produces `JparCurrentMetric(current_column: int, current_unit: str, collision_model: str)` and frozen `KimCondorPlan` with explicit `backend: Literal["kamel_kim_python", "kim_x_namelist"]`, solver/Python paths, optional KIM source path, `CondorToolConfig`, shared prefixes, CPU/OpenMP/memory, timeout/poll/deadline, exclusions, optional reviewed namelist, optional direct-output metric, and optional API-current component/unit.
- Produces `KimCondorJob(job: str, directory: Path, profile_scale_factor: float, scan_order_index: int)` and `stage_condor_sweep(root: str | Path, *, spec: SweepSpec, plan: KimCondorPlan) -> tuple[KimCondorJob, ...]`.
- Each job's `condor_job.json` contains the normalized base `SimulationConfig`, backend, exact scan index/factor, and execution bounds; the staged request points only to files in that job directory.
- The API backend requires a periodic `SweepSpec` when requesting its integrated-current metric. API current scalar selection is explicitly `real`, `imag`, or `magnitude`, with a caller-supplied unit; neither is inferred. Namelist-column metrics use zero-based file columns, with column 0 reserved for radius.

- [ ] **Step 1: Write failing tests** named `test_plan_requires_explicit_backend_and_valid_resources`, `test_stage_creates_one_job_in_scan_order`, `test_api_backend_stages_unscaled_source_profiles`, `test_namelist_backend_scales_only_er_copy`, `test_namelist_overrides_must_exist_in_reviewed_base`, and `test_stage_refuses_existing_job_directory`. For factors `(-1, 0, 1)`, assert three ordered directories; API staged Er values equal the source and namelist Er values are respectively negated/zero/original, with every non-Er profile unchanged and the source tree untouched.
- [ ] **Step 2: Run `python -m pytest KIM/python/tests/unit/test_condor_staging.py -q`** and confirm the new API is missing.
- [ ] **Step 3: Implement the plan/metric validation and staging** using existing `SweepSpec`, `ProfileSet`, and `render_namelist`-equivalent logic local to `kim`; never alter source profiles or silently insert solver-tuning keys.
- [ ] **Step 4: Run staging tests** and confirm both backends produce isolated job payloads with explicit provenance.
- [ ] **Step 5: Commit** as `feat(KIM): stage Condor scan jobs`.

### Task 3: Bounded workers for both explicit KIM backends

**Files:**
- Create: `KIM/python/src/kim/condor_worker.py`
- Create: `KIM/python/tests/unit/test_condor_worker.py`

**Interfaces:**
- Produces `main(argv: list[str] | None = None) -> int` and worker functions for the `kim-sweep` API backend and `kim-run` namelist backend.
- Each worker reads its local `condor_job.json`, applies `OMP_NUM_THREADS`, runs within the staged job directory, writes `condor_run_record.json`, and returns a meaningful exit code. API worker calls `run_sweep(..., timeout=job_timeout_s, environment={"OMP_NUM_THREADS": ...})`; namelist worker uses a bounded process group for the declared `KIM.x` executable.

- [ ] **Step 1: Write failing tests** named `test_api_worker_runs_one_profile_scale_with_timeout`, `test_namelist_worker_runs_staged_config`, `test_worker_records_nonzero_exit_and_timeout`, `test_worker_records_missing_api_dependency`, and `test_worker_rejects_missing_or_invalid_job_input`. Assert the API worker passes a singleton `ProfileScale`, exact timeout, and OMP environment to `run_sweep`; reads integrated current from `result.children[0].periodic.integrated_parallel_current(region="as_is")`; and writes distinct terminal records for nonzero exit, timeout, import/launch failure, and invalid payload.
- [ ] **Step 2: Run `python -m pytest KIM/python/tests/unit/test_condor_worker.py -q`** and confirm the worker module is missing.
- [ ] **Step 3: Implement both worker modes** with a terminal record for success, timeout, solver error, or launch failure; keep imports of `kim` inside the API-mode boundary so a missing install is an explicit recorded failure.
- [ ] **Step 4: Run worker tests** with the existing fake KIM executable and confirm the solver timeout is also passed to `run_sweep`.
- [ ] **Step 5: Commit** as `feat(KIM): add bounded Condor workers`.

### Task 4: Resumable submission and bounded queue monitoring

**Files:**
- Modify: `KIM/python/src/kim/condor.py`
- Create: `KIM/python/tests/unit/test_condor_queue.py`

**Interfaces:**
- Produces `submit_condor_sweep(root, *, plan: KimCondorPlan, dry_run: bool = False, adopt_existing: bool = True) -> dict[str, object]`.
- Produces `wait_condor_sweep(root, *, plan: KimCondorPlan) -> tuple[dict[str, object], ...]` and writes versioned `condor_status.json`.
- Persist manifest state before each submit. A failed/lost acknowledgement becomes `ambiguous`; never automatically resubmit that point. Adopt only a recorded cluster whose staged job identity still matches.

- [ ] **Step 1: Write failing tests** named `test_submit_tracks_each_scale_factor_cluster`, `test_submit_adopts_known_live_clusters`, `test_submit_fails_closed_after_ambiguous_ack`, `test_dry_run_does_not_submit`, `test_wait_distinguishes_condor_terminal_states`, and `test_deadline_removes_pending_clusters_and_persists_status`. Assert one manifest entry per scale, no second submit for an adopted cluster, ambiguous state blocks retry, dry-run invokes no tool, and deadline removals/status survive a fresh manifest read.
- [ ] **Step 2: Run `python -m pytest KIM/python/tests/unit/test_condor_queue.py -q`** and confirm queue behavior is absent.
- [ ] **Step 3: Implement manifest-first submit/adopt/poll** with bounded Condor calls and a driver deadline that records removed jobs and reasons.
- [ ] **Step 4: Run queue tests** and verify restart behavior does not create a duplicate cluster for known in-flight work.
- [ ] **Step 5: Commit** as `feat(KIM): manage Condor scan queue`.

### Task 5: Collect solver results and gate resonance metrics

**Files:**
- Modify: `KIM/python/src/kim/condor.py`
- Modify: `KIM/python/src/kim/__init__.py`
- Create: `KIM/python/tests/unit/test_condor_collection.py`
- Modify: `KIM/python/tests/unit/test_package.py`

**Interfaces:**
- Produces frozen `KimCondorScanResults(curve: tuple[dict[str, object], ...], resonance: dict[str, object], notes: tuple[str, ...])` and `collect_condor_sweep(root: str | Path, *, plan: KimCondorPlan) -> KimCondorScanResults`.
- Export `KimCondorPlan`, `JparCurrentMetric`, `KimCondorJob`, `KimCondorScanResults`, `stage_condor_sweep`, `submit_condor_sweep`, `wait_condor_sweep`, and `collect_condor_sweep` from `kim`.
- Preserve complex API current as separate real/imag values. Only calculate a minimum after the configured API scalar component or namelist current-column metric yields a complete finite curve; otherwise return `insufficient_scalars` and no selected minimum.

- [ ] **Step 1: Write failing tests** named `test_collection_preserves_complex_api_current`, `test_collection_withholds_minimum_without_metric`, `test_collection_calculates_minimum_only_with_declared_metric`, `test_collection_integrates_declared_namelist_current_column`, `test_collection_rejects_pending_missing_or_invalid_output`, and `test_package_exports_condor_api`. Assert raw API `(real, imag)` values survive exactly, no metric yields `insufficient_scalars` and a null selected factor, a declared component uses its corresponding scalar values only, and an explicit namelist column is integrated over its validated increasing radius grid.
- [ ] **Step 2: Run `python -m pytest KIM/python/tests/unit/test_condor_collection.py KIM/python/tests/unit/test_package.py -q`** and confirm the API and collection behavior are missing.
- [ ] **Step 3: Implement collection and public exports**; validate each worker record against staged backend/scale identity and retain one result record per requested scan point.
- [ ] **Step 4: Run collection/package tests** and verify no resonance factor is selected when scalar interpretation is absent, incomplete, or non-finite.
- [ ] **Step 5: Commit** as `feat(KIM): collect Condor scan results`.

### Task 6: Document the KIM Python API and run the full focused suite

**Files:**
- Modify: `KIM/python/README.md`
- Test: `KIM/python/tests/unit/test_condor_*.py`

**Interfaces:**
- Document explicit `stage_condor_sweep` → `submit_condor_sweep` → `wait_condor_sweep` → `collect_condor_sweep` usage, backend prerequisites, shared-filesystem/transfer requirements, worker Python environment, and the explicit scalar-metric behavior.
- Do not add a Condor CLI or alter existing CLI commands in this change.

- [ ] **Step 1: Add a Python API example** showing `SweepSpec`, `KimCondorPlan`, and explicit backend/metric selection without site credentials or machine-specific paths.
- [ ] **Step 2: Run `python -m pytest KIM/python/tests/unit/test_condor_client.py KIM/python/tests/unit/test_condor_staging.py KIM/python/tests/unit/test_condor_worker.py KIM/python/tests/unit/test_condor_queue.py KIM/python/tests/unit/test_condor_collection.py -q`** and confirm all focused Condor tests pass.
- [ ] **Step 3: Run `python -m pytest KIM/python/tests/unit -q`** and confirm existing KIM Python API behavior remains unchanged.
- [ ] **Step 4: Commit** as `docs(KIM): document Condor scan API`.
