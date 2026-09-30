# NEO-2 HTCondor Runner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a bounded, resumable HTCondor runner for already-staged NEO-2 surface jobs to the top-level KAMELpy package.

**Architecture:** Keep NEO-2 Condor code in `python/neo2_for_Er/`, separate from `kamel-kim`. The submit host stages a standard-library worker into each existing surface directory; the worker executes one surface, while the driver owns submit/adopt/poll/collect state. Preserve the existing `surfaces.dat`, `jobs_list.txt`, solver inputs, and local execution behavior.

**Tech Stack:** Python 3, `dataclasses`, `subprocess`, `numpy`, `h5py`, HTCondor command-line tools, `pytest`.

**Branch:** Execute this plan on `feature/neo2-condor-runner`, based on `main`.

**Spec:** `docs/superpowers/specs/2026-09-30-kamel-condor-integration-design.md`

## Global Constraints

- Keep the NEO-2 runner in top-level KAMELpy (`python/`) and KIM in the separate `KIM/python/` distribution.
- KAMELpy must not import from `kim` or `massive_shot_kim`.
- Use one independently staged Condor job per NEO-2 surface; do not alter surface selection or `neo2.in`.
- Make paths, resources, timeouts, polling, machine exclusions, transfer policy, and shared-filesystem prefixes explicit plan inputs.
- Version machine-readable manifests and worker records; refuse overwrite and fail closed on ambiguous submit outcomes.
- Never write credentials or environment secrets into job records.
- Validate `neo2_config.h5` and `fulltransp.h5`; retain a status record for every surface, including partial failures.
- Do not change solver physics, the local NEO-2 runner, or legacy `neo2_for_Er` execution semantics.

## Review Focus

- Invalid or unsafe `jobs_list.txt` entries, surface-count mismatches, and pre-existing job payloads must be rejected before submission (Task 3 tests).
- A lost/ambiguous `condor_submit` acknowledgement must not cause duplicate work on restart (Task 4 tests).
- OpenMP-linked NEO-2 with an under-requested CPU count, timeout, and process descendants must be handled explicitly (Tasks 2–3 tests).
- Missing, malformed, or mismatched HDF5 output must fail that surface without discarding other valid points (Task 5 tests).
- Queue query failures and a deadline with jobs still pending must yield bounded cleanup and durable per-job status (Task 4 tests).

---

### Task 1: Package-local Condor command and submit primitives

**Files:**
- Create: `python/neo2_for_Er/condor.py`
- Create: `test/test_neo2_condor.py`

**Interfaces:**
- Produces `CondorError`, frozen `CondorToolConfig(condor_bin_directory: Path = Path("/usr/bin"), command_timeout_s: float = 120.0)`, and frozen `CondorSubmitSpec` with absolute `executable`/`initialdir`, resource requests, relative output paths, transfer policy, exclusions, and accounting attributes.
- Produces `run_condor(arguments: Sequence[str], *, config: CondorToolConfig, tool: str) -> subprocess.CompletedProcess[str]`, `parse_condor_submit_output(text: str) -> int`, `query_job_ads(clusters: Sequence[int], *, config: CondorToolConfig, include_history: bool = True) -> tuple[dict[str, object], ...]`, `remove_clusters(clusters: Sequence[int], *, config: CondorToolConfig, reason: str) -> None`, and `verify_shared_filesystem(root, *, allowed_prefixes, should_transfer_files) -> dict[str, object]`.

- [ ] **Step 1: Write failing tests** asserting (a) a two-CPU spec renders `request_cpus = 2`, correct escaping, and its machine exclusion; (b) `Submitted to cluster 418` parses as `418`; (c) a JSON running ad normalizes to `running`; (d) an expired command raises `CondorError`; and (e) a run root outside the allowed prefix is rejected when transfer is `NEVER`.
- [ ] **Step 2: Run `PYTHONPATH=python python -m pytest test/test_neo2_condor.py -q`** and confirm the new imports/assertions fail.
- [ ] **Step 3: Implement the package-local Condor primitives** in `python/neo2_for_Er/condor.py`; call tools with argument arrays (`shell=False`) and finite timeouts.
- [ ] **Step 4: Re-run the focused test file** and confirm all cases pass.
- [ ] **Step 5: Commit** as `feat: add KAMELpy Condor primitives`.

### Task 2: Bounded standard-library NEO-2 worker

**Files:**
- Create: `python/neo2_for_Er/condor_worker.py`
- Create: `test/test_neo2_condor_worker.py`

**Interfaces:**
- Produces `main(argv: list[str] | None = None) -> int` and `run_bounded(command: list[str], *, timeout_s: float, omp_threads: int, record_path: Path, job_input_path: Path) -> int`.
- Worker payload contract: read `condor_job.json` from the current surface directory; run the declared absolute executable; write `condor_run_record.json`; preserve solver stdout/stderr for Condor capture.

- [ ] **Step 1: Write failing tests** asserting success returns `0`, solver failure returns its code, launch failure returns `127`, timeout returns `124` and terminates its process group, `OMP_NUM_THREADS` equals the requested value, and each record includes mode/status/exit code/host/executable hash/thread count.
- [ ] **Step 2: Run `PYTHONPATH=python python -m pytest test/test_neo2_condor_worker.py -q`** and confirm the tests fail because the worker is absent.
- [ ] **Step 3: Implement the worker** using only the Python standard library; terminate the process group on timeout and write a terminal record for every launch outcome.
- [ ] **Step 4: Run the worker tests** and confirm expected exit codes and records.
- [ ] **Step 5: Commit** as `feat: add bounded NEO-2 Condor worker`.

### Task 3: Validate plans and stage one Condor payload per existing surface

**Files:**
- Create: `python/neo2_for_Er/condor_runner.py`
- Modify: `python/neo2_for_Er/__init__.py`
- Create: `test/test_neo2_condor_staging.py`

**Interfaces:**
- Produces frozen `Neo2CondorPlan` with executable/Python paths, `CondorToolConfig`, shared prefixes, transfer policy, CPU/OpenMP counts, memory, per-surface timeout, closure-period limit, poll interval, wall-clock limit, exclusions, and optional shot/time provenance.
- Produces `stage_neo2_condor_jobs(work_directory: str | Path, *, plan: Neo2CondorPlan) -> tuple[Path, ...]`.
- Reads only the existing `jobs_list.txt` and `surfaces.dat` contract. Job paths must be safe relative child directories and match each finite seven-column surface row exactly once.

- [ ] **Step 1: Write failing tests** named `test_plan_rejects_oversubscribed_openmp`, `test_stage_writes_payload_and_submit_description`, `test_stage_rejects_unsafe_or_mismatched_jobs_list`, and `test_stage_refuses_existing_payload`. Assert staging preserves each existing `neo2.in` byte-for-byte, writes exactly one payload per seven-column row, rejects `../` job names and count mismatches, and refuses a pre-existing worker/input/submit file.
- [ ] **Step 2: Run `PYTHONPATH=python python -m pytest test/test_neo2_condor_staging.py -q`** and confirm the new plan/staging APIs are missing.
- [ ] **Step 3: Implement plan validation and staging**; add the worker and versioned job input with solver hash and physical surface fields, and do not overwrite any existing artifact.
- [ ] **Step 4: Run the staging tests**; confirm malformed surfaces and filesystem-policy violations fail before any submit command is called.
- [ ] **Step 5: Commit** as `feat: stage NEO-2 Condor surface jobs`.

### Task 4: Resumable submission and bounded queue monitoring

**Files:**
- Modify: `python/neo2_for_Er/condor_runner.py`
- Create: `test/test_neo2_condor_queue.py`

**Interfaces:**
- Produces `submit_neo2_condor_jobs(work_directory, *, plan, dry_run: bool = False, adopt_existing: bool = True) -> dict[str, object]`.
- Produces `wait_neo2_condor_jobs(work_directory, *, plan) -> tuple[dict[str, object], ...]` and writes versioned `condor_status.json`.

- [ ] **Step 1: Write failing tests** named `test_submit_records_clusters_and_adopts_live_jobs`, `test_submit_fails_closed_on_ambiguous_acknowledgement`, `test_dry_run_does_not_call_condor_submit`, `test_wait_classifies_terminal_job_states`, and `test_wait_deadline_removes_pending_jobs_and_persists_status`. Assert a known live cluster is adopted without another submit, ambiguous state persists and blocks retries, dry-run records a null cluster, each Condor terminal state maps distinctly, and deadline removal is recorded in `condor_status.json`.
- [ ] **Step 2: Run `PYTHONPATH=python python -m pytest test/test_neo2_condor_queue.py -q`** and confirm the lifecycle tests fail.
- [ ] **Step 3: Implement manifest-first submission, safe adoption, durable state updates, polling, and bounded driver-deadline removal** using the Task 1 client.
- [ ] **Step 4: Run the queue tests** and confirm ambiguous submissions are not automatically resubmitted.
- [ ] **Step 5: Commit** as `feat: manage NEO-2 Condor queue lifecycle`.

### Task 5: Validate and collect per-surface solver outcomes

**Files:**
- Modify: `python/neo2_for_Er/condor_runner.py`
- Create: `test/test_neo2_condor_collection.py`

**Interfaces:**
- Produces frozen `Neo2CondorResults(profile: np.ndarray | None, surface_records: tuple[dict[str, object], ...])` and `collect_neo2_condor_results(work_directory: str | Path, *, plan: Neo2CondorPlan) -> Neo2CondorResults`.
- Validate worker identity/hash/thread records, optional `period:` closure limit, `neo2_config.h5` `settings/boozer_s`, and finite `fulltransp.h5` `k_cof`; report profile points as sorted `(r_eff_cm, k_cof)` pairs.

- [ ] **Step 1: Write failing tests** named `test_collects_valid_points_in_radius_order`, `test_collection_retains_partial_surface_failure`, `test_collection_rejects_mismatched_worker_provenance`, `test_collection_rejects_invalid_hdf5_output`, and `test_collection_marks_excess_closure_period`. Assert output pairs are sorted by `r_eff_cm`, a failed surface remains in `surface_records` while successful points remain in `profile`, altered executable/thread provenance is rejected, missing/non-finite HDF5 values are not accepted, and a reported period above the configured maximum is a failure.
- [ ] **Step 2: Run `PYTHONPATH=python python -m pytest test/test_neo2_condor_collection.py -q`** and confirm missing collection behavior.
- [ ] **Step 3: Implement per-surface validation and collection** without turning one failed surface into a loss of valid points from the same batch.
- [ ] **Step 4: Run the collection tests** and confirm each record retains a failure reason and valid points match the existing local-runner `(r_eff_cm, k_cof)` contract.
- [ ] **Step 5: Commit** as `feat: collect NEO-2 Condor surface results`.

### Task 6: Document usage and include focused tests in KAMELpy verification

**Files:**
- Modify: `python/README.md`
- Modify: `Makefile`
- Test: `test/test_neo2_condor*.py`, `test/test_neo2_local_runner.py`

**Interfaces:**
- Document explicit `stage_neo2_condor_jobs` → `submit_neo2_condor_jobs` → `wait_neo2_condor_jobs` → `collect_neo2_condor_results` usage and the pool/shared-filesystem prerequisites.
- Make the top-level `pytest` target include the NEO-2 local and Condor tests with `PYTHONPATH=python`.

- [ ] **Step 1: Update the README and Makefile target**; do not change existing local-runner defaults.
- [ ] **Step 2: Run `PYTHONPATH=python python -m pytest test/test_neo2_condor.py test/test_neo2_condor_worker.py test/test_neo2_condor_staging.py test/test_neo2_condor_queue.py test/test_neo2_condor_collection.py test/test_neo2_local_runner.py -q`** and confirm all focused tests pass.
- [ ] **Step 3: Run `make pytest`** and confirm the repository's configured Python regression target passes.
- [ ] **Step 4: Commit** as `docs: document NEO-2 Condor runner`.
