# KAMEL NEO-2 and KIM HTCondor integration

## Goal

Port the NEO-2 and KIM HTCondor workflows from the massive-shot analysis project
into KAMEL, where they integrate with KAMEL's existing Python interfaces. Keep
the NEO-2 runner in the top-level KAMELpy distribution (`python/`) and the KIM
runner inside the independently installable KIM Python distribution
(`KIM/python/`). Develop them on separate feature branches from `main`.

This adds distributed execution without changing solver physics, profile
interpretation, surface selection, or the existing local execution paths.

## Constraints and decisions

- KAMEL has two Python distributions. `kamel-kim` must remain installable and
  usable without importing KAMELpy or adding the repository's top-level
  `python/` directory to `PYTHONPATH`.
- Therefore each distribution owns its Condor adapter, worker payload, and the
  small HTCondor command/submit helpers it needs. They share documented
  behavior and test conventions, not Python imports.
- Use one independently staged Condor job per NEO-2 surface and per KIM scan
  point. Never change the input surface selection or reorder scan parameters
  to control scheduling.
- Make resource requests, executable paths, timeouts, polling bounds, machine
  exclusions, transfer policy, and shared-filesystem prefixes explicit plan
  inputs. Do not persist credentials or copy secrets into job payloads.
- Keep manifests and worker records versioned and machine-readable. Write
  staged payloads without overwriting an existing run directory; record
  ambiguous submission outcomes and fail closed rather than risk duplicate
  jobs. On restart, adopt jobs still known to Condor when their identity can be
  matched safely.
- Do not infer a real-valued KIM shielding-current metric from a complex API
  result. Preserve raw values and only claim a scalar minimum when the caller
  supplies an explicit metric definition.

## Architecture

### NEO-2: top-level KAMELpy

Add a Condor orchestration module and standard-library worker under
`python/neo2_for_Er/`. The public orchestration flow accepts an already prepared
NEO-2 work directory using the current `surfaces.dat`, per-surface directory,
and `jobs_list.txt` contract. It stages only the Condor metadata and worker into
those jobs; it does not select surfaces or rewrite `neo2.in`.

The plan validates the absolute NEO-2 and Python executable paths, Condor tool
configuration, CPU/OpenMP consistency, run bounds, and shared filesystem
policy. Each job records its surface values and executable fingerprint, runs
NEO-2 with the requested OpenMP setting and a hard timeout, and leaves a run
record plus the established stdout/stderr logs. The driver submits, adopts,
polls, and records scheduler state under the work directory. Collection
validates worker provenance and the existing `neo2_config.h5` /
`fulltransp.h5` outputs, retains an outcome per surface, and returns successful
`(r_eff_cm, k_cof)` points without hiding failed surfaces.

The existing local runner and legacy `neo2_for_Er` class remain usable as-is;
Condor orchestration is an explicit API rather than a silent backend switch.

### KIM: standalone `kamel-kim` package

Add the Condor API and worker within `KIM/python/src/kim/`. Expose typed plans,
job/results records, and stage/submit/wait/collect operations from
`kim.condor`. Preserve the two explicit execution backends in the source
integration:

1. `kamel_kim_python`: a one-point `kim` API sweep per job, using the existing
   `SimulationConfig`, `ProfileScale`, `SweepSpec`, and `run_sweep` contracts.
2. `kim_x_namelist`: run `KIM.x` against a reviewed base namelist, applying the
   Er scale only to the staged `Er.dat` profile and only overriding declared
   namelist keys.

Do not guess between backends or silently replace solver settings. Record
backend, scale factor, profile provenance, run status, logs, and Condor
identity. Preserve complex API current results as complex data. A scalar
shielding-current curve and resonance minimum are only produced when an
explicit, validated scalar metric is supplied; otherwise results state that
the scalar/minimum is unreviewed. The KIM module remains usable through the
Python API without requiring a new Condor CLI.

### Shared scheduling behavior

Both adapters render self-contained vanilla-universe submissions and use
bounded Condor commands. When file transfer is disabled, require the run root
to be within an explicitly allowed shared-filesystem prefix. Record scheduler
status, exit code, hold reason, execute host, and wall time when available.
Driver deadlines trigger bounded removal of pending clusters and are reflected
in status records. Network/command errors must not erase previously written
manifests or completed outputs.

## Data flow

For NEO-2, callers continue to prepare physics inputs through KAMELpy and its
existing work-directory contract, then explicitly stage, submit, wait, and
collect the Condor jobs. For KIM, callers pass a validated scan and explicit
Condor plan; staging creates one isolated job directory per requested scale
factor, and the selected backend runs only that point. Both workflows keep
queue status separate from solver completion records so Condor failures and
solver failures remain distinguishable.

## Failure handling and provenance

- Validate plans and all required inputs before submission.
- Worker payloads use bounded solver execution and write a terminal record for
  success, solver exit, timeout, or launch failure.
- Keep submit acknowledgement ambiguity explicit; do not retry an ambiguous
  submission automatically.
- Treat held, removed, missing, pending, and completed-with-error Condor jobs as
  distinct states. Do not collect a pending job as a completed result.
- Record paths, selected backend, requested resources, solver executable hash,
  surface/scan identity, scheduler cluster, worker host, timing, and failure
  details, but never environment secrets.

## Testing and acceptance

NEO-2 tests live with the top-level tests; KIM Condor tests live under
`KIM/python/tests/unit`. Use temporary directories, fake solver executables,
and mocked Condor commands. No unit test submits to a real pool.

Cover plan validation, exact submit-description rendering, payload staging and
refusal to overwrite, dry runs, submit/adopt behavior, ambiguous submission,
queue-state parsing, bounded deadlines/removal, worker timeout and launch
failure records, provenance validation, successful output collection, and
partial surface failure. KIM tests additionally cover both explicit backends,
Er-only pre-scaling for the namelist backend, unchanged profile scaling for the
Python API backend, complex-current preservation, and withholding a resonance
claim when no scalar metric is configured.

Run the focused test files first, then the existing NEO-2 local-runner tests
and KIM Python unit suite. A real Condor pool run is an opt-in operational
check, not a prerequisite for the unit suite.

## Out of scope

- Changes to NEO-2 or KIM Fortran solver behavior or scientific configuration.
- Changes to local NEO-2 execution, existing KIM run/sweep semantics, or default
  CLI behavior.
- Automatic data/profile acquisition or choosing NEO-2 surfaces and KIM scan
  values.
- A common third Python distribution or cross-import between KAMELpy and
  `kamel-kim`.
- Committing live-pool run outputs, credentials, or machine-specific secrets.
