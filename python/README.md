# Python modules for KiLCA/QL-Balance and KIM analysis

Production KIM orchestration is provided as the separate `kamel-kim` distribution under
`KIM/python`. The remaining `KIMpy` package contains experimental dispersion and field-solver
analysis modules.

To install:
```
pip install -e .
```

## NEO-2 HTCondor runs

NEO-2 Condor execution is an explicit KAMELpy workflow; it does not change the
existing local runner or select surfaces. Prepare the run directory with the
existing `surfaces.dat`, per-surface `neo2.in` files, and `jobs_list.txt`
contract (for example, with `stage_surfaces`). Then stage, submit, wait, and
collect:

```python
from pathlib import Path

from neo2_for_Er import (
    Neo2CondorPlan,
    collect_neo2_condor_results,
    stage_neo2_condor_jobs,
    submit_neo2_condor_jobs,
    wait_neo2_condor_jobs,
)
from neo2_for_Er.local_runner import stage_surfaces

run_directory = Path("/shared/neo2/run-001").resolve()
template_directory = Path("/shared/neo2/template")
# run_directory already contains surfaces.dat from the caller's surface selection.
stage_surfaces(run_directory, template_directory)

plan = Neo2CondorPlan(
    executable=Path("/shared/neo2/bin/neo_2_par.x"),
    shared_filesystem_prefixes=(Path("/shared/neo2"),),
    should_transfer_files="NEVER",
    request_cpus=4,
    omp_threads_per_process=4,
    request_memory_mb=30_720,
    surface_timeout_s=300,
    poll_interval_s=30,
    max_wall_clock_s=6 * 60 * 60,
)

stage_neo2_condor_jobs(run_directory, plan=plan)
submission = submit_neo2_condor_jobs(run_directory, plan=plan, dry_run=True)
print(submission["jobs"])  # dry run has no scheduler side effects
submit_neo2_condor_jobs(run_directory, plan=plan)
statuses = wait_neo2_condor_jobs(run_directory, plan=plan)
results = collect_neo2_condor_results(run_directory, plan=plan)

print(statuses)
print(results.profile)  # sorted (r_eff_cm, k_cof) pairs, or None if none succeeded
print(results.surface_records)  # one success/failure record for every surface
```

The submit host and execute nodes must share the run directory and the absolute
Python/NEO-2 executable paths. The NEO-2 adapter currently supports only
`should_transfer_files="NEVER"`; it rejects transfer modes because it does not
stage a complete Python/solver runtime bundle. Declare the shared filesystem
prefix explicitly, and keep the run directory within it. Condor tools default
to `/usr/bin`; set `condor=CondorToolConfig(...)` from
`neo2_for_Er.condor` when they are installed elsewhere.

Jobs do not inherit the submit host environment (`Getenv = false`). Ensure the
execute-node Python and NEO-2 runtime dependencies are available through the
installed runtime/RPATH or site configuration. Resource values, timeouts,
polling, and machine exclusions are explicit `Neo2CondorPlan` inputs.

Staging refuses existing Condor payload files. Use a fresh staged run directory
for a new submission. `condor_jobs.json` records cluster IDs and per-surface
identity; restart submission adopts only scheduler ads whose identity matches.
Ambiguous acknowledgements and mismatches fail closed instead of risking
duplicate work. `wait_neo2_condor_jobs` writes `condor_status.json`, distinguishes
Condor failures from worker failures, and removes jobs still pending at the
driver deadline. Collection validates worker provenance, closure periods,
`neo2_config.h5`, and finite `fulltransp.h5` `k_cof` values; individual failures
remain in `surface_records` without discarding successful profile points.
