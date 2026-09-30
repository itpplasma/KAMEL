# KIM Python API

`kamel-kim` provides the supported Python API and command-line interface for the KIM
plasma-response solver. The scientific implementation remains in the `KIM.x` Fortran executable;
this package validates inputs, stages reproducible runs, invokes the executable, reads structured
results, and runs one-dimensional parameter scans.

## Requirements

- Python 3.10 or newer
- A separately built `KIM.x` executable for simulation runs

## Development installation

From the KAMEL repository root:

```bash
python -m pip install -e './KIM/python[test]'
```

Check the installed command:

```bash
kim --help
```

## CLI overview

The CLI accepts a complete JSON request or a supported KIM namelist:

```bash
kim validate request.json
kim run request.json --executable /path/to/KIM.x --runs-dir runs
kim sweep request.json --parameter periodic.n_rg --values 512 --values 1024
kim sweep request.json --scale-profile Er --values -1 --values 0 --values 1
kim status RUN_ID --runs-dir runs
kim inspect RUN_ID --runs-dir runs --format json
kim result RUN_ID --runs-dir runs --list
```

Use `kim parameters --format json-schema` to obtain the structured configuration schema used by
the API and future automation tools.

See the [request JSON guide](docs/request-json.md) for a complete field description, profile
requirements, and validated examples for every supported run type. Copyable requests are available
under [`examples/`](examples/).

## HTCondor scan API

HTCondor orchestration is a Python API; it does not add or change a `kim` CLI command. Prepare a
validated `SweepSpec`, then explicitly stage, submit, wait, and collect it:

```python
import sys
from pathlib import Path

from kim import (
    JparCurrentMetric,
    KimCondorPlan,
    ProfileScale,
    SimulationConfig,
    SweepSpec,
    collect_condor_sweep,
    stage_condor_sweep,
    submit_condor_sweep,
    wait_condor_sweep,
)
from kim.condor_client import CondorToolConfig

# The request must name real, validated profiles accessible from this driver.
config = SimulationConfig.model_validate_json(
    Path("request-periodic.json").read_text(encoding="utf-8")
)
spec = SweepSpec(base=config, variation=ProfileScale(values=(-1.0, 0.0, 1.0)))
root = Path("/shared/kamel/runs/er-scan").resolve()

api_plan = KimCondorPlan(
    backend="kamel_kim_python",
    executable=Path("/shared/kamel/bin/KIM.x"),
    python_executable=Path(sys.executable),
    condor=CondorToolConfig(condor_bin_directory=Path("/path/to/condor/bin")),
    shared_filesystem_prefixes=(Path("/shared/kamel"),),
    should_transfer_files="NEVER",
    request_cpus=1,
    omp_threads_per_process=1,
    job_timeout_s=3600.0,
    poll_interval_s=30.0,
    max_wall_clock_s=21600.0,
    # Example explicit scalar choice; omit both fields to retain complex values only.
    api_current_component="imag",
    api_current_unit="statA",
)

staged = stage_condor_sweep(root, spec=spec, plan=api_plan)
preview = submit_condor_sweep(root, plan=api_plan, dry_run=True)
submitted = submit_condor_sweep(root, plan=api_plan)
statuses = wait_condor_sweep(root, plan=api_plan)
results = collect_condor_sweep(root, plan=api_plan)
```

Replace the paths and resource values with the local site configuration. `dry_run=True` renders the
submit descriptions without invoking Condor; it does not submit jobs. Staging refuses to reuse an
existing job directory. Each scan point runs independently, and the worker applies the configured
solver timeout and `OMP_NUM_THREADS`. The driver deadline removes pending clusters and records the
outcome in `condor_status.json`. Inspect that status document and the per-job
`condor_run_record.json` files when a job fails; collection requires a terminal successful scheduler
status and successful worker record for every point.

Choose a backend explicitly:

- **`kamel_kim_python`** runs one `kamel-kim` API sweep point per job. The execute-node Python
  interpreter must have `kamel-kim` and its dependencies installed. If using a source checkout
  instead, set `kim_source_path` to its shared `KIM/python/src` directory. `getenv` is disabled, so
  do not rely on the submitter's `PYTHONPATH` or other environment being copied to workers.
- **`kim_x_namelist`** runs the staged `KIM.x` with a reviewed base namelist. It copies the source
  profiles per point and applies the scale only to the staged `Er.dat`; it does not rewrite the
  reviewed physics settings. When selecting a text current-column metric, the reviewed namelist must
  use `electrostatic`, `electromagnetic`, `flr2`, or `flr2_benchmark`, set
  `kim_io.hdf5_output = .false.`, and have a `kim_config.collision_model` matching the metric. KIM
  then writes `fields/jpar.dat`; column 0 is
  radius and columns 1–3 are the real, imaginary, and magnitude values written by KIM. The
  metric's zero-based `current_column` selects one column, which is trapezoid-integrated against the
  radius grid. `current_unit` is the caller-declared unit of that integral; this operation does not
  add a `2*pi*r` factor or perform a unit conversion. Periodic namelist runs write their fields to
  HDF5 and cannot use this text-column metric.

  For example, after reviewing the output-column meaning and the resulting integral unit:

  ```python
  namelist_metric = JparCurrentMetric(
      current_column=1,  # Example only: the real column in KIM's text field output.
      current_unit="replace-with-reviewed-integral-unit",
      collision_model="FokkerPlanck",  # Must match the reviewed namelist.
  )
  namelist_plan = KimCondorPlan(
      backend="kim_x_namelist",
      executable=Path("/shared/kamel/bin/KIM.x"),
      python_executable=Path(sys.executable),
      condor=CondorToolConfig(condor_bin_directory=Path("/path/to/condor/bin")),
      shared_filesystem_prefixes=(Path("/shared/kamel"),),
      should_transfer_files="NEVER",
      base_namelist=Path("/path/to/reviewed/KIM_config.nml"),
      jpar_current_metric=namelist_metric,
  )
  ```

  Replace the example column/unit/model and paths with independently reviewed values. Pass
  `namelist_plan` to the same stage-submit-wait-collect calls shown above.

The staged run tree, solver executable, and worker Python runtime must be visible to execute nodes at
the same paths. The adapter currently relies on those shared paths; setting `should_transfer_files`
does not package the staged profiles, job input, or outputs. With `NEVER`,
`shared_filesystem_prefixes` is checked before staging. Avoid environment secrets: workers run with
`getenv` disabled and the manifests do not copy the driver's environment.

Current interpretation is deliberately opt-in. API periodic results always retain separate complex
real and imaginary current values. To request a scalar curve, declare `api_current_component` as
`real`, `imag`, or `magnitude` and supply its unit; otherwise the resonance result is
`insufficient_scalars` with no selected scale factor. For namelist runs, supply a reviewed
`JparCurrentMetric`; without it, collection does not select a scalar minimum. The selected minimum is
over the declared scalar values and does not infer a metric from complex current data.

## Tests

The unit suite is self-contained and uses a fake executable. It does not require a Fortran build:

```bash
python -m pytest tests/unit
```

CI runs this suite in a separate Python 3.10 job so the package's minimum supported Python version
and its orchestration behavior are checked independently of the Fortran toolchain.

The periodic integration reference runs a compact parabolic-profile case against a real `KIM.x`.
It is opt-in so ordinary unit tests do not depend on a compiled solver:

```bash
KIM_RUN_INTEGRATION=1 \
KIM_EXECUTABLE=/absolute/path/to/KIM.x \
python -m pytest tests/integration/test_periodic_run.py -v
```

CI enables this test in the existing Fortran build job, using the executable produced at
`build/install/bin/KIM.x`. Running `python -m pytest` locally executes the unit suite and skips the
real-executable case unless both integration environment variables are set.

The reference checks mandatory HDF5 datasets, field shapes, finite nonzero fields, the derived
Fourier mode count, the q-crossing radius, and both campaign-compatible parallel-current
integrals. Integral tolerances allow small compiler and linear-algebra differences while still
detecting scientific drift.

The `as_is` integral uses the campaign convention: trapezoidal quadrature over grid samples whose
centres lie inside the requested as-is interval. It does not interpolate fractional cells at the
interval boundaries. The `full_window` integral uses every stored point of the endpoint-exclusive
periodic grid.

Update `tests/fixtures/periodic_reference.json` only after intentionally changing KIM physics or
numerics. Run the case with one OpenMP thread, inspect the complete HDF5 schema and logs, record the
new values, justify the change in review, and keep tolerances no wider than cross-platform evidence
requires.

## Legacy interface migration

This package replaces the former `KIMpy` runner, `KIMData` reader, and Tk-based KIM GUIs. Use
`Simulation` or `kim run` instead of creating executable symlinks and changing the process working
directory. Supply profiles explicitly; the API copies them into every run directory rather than
linking mutable external directories. Use `Result`, `kim inspect`, and `kim result` for production
HDF5 output instead of converting legacy text field files.

The experimental dispersion-relation, Poisson, and Poisson-Ampere analysis modules remain in the
top-level KAMELpy distribution because they perform separate scientific analysis and are not
replaced by this orchestration API.
