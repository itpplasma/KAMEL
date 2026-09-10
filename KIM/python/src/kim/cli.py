"""Thin command-line adapters for the public KIM Python API."""

from __future__ import annotations

import json
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
import typer
from kim.config import ProfileConfig, SimulationConfig
from kim.errors import ConfigurationError, KimError
from kim.profiles import ProfileSet
from kim.results import Result
from kim.runs import RunManifest, RunRepository, RunStatus
from kim.simulation import Simulation
from kim.sweep import ParameterSweep, ProfileScale, SweepSpec, run_sweep
from pydantic import ValidationError

app = typer.Typer(help="Run and inspect KIM plasma simulations.", no_args_is_help=True)


class OutputFormat(str, Enum):
    TABLE = "table"
    JSON = "json"


class ParameterFormat(str, Enum):
    TABLE = "table"
    JSON_SCHEMA = "json-schema"


@app.callback()
def main(
    context: typer.Context,
    debug: bool = typer.Option(False, "--debug", help="Show Python exceptions."),
) -> None:
    """Run and inspect KIM plasma simulations."""

    context.ensure_object(dict)
    context.obj["debug"] = debug


@app.command("parameters")
def parameters_command(
    output_format: ParameterFormat = typer.Option(
        ParameterFormat.TABLE,
        "--format",
        help="Render a readable table or the complete JSON Schema.",
    ),
) -> None:
    """Show validated KIM parameters and their automation schema."""

    schema = SimulationConfig.model_json_schema()
    if output_format is ParameterFormat.JSON_SCHEMA:
        _json_output(schema)
        return
    typer.echo("PARAMETER                         TYPE          UNITS       SWEEPABLE")
    for path, details in _parameter_rows(schema):
        typer.echo(
            f"{path:<33} {details['type']:<13} " f"{details['units']:<11} {details['sweepable']}"
        )


@app.command("validate")
def validate_command(
    context: typer.Context,
    config_file: Path = typer.Argument(..., help="JSON request or supported KIM namelist."),
    profiles: Path | None = typer.Option(None, "--profiles", help="Override profile directory."),
    output_format: OutputFormat = typer.Option(OutputFormat.TABLE, "--format"),
) -> None:
    """Validate a configuration and its r_eff input profiles."""

    try:
        config = _load_config(config_file, profiles)
        report = ProfileSet.from_simulation(config).validate_for(config)
    except (KimError, ValidationError, OSError, ValueError) as error:
        _fail(context, error)
    payload = {
        "valid": True,
        "run_type": config.run.run_type.value,
        "profile_points": report.point_count,
        "radial_minimum": report.radial_minimum,
        "radial_maximum": report.radial_maximum,
        "resonance_radius": report.resonance_radius,
    }
    _render(payload, output_format)


@app.command("run")
def run_command(
    context: typer.Context,
    config_file: Path = typer.Argument(..., help="JSON request or supported KIM namelist."),
    profiles: Path | None = typer.Option(None, "--profiles", help="Override profile directory."),
    executable: Path | None = typer.Option(None, "--executable", help="Path to KIM.x."),
    runs_directory: Path = typer.Option(Path("runs"), "--runs-dir"),
    label: str | None = typer.Option(None, "--label"),
    timeout: float | None = typer.Option(None, "--timeout", min=0.0),
    omp_threads: int | None = typer.Option(None, "--omp-threads", min=1),
    output_format: OutputFormat = typer.Option(OutputFormat.TABLE, "--format"),
) -> None:
    """Execute one validated, reproducible KIM simulation."""

    try:
        config = _load_config(config_file, profiles)
        result = Simulation(
            config,
            executable=executable,
            runs_directory=runs_directory,
            label=label,
            timeout=timeout,
        ).run(environment=_execution_environment(omp_threads))
    except (KimError, ValidationError, OSError, ValueError) as error:
        _fail(context, error)
    payload = {
        "run_id": result.run_id,
        "status": result.status.value,
        "run_directory": str(result.run_directory),
        "output_file": str(result.output_file),
        "exit_code": result.manifest.exit_code,
        "failure": (
            result.manifest.failure.model_dump(mode="json")
            if result.manifest.failure is not None
            else None
        ),
    }
    _render(payload, output_format)
    if result.status is not RunStatus.SUCCEEDED:
        raise typer.Exit(3)


@app.command("sweep")
def sweep_command(
    context: typer.Context,
    config_file: Path = typer.Argument(..., help="JSON request or supported KIM namelist."),
    parameter: str | None = typer.Option(None, "--parameter"),
    scale_profile: str | None = typer.Option(None, "--scale-profile"),
    values: list[str] = typer.Option([], "--values", help="Repeat for each ordered value."),
    profiles: Path | None = typer.Option(None, "--profiles"),
    executable: Path | None = typer.Option(None, "--executable", help="Path to KIM.x."),
    runs_directory: Path = typer.Option(Path("runs"), "--runs-dir"),
    label: str | None = typer.Option(None, "--label"),
    timeout: float | None = typer.Option(None, "--timeout", min=0.0),
    omp_threads: int | None = typer.Option(None, "--omp-threads", min=1),
    stop_on_failure: bool = typer.Option(False, "--stop-on-failure"),
    output_format: OutputFormat = typer.Option(OutputFormat.TABLE, "--format"),
) -> None:
    """Run an ordered scalar or Er-profile parameter sweep."""

    try:
        if (parameter is None) == (scale_profile is None):
            raise ConfigurationError("select exactly one of --parameter or --scale-profile")
        if not values:
            raise ConfigurationError("provide at least one --values entry")
        config = _load_config(config_file, profiles)
        parsed_values = tuple(_parse_number(value) for value in values)
        if parameter is not None:
            spec: ParameterSweep | SweepSpec = ParameterSweep(
                base=config,
                parameter=parameter,
                values=parsed_values,
                continue_on_failure=not stop_on_failure,
            )
        else:
            spec = SweepSpec(
                base=config,
                variation=ProfileScale(profile=scale_profile, values=parsed_values),
                continue_on_failure=not stop_on_failure,
            )
        result = run_sweep(
            spec,
            executable=executable,
            runs_directory=runs_directory,
            label=label,
            timeout=timeout,
            environment=_execution_environment(omp_threads),
        )
    except (KimError, ValidationError, OSError, ValueError) as error:
        _fail(context, error)
    payload = {
        "sweep_id": result.sweep_id,
        "status": result.status.value,
        "manifest": str(result.manifest),
        "child_run_ids": [child.run_id for child in result.children],
        "child_statuses": [child.status.value for child in result.children],
    }
    _render(payload, output_format)
    if result.status.value != "succeeded":
        raise typer.Exit(3)


@app.command("status")
def status_command(
    context: typer.Context,
    run_id: str = typer.Argument(...),
    runs_directory: Path = typer.Option(Path("runs"), "--runs-dir"),
    output_format: OutputFormat = typer.Option(OutputFormat.TABLE, "--format"),
) -> None:
    """Show the lifecycle state of an existing run."""

    try:
        manifest = RunRepository(runs_directory).inspect(run_id)
    except (KimError, ValidationError, OSError, ValueError) as error:
        _fail(context, error)
    _render({"run_id": manifest.run_id, "status": manifest.status.value}, output_format)


@app.command("inspect")
def inspect_command(
    context: typer.Context,
    run_id: str = typer.Argument(...),
    runs_directory: Path = typer.Option(Path("runs"), "--runs-dir"),
    output_format: OutputFormat = typer.Option(OutputFormat.TABLE, "--format"),
) -> None:
    """Inspect the complete persisted manifest for one run."""

    try:
        manifest = RunRepository(runs_directory).inspect(run_id)
    except (KimError, ValidationError, OSError, ValueError) as error:
        _fail(context, error)
    payload = manifest.model_dump(mode="json")
    if output_format is OutputFormat.JSON:
        _json_output(payload)
    else:
        for name, value in payload.items():
            typer.echo(f"{name}: {_compact(value)}")


@app.command("result")
def result_command(
    context: typer.Context,
    run_id: str = typer.Argument(...),
    list_contents: bool = typer.Option(False, "--list", help="List available datasets."),
    dataset: str | None = typer.Option(None, "--dataset", help="Read a bounded dataset path."),
    runs_directory: Path = typer.Option(Path("runs"), "--runs-dir"),
    output_format: OutputFormat = typer.Option(OutputFormat.TABLE, "--format"),
) -> None:
    """List or read structured HDF5 output from an existing run."""

    try:
        if list_contents and dataset is not None:
            raise ConfigurationError("select only one of --list or --dataset")
        repository = RunRepository(runs_directory)
        manifest = repository.inspect(run_id)
        output_path = _primary_output(repository, run_id, manifest)
        result = Result(output_path)
        if dataset is None:
            payload: dict[str, Any] = {
                "run_id": run_id,
                "output_file": str(output_path),
                "datasets": list(result.list_datasets()),
            }
        else:
            payload = {
                "run_id": run_id,
                "path": dataset,
                "value": _json_value(result.read_dataset(dataset)),
                "metadata": result.dataset_metadata(dataset).__dict__,
            }
    except (KimError, ValidationError, OSError, ValueError) as error:
        _fail(context, error)
    if output_format is OutputFormat.JSON:
        _json_output(payload)
    elif dataset is None:
        typer.echo("\n".join(payload["datasets"]))
    else:
        typer.echo(_compact(payload["value"]))


def _load_config(path: Path, profiles: Path | None) -> SimulationConfig:
    try:
        if path.suffix.casefold() == ".json":
            config = SimulationConfig.model_validate_json(path.read_text(encoding="utf-8"))
        else:
            config = SimulationConfig.from_namelist(path)
        if profiles is not None:
            profile_config = ProfileConfig.model_validate(
                {**config.profiles.model_dump(mode="python"), "directory": profiles}
            )
            config = SimulationConfig.model_validate(
                {**config.model_dump(mode="python"), "profiles": profile_config}
            )
        return config
    except (ValidationError, OSError, ValueError) as error:
        raise ConfigurationError(f"could not load configuration {path}: {error}") from error


def _primary_output(repository: RunRepository, run_id: str, manifest: RunManifest) -> Path:
    outputs = manifest.discovered_outputs or manifest.expected_outputs
    if not outputs:
        raise ConfigurationError(f"run {run_id} has no declared result output")
    return repository.paths(run_id).root / outputs[0]


def _execution_environment(omp_threads: int | None) -> dict[str, str]:
    return {"OMP_NUM_THREADS": str(omp_threads)} if omp_threads is not None else {}


def _parse_number(value: str) -> int | float:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise ConfigurationError(f"invalid numeric sweep value: {value!r}") from error
    if isinstance(parsed, bool) or not isinstance(parsed, (int, float)):
        raise ConfigurationError(f"sweep value must be a JSON number: {value!r}")
    return parsed


def _fail(context: typer.Context, error: Exception) -> None:
    if context.obj and context.obj.get("debug"):
        raise error
    typer.echo(f"Error: {_error_message(error)}", err=True)
    raise typer.Exit(1)


def _error_message(error: Exception) -> str:
    validation: ValidationError | None = None
    if isinstance(error, ValidationError):
        validation = error
    elif isinstance(error.__cause__, ValidationError):
        validation = error.__cause__
    if validation is None:
        return str(error)
    details = validation.errors(include_url=False, include_context=False, include_input=False)
    messages = [
        f"{'.'.join(str(part) for part in item['loc'])}: {item['msg']}" for item in details[:6]
    ]
    if len(details) > len(messages):
        messages.append(f"and {len(details) - len(messages)} more error(s)")
    return "validation failed: " + "; ".join(messages)


def _render(payload: dict[str, Any], output_format: OutputFormat) -> None:
    if output_format is OutputFormat.JSON:
        _json_output(payload)
    else:
        for name, value in payload.items():
            typer.echo(f"{name}: {_compact(value)}")


def _json_output(value: Any) -> None:
    typer.echo(json.dumps(_json_value(value), indent=2, sort_keys=True))


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if np.issubdtype(value.dtype, np.complexfloating):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(name): _json_value(item) for name, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _compact(value: Any) -> str:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(_json_value(value), sort_keys=True)
    return str(value)


def _parameter_rows(schema: dict[str, Any]) -> list[tuple[str, dict[str, str]]]:
    definitions = schema.get("$defs", {})
    rows: dict[str, dict[str, str]] = {}

    def resolve(node: dict[str, Any]) -> dict[str, Any]:
        reference = node.get("$ref")
        if isinstance(reference, str) and reference.startswith("#/$defs/"):
            return definitions[reference.rsplit("/", 1)[-1]]
        return node

    def walk(node: dict[str, Any], prefix: str) -> None:
        node = resolve(node)
        alternatives = node.get("anyOf") or node.get("oneOf")
        if alternatives:
            for alternative in alternatives:
                if isinstance(alternative, dict) and alternative.get("type") != "null":
                    walk(alternative, prefix)
            return
        properties = node.get("properties")
        if isinstance(properties, dict):
            for name, child in properties.items():
                public_prefix = f"{prefix}.{name}" if prefix else name
                if public_prefix.startswith("run.periodic."):
                    public_prefix = public_prefix.removeprefix("run.")
                walk(child, public_prefix)
            return
        if prefix and prefix != "run.run_type":
            rows[prefix] = {
                "type": str(node.get("type", "choice")),
                "units": str(node.get("units", "")),
                "sweepable": "yes" if node.get("sweepable") else "no",
            }

    walk(schema, "")
    return sorted(rows.items())
