"""Small HTCondor client and submit-description helpers for KAMELpy."""

from __future__ import annotations

import json
import math
import re
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

_CLUSTER_ID_RE = re.compile(r"submitted to cluster\s+(?P<cluster>\d+)", re.IGNORECASE)
_JOB_STATUS_NAMES = {
    1: "idle",
    2: "running",
    3: "removed",
    4: "completed",
    5: "held",
    6: "transferring_output",
    7: "suspended",
}
_QUERY_ATTRIBUTES = (
    "ClusterId",
    "ProcId",
    "JobStatus",
    "ExitCode",
    "RemoteHost",
    "RemoteWallClockTime",
    "NumHolds",
    "HoldReason",
    "RequestCpus",
    "RequestMemory",
    "ImageSize",
    "KAMELJobIdentity",
)
_ACCOUNTING_KEY_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


class CondorError(RuntimeError):
    """An HTCondor command, submit description, or job record is invalid."""


@dataclass(frozen=True, slots=True)
class CondorToolConfig:
    """Absolute HTCondor tool directory and bounded command timeout."""

    condor_bin_directory: Path = Path("/usr/bin")
    command_timeout_s: float = 120.0

    def __post_init__(self) -> None:
        directory = Path(self.condor_bin_directory)
        if not directory.is_absolute():
            raise CondorError("condor_bin_directory must be an absolute path")
        timeout = float(self.command_timeout_s)
        if not math.isfinite(timeout) or timeout <= 0.0:
            raise CondorError("command_timeout_s must be strictly positive and finite")
        object.__setattr__(self, "condor_bin_directory", directory)
        object.__setattr__(self, "command_timeout_s", timeout)

    def tool(self, name: str) -> Path:
        """Resolve one executable directly below the configured tool directory."""

        if not name or Path(name).name != name:
            raise CondorError(f"invalid Condor tool name: {name!r}")
        path = self.condor_bin_directory / name
        if not path.is_file():
            raise CondorError(f"Condor tool is missing: {path}")
        return path


@dataclass(frozen=True, slots=True)
class CondorSubmitSpec:
    """One vanilla-universe submit description for a staged solver job."""

    executable: Path
    initialdir: Path
    arguments: tuple[str, ...] = ()
    universe: str = "vanilla"
    request_cpus: int = 6
    request_memory_mb: int = 30_720
    request_disk_mb: int | None = None
    output_file: str = "NEO2-0.o$(Cluster)"
    error_file: str = "NEO2-0.e$(Cluster)"
    log_file: str = "NEO2-0.l$(Cluster)"
    notification: str = "never"
    getenv: bool = True
    should_transfer_files: str = "NEVER"
    run_as_owner: bool = True
    exclude_machines: tuple[str, ...] = ("faepop43",)
    requirements: tuple[str, ...] = ()
    queue_count: int = 1
    accounting: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        executable = Path(self.executable)
        initialdir = Path(self.initialdir)
        if not executable.is_absolute():
            raise CondorError("Executable must be an absolute path")
        if not initialdir.is_absolute():
            raise CondorError("Initialdir must be an absolute path")
        if self.universe not in {"vanilla", "local", "scheduler"}:
            raise CondorError(f"unsupported Condor universe: {self.universe}")
        if self.should_transfer_files not in {"NEVER", "IF_NEEDED", "ALWAYS"}:
            raise CondorError("should_transfer_files must be NEVER, IF_NEEDED, or ALWAYS")
        if self.notification not in {"never", "always", "complete", "error"}:
            raise CondorError(f"unsupported notification value: {self.notification}")
        for name in ("request_cpus", "request_memory_mb", "queue_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise CondorError(f"{name} must be a strictly positive integer")
        disk = self.request_disk_mb
        if disk is not None and (isinstance(disk, bool) or not isinstance(disk, int) or disk <= 0):
            raise CondorError("request_disk_mb must be a strictly positive integer or None")
        for name in ("output_file", "error_file", "log_file"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value or Path(value).is_absolute():
                raise CondorError(f"{name} must be a non-empty relative path")
            if "\n" in value or "\r" in value:
                raise CondorError(f"{name} must not contain a newline")
        if any(not isinstance(item, str) or not item for item in self.arguments):
            raise CondorError("arguments must contain non-empty strings")
        if any(not isinstance(item, str) or not item.strip() for item in self.requirements):
            raise CondorError("requirements must contain non-empty expressions")
        if any(not isinstance(item, str) or not item for item in self.exclude_machines):
            raise CondorError("exclude_machines must contain non-empty machine names")
        for key in self.accounting:
            if not _ACCOUNTING_KEY_RE.fullmatch(key):
                raise CondorError(f"invalid accounting attribute name: {key!r}")
        object.__setattr__(self, "executable", executable)
        object.__setattr__(self, "initialdir", initialdir)
        object.__setattr__(self, "arguments", tuple(self.arguments))
        object.__setattr__(self, "exclude_machines", tuple(self.exclude_machines))
        object.__setattr__(self, "requirements", tuple(self.requirements))
        object.__setattr__(self, "accounting", dict(self.accounting))

    def requirements_expression(self) -> str | None:
        """Combine host exclusions and caller-supplied Condor expressions."""

        clauses = [
            f"(TARGET.Machine != {json.dumps(machine)})" for machine in self.exclude_machines
        ]
        clauses.extend(self.requirements)
        return " && ".join(clauses) if clauses else None

    def render(self) -> str:
        """Render a submit description with deterministic attribute ordering."""

        lines = [
            "# Generated by neo2_for_Er.condor - do not edit in place.",
            f"Executable = {_render_path(self.executable)}",
            f"Initialdir = {_render_path(self.initialdir)}",
            f"Universe = {self.universe}",
            f"Output = {json.dumps(self.output_file)}",
            f"Error = {json.dumps(self.error_file)}",
            f"Log = {json.dumps(self.log_file)}",
            f"Notification = {self.notification}",
            f"request_cpus = {self.request_cpus}",
            f"request_memory = {self.request_memory_mb}",
            f"Getenv = {'true' if self.getenv else 'false'}",
            f"should_transfer_files = {self.should_transfer_files}",
            f"run_as_owner = {'true' if self.run_as_owner else 'false'}",
        ]
        if self.request_disk_mb is not None:
            lines.append(f"request_disk = {self.request_disk_mb}")
        if self.arguments:
            lines.append("Arguments = " + " ".join(json.dumps(arg) for arg in self.arguments))
        requirements = self.requirements_expression()
        if requirements:
            lines.append(f"requirements = {requirements}")
        for key, value in sorted(self.accounting.items()):
            lines.append(f"+{key} = {_render_attribute(value)}")
        lines.append(f"queue {self.queue_count}")
        return "\n".join(lines) + "\n"

    def write(self, directory: str | Path, *, filename: str = "condor.submit") -> Path:
        """Write a new submit file without replacing an existing one."""

        destination = Path(directory)
        if not destination.is_dir():
            raise CondorError(f"submit directory is missing: {destination}")
        if Path(filename).name != filename:
            raise CondorError("submit filename must be a basename")
        path = destination / filename
        if path.exists():
            raise CondorError(f"refusing to overwrite existing submit description: {path}")
        path.write_text(self.render(), encoding="utf-8")
        return path


def _render_path(path: Path) -> str:
    value = str(path)
    return json.dumps(value) if any(character.isspace() for character in value) else value


def _render_attribute(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float) and math.isfinite(value):
        return repr(value)
    raise CondorError(f"unsupported accounting value: {value!r}")


def parse_condor_submit_output(text: str) -> int:
    """Return the cluster identifier acknowledged by ``condor_submit``."""

    match = _CLUSTER_ID_RE.search(text)
    if match is None:
        raise CondorError(f"could not parse a cluster id from condor_submit output: {text!r}")
    return int(match.group("cluster"))


def run_condor(
    arguments: Sequence[str],
    *,
    config: CondorToolConfig,
    tool: str,
) -> subprocess.CompletedProcess[str]:
    """Run one Condor command with a finite timeout and captured output."""

    argv = [str(config.tool(tool)), *(str(argument) for argument in arguments)]
    try:
        return subprocess.run(
            argv,
            check=False,
            capture_output=True,
            text=True,
            timeout=config.command_timeout_s,
        )
    except subprocess.TimeoutExpired as error:
        raise CondorError(f"{tool} exceeded {config.command_timeout_s:g} s") from error
    except OSError as error:
        raise CondorError(f"could not start {tool}: {argv[0]}") from error


def _normalize_job_ad(ad: Mapping[str, object], *, source: str) -> dict[str, object]:
    status = ad.get("JobStatus")
    if isinstance(status, bool) or not isinstance(status, int) or status not in _JOB_STATUS_NAMES:
        raise CondorError(f"unrecognised Condor JobStatus: {status!r}")
    try:
        cluster = _integer_attribute(ad, "ClusterId")
        proc = _integer_attribute(ad, "ProcId")
    except (KeyError, TypeError, ValueError) as error:
        raise CondorError(f"Condor job ad is missing a valid cluster/process id: {ad!r}") from error
    exit_code = ad.get("ExitCode")
    if exit_code is not None:
        try:
            exit_code = _integer_value(exit_code)
        except (TypeError, ValueError) as error:
            raise CondorError(f"Condor job ad has an invalid ExitCode: {exit_code!r}") from error
    return {
        "cluster": cluster,
        "proc": proc,
        "job_status_code": status,
        "job_status": _JOB_STATUS_NAMES[status],
        "exit_code": exit_code,
        "remote_host": ad.get("RemoteHost"),
        "remote_wall_clock_s": ad.get("RemoteWallClockTime"),
        "num_holds": ad.get("NumHolds"),
        "hold_reason": ad.get("HoldReason"),
        "request_cpus": ad.get("RequestCpus"),
        "request_memory_mb": ad.get("RequestMemory"),
        "job_identity": ad.get("KAMELJobIdentity"),
        "ad_source": source,
    }


def _integer_value(value: object) -> int:
    if isinstance(value, bool):
        raise TypeError("boolean is not an integer identifier")
    return int(value)


def _integer_attribute(ad: Mapping[str, object], name: str) -> int:
    return _integer_value(ad[name])


def query_job_ads(
    clusters: Sequence[int],
    *,
    config: CondorToolConfig,
    include_history: bool = True,
) -> tuple[dict[str, object], ...]:
    """Read and normalize queue ads, then history ads, for the selected clusters."""

    if not clusters:
        raise CondorError("no clusters supplied")
    if any(
        isinstance(cluster, bool) or not isinstance(cluster, int) or cluster <= 0
        for cluster in clusters
    ):
        raise CondorError("cluster identifiers must be positive integers")
    records: list[dict[str, object]] = []
    tools = ("condor_q", "condor_history") if include_history else ("condor_q",)
    for tool in tools:
        result = run_condor(
            ["-json", "-attributes", ",".join(_QUERY_ATTRIBUTES), *(str(c) for c in clusters)],
            config=config,
            tool=tool,
        )
        if result.returncode != 0:
            if tool == "condor_history":
                continue
            raise CondorError(
                f"{tool} failed with exit code {result.returncode}: {result.stderr.strip()}"
            )
        for ad in _parse_ads(result.stdout, source=tool):
            records.append(_normalize_job_ad(ad, source=tool))
    return tuple(records)


def _parse_ads(text: str, *, source: str) -> tuple[dict[str, object], ...]:
    if not text.strip():
        return ()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as error:
        raise CondorError(f"{source} returned unparsable JSON") from error
    if isinstance(payload, dict):
        return (payload,)
    if isinstance(payload, list) and all(isinstance(ad, dict) for ad in payload):
        return tuple(payload)
    raise CondorError(f"{source} returned an unexpected JSON structure")


def remove_clusters(clusters: Sequence[int], *, config: CondorToolConfig, reason: str) -> None:
    """Remove selected clusters and require an auditable caller-supplied reason."""

    if not clusters:
        raise CondorError("no clusters supplied for removal")
    if not isinstance(reason, str) or not reason.strip():
        raise CondorError("a removal reason is required")
    result = run_condor([str(cluster) for cluster in clusters], config=config, tool="condor_rm")
    if result.returncode != 0:
        raise CondorError(
            f"condor_rm failed with exit code {result.returncode}: {result.stderr.strip()}"
        )


def verify_shared_filesystem(
    root: str | Path,
    *,
    allowed_prefixes: Sequence[str | Path],
    should_transfer_files: str,
) -> dict[str, object]:
    """Require an explicitly allowed shared root when transfer is disabled."""

    if should_transfer_files not in {"NEVER", "IF_NEEDED", "ALWAYS"}:
        raise CondorError("should_transfer_files must be NEVER, IF_NEEDED, or ALWAYS")
    resolved = Path(root).resolve()
    if should_transfer_files != "NEVER":
        return {
            "shared_filesystem_required": False,
            "should_transfer_files": should_transfer_files,
            "run_root": str(resolved),
            "status": "not_applicable",
        }
    prefixes = tuple(Path(prefix) for prefix in allowed_prefixes)
    if not prefixes:
        raise CondorError("NEVER transfer requires explicit shared_filesystem_prefixes")
    for prefix in prefixes:
        if not prefix.is_absolute():
            raise CondorError(f"shared filesystem prefix must be absolute: {prefix}")
        normalized = prefix.resolve()
        if resolved == normalized or normalized in resolved.parents:
            return {
                "shared_filesystem_required": True,
                "should_transfer_files": "NEVER",
                "run_root": str(resolved),
                "matched_prefix": str(normalized),
                "status": "accepted",
            }
    raise CondorError(
        f"run root {resolved} is outside the declared shared filesystem prefixes "
        f"{[str(prefix) for prefix in prefixes]}; enable transfer or move the run"
    )
