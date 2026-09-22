"""Bounded local NEO-2 execution using KAMEL's existing real-space template."""

import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


class Neo2LocalError(RuntimeError):
    """A staged local NEO-2 run is invalid or failed."""


@dataclass(frozen=True)
class LocalNeo2Plan:
    executable: Path
    physical_cores: int = 28
    omp_threads: int = 1
    max_workers: int | None = None

    def __post_init__(self):
        exe = Path(self.executable)
        if not exe.is_absolute() or not exe.is_file():
            raise Neo2LocalError("executable must be an existing absolute path")
        workers = self.max_workers or self.physical_cores // self.omp_threads
        if (
            min(self.physical_cores, self.omp_threads, workers) <= 0
            or workers * self.omp_threads > self.physical_cores
        ):
            raise Neo2LocalError("invalid or oversubscribed NEO-2 resource plan")
        object.__setattr__(self, "executable", exe)
        object.__setattr__(self, "max_workers", workers)


def _fmt(value: float) -> str:
    return f"{value:.9e}".replace("e", "d")


def stage_surfaces(workdir: str | Path, template_dir: str | Path) -> tuple[Path, ...]:
    """Copy a template into one directory per explicit seven-column surface.

    The caller chooses surfaces. Existing directories are never overwritten.
    """
    root, template = Path(workdir), Path(template_dir)
    rows = np.loadtxt(root / "surfaces.dat", ndmin=2)
    if (
        rows.ndim != 2
        or rows.shape[1] != 7
        or not np.all(np.isfinite(rows))
        or not (template / "neo2.in").is_file()
    ):
        raise Neo2LocalError("need finite seven-column surfaces.dat and template neo2.in")
    jobs = []
    for s, _r, rbeg, zbeg, _ti, _ne, kappa in rows:
        job = root / ("s" + _fmt(s).replace("d-", "m"))
        if job.exists():
            raise Neo2LocalError(f"refusing to overwrite {job}")
        shutil.copytree(template, job, symlinks=True)
        text = (job / "neo2.in").read_text()
        for token, value in {
            "<boozer_s>": s,
            "<conl_over_mfp>": kappa,
            "<rbeg>": rbeg,
            "<zbeg>": zbeg,
        }.items():
            if token not in text:
                raise Neo2LocalError(f"missing template token {token}")
            text = text.replace(token, _fmt(value))
        (job / "neo2.in").write_text(text)
        jobs.append(job)
    (root / "jobs_list.txt").write_text("".join(job.name + "\n" for job in jobs))
    return tuple(jobs)


def run_staged_surfaces(workdir: str | Path, plan: LocalNeo2Plan) -> np.ndarray:
    """Run staged jobs locally and return sorted ``r_eff_cm, k_cof`` values."""
    root = Path(workdir)
    rows = np.loadtxt(root / "surfaces.dat", ndmin=2)
    jobs = [root / line for line in (root / "jobs_list.txt").read_text().splitlines()]
    if len(jobs) != len(rows) or any(not (job / "neo2.in").is_file() for job in jobs):
        raise Neo2LocalError("staged jobs do not match surfaces.dat")

    def run(job):
        env = os.environ | {"OMP_NUM_THREADS": str(plan.omp_threads)}
        result = subprocess.run(
            [str(plan.executable)], cwd=job, capture_output=True, text=True, env=env
        )
        (job / "neo2.stdout.log").write_text(result.stdout)
        (job / "neo2.stderr.log").write_text(result.stderr)
        return job, result.returncode

    with ThreadPoolExecutor(max_workers=plan.max_workers) as pool:
        failed = [
            job.name
            for job, code in (
                future.result() for future in as_completed([pool.submit(run, job) for job in jobs])
            )
            if code
        ]
    if failed:
        raise Neo2LocalError("failed NEO-2 jobs: " + ", ".join(sorted(failed)))
    result = []
    for row, job in zip(rows, jobs, strict=True):
        with h5py.File(job / "neo2_config.h5") as config, h5py.File(
            job / "fulltransp.h5"
        ) as transport:
            if not np.isclose(float(config["settings/boozer_s"][()]), row[0]) or not np.isfinite(
                float(transport["k_cof"][()])
            ):
                raise Neo2LocalError(f"invalid NEO-2 output in {job}")
            result.append((row[1], float(transport["k_cof"][()])))
    return np.asarray(sorted(result))
