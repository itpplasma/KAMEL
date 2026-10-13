"""Bounded local NEO-2 execution using KAMEL's existing real-space template."""

import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import f90nml
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


def stage_surfaces(
    workdir: str | Path, template_dir: str | Path, surfaces_per_run: int = 1
) -> tuple[Path, ...]:
    """Stage explicit seven-column surfaces, optionally in native PAR chunks.

    The caller chooses surfaces. Existing directories are never overwritten.
    ``surfaces_per_run=1`` preserves scalar runs; zero runs all surfaces together.
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
    if not isinstance(surfaces_per_run, int) or surfaces_per_run < 0:
        raise Neo2LocalError("surfaces_per_run must be a nonnegative integer")
    batch = surfaces_per_run != 1
    chunk_size = surfaces_per_run or len(rows)
    table = root / "k_surfaces.dat"
    if batch:
        if table.exists() or (root / "jobs_list.txt").exists():
            raise Neo2LocalError("refusing to overwrite an existing staged batch")
        for first in range(0, len(rows), chunk_size):
            if (root / f"chunk_{first + 1:06d}_{min(first + chunk_size, len(rows)):06d}").exists():
                raise Neo2LocalError("refusing to overwrite an existing chunk directory")
        np.savetxt(table, rows[:, [0, 2, 3, 6]], fmt="%.9e")
    jobs = []
    for first in range(0, len(rows), chunk_size):
        last = min(first + chunk_size, len(rows))
        s, _r, rbeg, zbeg, _ti, _ne, kappa = rows[first]
        name = f"chunk_{first + 1:06d}_{last:06d}" if batch else "s" + _fmt(s).replace("d-", "m")
        job = root / name
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
        if batch:
            nml = f90nml.read(job / "neo2.in")
            nml["settings"].update(
                k_surface_file="../k_surfaces.dat",
                k_output_file="k_coefficients.h5",
                surface_first=first + 1,
                surface_last=last,
            )
            nml["collision"]["k_coefficient_only"] = True
            nml.write(job / "neo2.in", force=True)
        jobs.append(job)
    (root / "jobs_list.txt").write_text("".join(job.name + "\n" for job in jobs))
    return tuple(jobs)


def run_staged_surfaces(workdir: str | Path, plan: LocalNeo2Plan) -> np.ndarray:
    """Run staged jobs locally and return sorted ``r_eff_cm, k_cof`` values."""
    root = Path(workdir)
    rows = np.loadtxt(root / "surfaces.dat", ndmin=2)
    jobs = [root / line for line in (root / "jobs_list.txt").read_text().splitlines()]
    if (not (root / "k_surfaces.dat").is_file() and len(jobs) != len(rows)) or any(
        not (job / "neo2.in").is_file() for job in jobs
    ):
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
    return collect_staged_surfaces(root)


def collect_staged_surfaces(workdir: str | Path) -> np.ndarray:
    """Collect scalar or native chunk results using explicit table row ranges.

    Every native ``integral_converged`` flag must be one to return a profile.
    """
    root = Path(workdir)
    rows = np.loadtxt(root / "surfaces.dat", ndmin=2)
    jobs = [root / line for line in (root / "jobs_list.txt").read_text().splitlines()]
    if (root / "k_surfaces.dat").is_file():
        native = np.loadtxt(root / "k_surfaces.dat", ndmin=2)
        expected = np.array(
            [[float(f"{value:.9e}") for value in row] for row in rows[:, [0, 2, 3, 6]]]
        )
        if not np.array_equal(native, expected):
            raise Neo2LocalError("native input table does not match surfaces.dat")
        covered = np.zeros(len(rows), dtype=bool)
        result = np.column_stack((rows[:, 1], np.empty(len(rows))))
        for job in jobs:
            try:
                with h5py.File(job / "k_coefficients.h5") as output:
                    bounds = np.array(
                        [output[key][()].item() for key in ("surface_first", "surface_last")]
                    )
                    first, last = bounds.astype(int)
                    if (
                        not np.array_equal(bounds, [first, last])
                        or not 1 <= first <= last <= len(rows)
                        or covered[first - 1 : last].any()
                    ):
                        raise Neo2LocalError(f"invalid or overlapping NEO-2 row range in {job}")
                    keys = (
                        "boozer_s",
                        "rbeg",
                        "zbeg",
                        "conl_over_mfp",
                        "k_cof",
                        "qcurr_g",
                        "qcurr_t",
                        "integral_converged",
                    )
                    values = [np.asarray(output[key][()]) for key in keys]
                    if (
                        output["k_coefficient_only"][()].item() != 1
                        or any(
                            value.shape != (last - first + 1,) or not np.isfinite(value).all()
                            for value in values
                        )
                        or not np.array_equal(np.column_stack(values[:4]), native[first - 1 : last])
                        or not (values[-1] == 1).all()
                    ):
                        raise Neo2LocalError(f"invalid NEO-2 batch output in {job}")
                    result[first - 1 : last, 1] = values[4]
                    covered[first - 1 : last] = True
            except (OSError, KeyError, TypeError, ValueError) as error:
                raise Neo2LocalError(f"invalid NEO-2 batch output in {job}: {error}") from error
        if not covered.all():
            raise Neo2LocalError("NEO-2 batches do not cover all surfaces.dat rows")
        return result[np.argsort(result[:, 0], kind="stable")]
    if len(jobs) != len(rows):
        raise Neo2LocalError("staged jobs do not match surfaces.dat")
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
