from __future__ import annotations

import runpy
import sys
from pathlib import Path

import f90nml
import numpy as np
import pytest
from kim import BuiltinPlasma, PlasmaIsotope, ProfileScale, SimulationConfig, SweepSpec
from kim.condor import KimCondorPlan, stage_condor_sweep
from kim.condor_client import CondorError, CondorToolConfig
from kim.namelist import dumps_namelist
from pydantic import ValidationError

FIXTURES = Path(__file__).parents[1] / "fixtures"


def _base_config(profile_directory: Path) -> SimulationConfig:
    runpy.run_path(str(FIXTURES / "generate_profiles.py"))["generate_profiles"](profile_directory)
    return SimulationConfig.electrostatic_periodic(
        profiles=profile_directory,
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-17_977.413,
        major_radius=165.0,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=1.0,
        plasma_radius=9.0,
    )


def _executable(path: Path) -> Path:
    path.write_text(f"#!{sys.executable}\npass\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _plan(
    *,
    backend: str,
    executable: Path,
    shared_prefix: Path,
    base_namelist: Path | None = None,
    **values,
) -> KimCondorPlan:
    return KimCondorPlan(
        backend=backend,
        executable=executable,
        python_executable=Path(sys.executable),
        condor=CondorToolConfig(shared_prefix / "condor-bin"),
        shared_filesystem_prefixes=(shared_prefix,),
        should_transfer_files="NEVER",
        base_namelist=base_namelist,
        **values,
    )


def _spec(config: SimulationConfig, factors: tuple[float, ...] = (-1.0, 0.0, 1.0)) -> SweepSpec:
    return SweepSpec(base=config, variation=ProfileScale(values=factors))


def _profile_bytes(directory: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in directory.glob("*.dat")}


def test_plan_requires_explicit_backend_and_valid_resources(tmp_path: Path) -> None:
    executable = _executable(tmp_path / "KIM.x")

    with pytest.raises(ValidationError, match="backend"):
        KimCondorPlan(
            executable=executable,
            python_executable=Path(sys.executable),
        )

    with pytest.raises(ValidationError, match="request_cpus"):
        KimCondorPlan(
            backend="kim_x_namelist",
            executable=executable,
            python_executable=Path(sys.executable),
            request_cpus=0,
        )


def test_stage_creates_one_job_in_scan_order(tmp_path: Path) -> None:
    source = tmp_path / "source-profiles"
    config = _base_config(source)
    executable = _executable(tmp_path / "KIM.x")
    root = tmp_path / "run"
    plan = _plan(
        backend="kamel_kim_python",
        executable=executable,
        shared_prefix=tmp_path,
    )

    jobs = stage_condor_sweep(root, spec=_spec(config), plan=plan)

    assert [job.profile_scale_factor for job in jobs] == [-1.0, 0.0, 1.0]
    assert [job.scan_order_index for job in jobs] == [0, 1, 2]
    assert [job.job for job in jobs] == ["scale-0000", "scale-0001", "scale-0002"]
    assert [job.directory for job in jobs] == [root / job.job for job in jobs]
    for job in jobs:
        payload = __import__("json").loads(
            (job.directory / "condor_job.json").read_text(encoding="utf-8")
        )
        assert payload["backend"] == "kamel_kim_python"
        assert payload["scan_order_index"] == job.scan_order_index
        assert payload["profile_scale_factor"] == job.profile_scale_factor
        assert payload["base_config"]["profiles"]["directory"] == str(job.directory / "profiles")
        assert payload["mode"] == "kim-sweep"
        assert payload["job_timeout_s"] == 3600.0
        assert payload["poll_interval_s"] == 30.0
        assert payload["max_wall_clock_s"] == 21600.0


def test_api_backend_stages_unscaled_source_profiles(tmp_path: Path) -> None:
    source = tmp_path / "source-profiles"
    config = _base_config(source)
    source_bytes = _profile_bytes(source)
    executable = _executable(tmp_path / "KIM.x")
    root = tmp_path / "api-run"
    plan = _plan(
        backend="kamel_kim_python",
        executable=executable,
        shared_prefix=tmp_path,
    )

    jobs = stage_condor_sweep(root, spec=_spec(config), plan=plan)

    source_er = np.loadtxt(source / "Er.dat")
    for job in jobs:
        staged = job.directory / "profiles"
        np.testing.assert_array_equal(np.loadtxt(staged / "Er.dat"), source_er)
        for name, original in source_bytes.items():
            assert (staged / name).read_bytes() == original
    assert _profile_bytes(source) == source_bytes


def test_namelist_backend_scales_only_er_copy(tmp_path: Path) -> None:
    source = tmp_path / "source-profiles"
    config = _base_config(source)
    source_bytes = _profile_bytes(source)
    reviewed = tmp_path / "reviewed.nml"
    reviewed.write_text(dumps_namelist(config), encoding="utf-8")
    reviewed_bytes = reviewed.read_bytes()
    executable = _executable(tmp_path / "KIM.x")
    root = tmp_path / "namelist-run"
    plan = _plan(
        backend="kim_x_namelist",
        executable=executable,
        shared_prefix=tmp_path,
        base_namelist=reviewed,
    )

    jobs = stage_condor_sweep(root, spec=_spec(config), plan=plan)

    source_er = np.loadtxt(source / "Er.dat")
    for job in jobs:
        staged = job.directory / "profiles"
        factor = job.profile_scale_factor
        staged_er = np.loadtxt(staged / "Er.dat")
        np.testing.assert_array_equal(staged_er[:, 0], source_er[:, 0])
        np.testing.assert_allclose(staged_er[:, 1], source_er[:, 1] * factor)
        for name, original in source_bytes.items():
            if name != "Er.dat":
                assert (staged / name).read_bytes() == original
        staged_namelist = f90nml.read(job.directory / "KIM_config.nml").todict()
        assert staged_namelist["kim_profiles"]["input_profile_dir"] == str(staged) + "/"
        payload = __import__("json").loads(
            (job.directory / "condor_job.json").read_text(encoding="utf-8")
        )
        assert payload["mode"] == "kim-run"
        assert payload["namelist_file"] == "KIM_config.nml"
        assert "KIM_config.nml" in payload["input_sha256"]
    assert _profile_bytes(source) == source_bytes
    assert reviewed.read_bytes() == reviewed_bytes


def test_namelist_overrides_must_exist_in_reviewed_base(tmp_path: Path) -> None:
    source = tmp_path / "source-profiles"
    config = _base_config(source)
    reviewed = tmp_path / "reviewed.nml"
    reviewed.write_text(dumps_namelist(config), encoding="utf-8")
    executable = _executable(tmp_path / "KIM.x")
    plan = _plan(
        backend="kim_x_namelist",
        executable=executable,
        shared_prefix=tmp_path,
        base_namelist=reviewed,
        namelist_overrides={"kim_grid.unreviewed_tuning": 12},
    )

    with pytest.raises(CondorError, match="not present in reviewed base namelist"):
        stage_condor_sweep(tmp_path / "run", spec=_spec(config, (1.0,)), plan=plan)

    assert not (tmp_path / "run").exists()


def test_stage_refuses_existing_job_directory(tmp_path: Path) -> None:
    source = tmp_path / "source-profiles"
    config = _base_config(source)
    executable = _executable(tmp_path / "KIM.x")
    root = tmp_path / "run"
    root.mkdir()
    existing = root / "scale-0001"
    existing.mkdir()
    (existing / "user-data.txt").write_text("keep", encoding="utf-8")
    plan = _plan(
        backend="kamel_kim_python",
        executable=executable,
        shared_prefix=tmp_path,
    )

    with pytest.raises(CondorError, match="refusing to overwrite existing job directory"):
        stage_condor_sweep(root, spec=_spec(config), plan=plan)

    assert list(root.iterdir()) == [existing]
    assert (existing / "user-data.txt").read_text(encoding="utf-8") == "keep"
