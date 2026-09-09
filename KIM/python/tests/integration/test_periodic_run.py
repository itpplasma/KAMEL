from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pytest
from kim import (
    BuiltinPlasma,
    ExecutableError,
    GridConfig,
    PeriodicConfig,
    PlasmaIsotope,
    Result,
    ResultError,
    RunStatus,
    Simulation,
    SimulationConfig,
    resolve_executable,
)

REFERENCE = Path(__file__).parents[1] / "fixtures" / "periodic_reference.json"


def _generate_profiles(directory: Path) -> None:
    """Write a compact, physical-scale parabolic r_eff test case in CGS units."""

    directory.mkdir(parents=True)
    radius = np.linspace(3.0, 67.0, 65)
    normalized = radius / 67.0
    parabola = 1.0 - normalized**2
    profiles = {
        "n": 5.0e11 + 4.6e13 * parabola,
        "Te": 50.0 + 3_800.0 * parabola,
        "Ti": 70.0 + 3_350.0 * parabola,
        "q": -(0.99 + 4.55 * normalized**2),
        "Er": 0.08 - 0.16 * normalized**2,
        "Vz": 1_500.0 + 1.15e7 * parabola,
    }
    for name, values in profiles.items():
        np.savetxt(
            directory / f"{name}.dat",
            np.column_stack((radius, values)),
            fmt="%.16e",
        )


def _configuration(profiles: Path) -> SimulationConfig:
    return SimulationConfig.electrostatic_periodic(
        profiles=profiles,
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-17_977.413,
        major_radius=165.0,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=3.0,
        plasma_radius=63.0,
        grid=GridConfig(
            radial_minimum=3.0,
            plasma_radius=63.0,
            l_space_dim=64,
            rg_space_dim=64,
            larmor_skip_factor=20.0,
            gauss_nodes_x=11,
            gauss_nodes_x_prime=10,
            gauss_nodes_theta=7,
        ),
        periodic=PeriodicConfig(
            as_is_width_scale=2.0,
            transition_width_scale=4.0,
            wavenumber_cutoff_scale=2.0,
            n_rg=64,
        ),
    )


def _reference() -> dict[str, Any]:
    return json.loads(REFERENCE.read_text(encoding="utf-8"))


def _integration_executable() -> Path:
    enabled = os.environ.get("KIM_RUN_INTEGRATION", "").casefold() in {"1", "true", "yes"}
    if not enabled:
        pytest.skip("set KIM_RUN_INTEGRATION=1 to run the real KIM integration reference")
    try:
        return resolve_executable()
    except ExecutableError as error:
        pytest.skip(f"KIM.x is unavailable: {error}")


def _assert_reference(path: Path, reference: dict[str, Any]) -> None:
    generic = Result(path)
    datasets = set(generic.list_datasets())
    missing = set(reference["required_datasets"]) - datasets
    assert not missing, f"missing mandatory datasets: {sorted(missing)}"

    periodic = generic.periodic
    expected = reference["expected"]
    tolerances = reference["tolerances"]
    assert list(periodic.radius.shape) == expected["field_shape"]
    assert np.all(np.isfinite(periodic.potential))
    assert np.all(np.isfinite(periodic.parallel_current_density))
    assert np.max(np.abs(periodic.potential)) > 0.0
    assert np.max(np.abs(periodic.parallel_current_density)) > 0.0
    assert periodic.resonance_radius == pytest.approx(
        expected["resonance_radius_cm"], abs=tolerances["resonance_absolute_cm"]
    )

    with h5py.File(path, "r") as handle:
        assert int(handle["setup/periodic_scale/M"][()]) == expected["periodic_modes"]
        assert int(handle["setup/periodic_scale/N_rg"][()]) == expected["periodic_grid_points"]

    for region, key in (
        ("as_is", "integrated_parallel_current_as_is_statA"),
        ("full_window", "integrated_parallel_current_full_window_statA"),
    ):
        target = complex(expected[key]["real"], expected[key]["imag"])
        actual = periodic.integrated_parallel_current(region)
        assert actual == pytest.approx(
            target,
            rel=tolerances["integral_relative"],
            abs=tolerances["integral_absolute_statA"],
        )


def test_reference_check_rejects_missing_mandatory_dataset(tmp_path: Path) -> None:
    partial = tmp_path / "partial.h5"
    with h5py.File(partial, "w") as handle:
        handle.create_dataset("backs/e/r", data=np.arange(4.0))

    with pytest.raises(AssertionError, match="missing mandatory datasets"):
        _assert_reference(partial, _reference())


def test_periodic_parabolic_reference_against_kim(tmp_path: Path) -> None:
    executable = _integration_executable()
    profiles = tmp_path / "profiles"
    _generate_profiles(profiles)
    q_profile = np.loadtxt(profiles / "q.dat")
    q_crossing = np.interp(3.5, np.abs(q_profile[:, 1]), q_profile[:, 0])
    reference = _reference()
    assert q_crossing == pytest.approx(
        reference["expected"]["profile_linear_q_crossing_cm"],
        abs=reference["tolerances"]["resonance_absolute_cm"],
    )

    run = Simulation(
        _configuration(profiles),
        executable=executable,
        runs_directory=tmp_path / "runs",
        timeout=60.0,
    ).run(environment={"OMP_NUM_THREADS": "1"})

    assert run.status is RunStatus.SUCCEEDED, run.manifest.failure
    try:
        _assert_reference(run.output_file, reference)
    except ResultError as error:
        pytest.fail(f"periodic output does not satisfy the typed result contract: {error}")
