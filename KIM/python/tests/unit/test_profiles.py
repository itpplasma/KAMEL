from __future__ import annotations

import hashlib
import runpy
from pathlib import Path

import numpy as np
import pytest
from kim.config import ProfileConfig
from kim.errors import ProfileError
from kim.profiles import ProfileSet

GENERATOR = Path(__file__).parents[1] / "fixtures" / "generate_profiles.py"


def generate(directory: Path) -> dict[str, np.ndarray]:
    return runpy.run_path(str(GENERATOR))["generate_profiles"](directory)


def profile_set(directory: Path) -> ProfileSet:
    return ProfileSet.from_config(ProfileConfig(directory=directory))


def test_generated_profiles_validate_and_have_known_resonant_values(tmp_path: Path) -> None:
    generated = generate(tmp_path / "profiles")

    report = profile_set(tmp_path / "profiles").validate(
        radial_minimum=1.0,
        plasma_radius=9.0,
        m_mode=7,
        n_mode=2,
    )

    assert report.point_count == 11
    assert report.resonance_radius == pytest.approx(5.0)
    resonance_index = np.flatnonzero(generated["radius"] == 5.0).item()
    assert generated["q"][resonance_index] == pytest.approx(-3.5)
    assert generated["n"][resonance_index] == pytest.approx(3.625e13)
    assert generated["Te"][resonance_index] == pytest.approx(3150.0)
    assert generated["Ti"][resonance_index] == pytest.approx(2400.0)


def test_generator_preserves_signed_er_and_vz_profiles(tmp_path: Path) -> None:
    generated = generate(tmp_path / "profiles")

    assert np.all(generated["Er"] < 0.0)
    assert np.all(generated["Vz"] < 0.0)
    assert generated["Er"][-1] == pytest.approx(-0.1)
    assert generated["Vz"][-1] == pytest.approx(-50_000.0)


def test_generator_output_is_byte_deterministic(tmp_path: Path) -> None:
    generate(tmp_path / "first")
    generate(tmp_path / "second")

    for name in ("n", "Te", "Ti", "q", "Er", "Vz"):
        assert (tmp_path / "first" / f"{name}.dat").read_bytes() == (
            tmp_path / "second" / f"{name}.dat"
        ).read_bytes()
    assert (tmp_path / "first" / "n.dat").read_text().splitlines()[0] == (
        "0.0000000000000000e+00 4.5000000000000000e+13"
    )


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("0 1 2\n1 2 3\n", "exactly two numeric columns"),
        ("0 1\n1 nan\n", "row 2.*finite"),
        ("0 1\n0 2\n", "row 2.*strictly increasing"),
        ("0 1\n", "at least two rows"),
    ],
)
def test_invalid_profile_file_reports_file_row_and_constraint(
    tmp_path: Path, contents: str, message: str
) -> None:
    generate(tmp_path / "profiles")
    path = tmp_path / "profiles" / "Te.dat"
    path.write_text(contents)

    with pytest.raises(ProfileError, match=rf"Te\.dat.*{message}"):
        profile_set(tmp_path / "profiles").validate(
            radial_minimum=1.0,
            plasma_radius=9.0,
            m_mode=7,
            n_mode=2,
        )


def test_required_profile_must_exist(tmp_path: Path) -> None:
    generate(tmp_path / "profiles")
    (tmp_path / "profiles" / "q.dat").unlink()

    with pytest.raises(ProfileError, match=r"q\.dat.*required"):
        profile_set(tmp_path / "profiles").validate(
            radial_minimum=1.0,
            plasma_radius=9.0,
            m_mode=7,
            n_mode=2,
        )


def test_core_profiles_require_the_density_grid(tmp_path: Path) -> None:
    generated = generate(tmp_path / "profiles")
    shifted = generated["radius"].copy()
    shifted[4] += 0.01
    np.savetxt(
        tmp_path / "profiles" / "Ti.dat",
        np.column_stack((shifted, generated["Ti"])),
        fmt="%.16e",
    )

    with pytest.raises(ProfileError, match=r"Ti\.dat.*grid.*n\.dat"):
        profile_set(tmp_path / "profiles").validate(
            radial_minimum=1.0,
            plasma_radius=9.0,
            m_mode=7,
            n_mode=2,
        )


def test_profile_roles_must_use_distinct_files(tmp_path: Path) -> None:
    generate(tmp_path / "profiles")
    profiles = ProfileSet.from_config(
        ProfileConfig(directory=tmp_path / "profiles", electron_temperature_file="n.dat")
    )

    with pytest.raises(ProfileError, match=r"n\.dat.*multiple profile roles"):
        profiles.validate(
            radial_minimum=1.0,
            plasma_radius=9.0,
            m_mode=7,
            n_mode=2,
        )


def test_profiles_must_cover_requested_radial_domain(tmp_path: Path) -> None:
    generate(tmp_path / "profiles")

    with pytest.raises(ProfileError, match=r"n\.dat.*cm.*cover.*11"):
        profile_set(tmp_path / "profiles").validate(
            radial_minimum=1.0,
            plasma_radius=11.0,
            m_mode=7,
            n_mode=2,
        )


def test_density_rejects_likely_si_units(tmp_path: Path) -> None:
    generated = generate(tmp_path / "profiles")
    np.savetxt(
        tmp_path / "profiles" / "n.dat",
        np.column_stack((generated["radius"], generated["n"] * 1.0e6)),
        fmt="%.16e",
    )

    with pytest.raises(ProfileError, match=r"n\.dat.*1/cm\^3.*SI"):
        profile_set(tmp_path / "profiles").validate(
            radial_minimum=1.0,
            plasma_radius=9.0,
            m_mode=7,
            n_mode=2,
        )


def test_requested_q_crossing_is_required_inside_domain(tmp_path: Path) -> None:
    generated = generate(tmp_path / "profiles")
    np.savetxt(
        tmp_path / "profiles" / "q.dat",
        np.column_stack((generated["radius"], np.full(11, 2.0))),
        fmt="%.16e",
    )

    with pytest.raises(ProfileError, match=r"q\.dat.*\|m/n\| = 3\.5.*crossing"):
        profile_set(tmp_path / "profiles").validate(
            radial_minimum=1.0,
            plasma_radius=9.0,
            m_mode=7,
            n_mode=2,
        )


def test_q_crossing_between_domain_boundary_and_first_grid_point_is_found(
    tmp_path: Path,
) -> None:
    generated = generate(tmp_path / "profiles")
    q = -(2.0 + generated["radius"])
    np.savetxt(
        tmp_path / "profiles" / "q.dat",
        np.column_stack((generated["radius"], q)),
        fmt="%.16e",
    )

    report = profile_set(tmp_path / "profiles").validate(
        radial_minimum=1.25,
        plasma_radius=9.0,
        m_mode=7,
        n_mode=2,
    )

    assert report.resonance_radius == pytest.approx(1.5)


def test_er_may_use_a_different_interpolation_grid(tmp_path: Path) -> None:
    generated = generate(tmp_path / "profiles")
    radius = np.linspace(0.0, 10.0, 21)
    values = np.interp(radius, generated["radius"], generated["Er"])
    np.savetxt(
        tmp_path / "profiles" / "Er.dat",
        np.column_stack((radius, values)),
        fmt="%.16e",
    )

    profile_set(tmp_path / "profiles").validate(
        radial_minimum=1.0,
        plasma_radius=9.0,
        m_mode=7,
        n_mode=2,
    )


def test_er_and_vz_are_optional(tmp_path: Path) -> None:
    generate(tmp_path / "profiles")
    (tmp_path / "profiles" / "Er.dat").unlink()
    (tmp_path / "profiles" / "Vz.dat").unlink()

    report = profile_set(tmp_path / "profiles").validate(
        radial_minimum=1.0,
        plasma_radius=9.0,
        m_mode=7,
        n_mode=2,
    )

    assert report.resonance_radius == pytest.approx(5.0)


def test_copy_to_stages_regular_files_and_records_sha256(tmp_path: Path) -> None:
    source = tmp_path / "source"
    generate(source)
    profiles = profile_set(source)
    profiles.validate(radial_minimum=1.0, plasma_radius=9.0, m_mode=7, n_mode=2)

    records = profiles.copy_to(tmp_path / "staged")

    assert {record.name for record in records} == {"n", "Te", "Ti", "q", "Er", "Vz"}
    for record in records:
        assert record.destination.is_file()
        assert not record.destination.is_symlink()
        assert record.source.resolve() != record.destination.resolve()
        assert record.sha256 == hashlib.sha256(record.destination.read_bytes()).hexdigest()


def test_copy_preflights_required_files_before_creating_destination(tmp_path: Path) -> None:
    source = tmp_path / "source"
    generate(source)
    (source / "Ti.dat").unlink()
    destination = tmp_path / "staged"

    with pytest.raises(ProfileError, match=r"Ti\.dat.*required"):
        profile_set(source).copy_to(destination)

    assert not destination.exists()
