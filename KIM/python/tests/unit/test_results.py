from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest
from kim import Result, ResultError


def complex_compound(values: np.ndarray) -> np.ndarray:
    encoded = np.empty(values.shape, dtype=[("real", "<f8"), ("imag", "<f8")])
    encoded["real"] = values.real
    encoded["imag"] = values.imag
    return encoded


def periodic_file(
    path: Path,
    *,
    radius: np.ndarray | None = None,
    potential: np.ndarray | None = None,
    current: np.ndarray | None = None,
) -> Path:
    radius = np.asarray(radius if radius is not None else [8.0, 9.0, 10.0, 11.0])
    potential = np.asarray(
        potential
        if potential is not None
        else np.arange(1, radius.size + 1) + 1j * np.arange(2, radius.size + 2)
    )
    current = np.asarray(current if current is not None else [1 + 1j] * radius.size)
    with h5py.File(path, "w") as handle:
        radial = handle.create_dataset("backs/e/r", data=radius)
        radial.attrs["unit"] = np.bytes_("cm")
        radial.attrs["comment"] = np.bytes_("Effective radius")
        handle.create_dataset("fields/Phi", data=complex_compound(potential))
        handle.create_dataset("fields/jpar", data=complex_compound(current))
        handle.create_dataset("fields/jpar_e", data=complex_compound(0.25 * current))
        scale = handle.create_dataset("setup/periodic_scale/dx_asis", data=0.6)
        scale.attrs["unit"] = np.bytes_("cm")
    return path


def test_generic_result_discovers_only_datasets_in_sorted_order(tmp_path: Path) -> None:
    path = periodic_file(tmp_path / "result.h5")
    with h5py.File(path, "a") as handle:
        handle.create_group("empty_group")
        handle.create_dataset("runtime", data=2.5)

    assert Result(path).list_datasets() == (
        "backs/e/r",
        "fields/Phi",
        "fields/jpar",
        "fields/jpar_e",
        "runtime",
        "setup/periodic_scale/dx_asis",
    )


def test_reads_are_owned_after_hdf5_file_is_closed(tmp_path: Path) -> None:
    result = Result(periodic_file(tmp_path / "result.h5"))

    radius = result.read_dataset("backs/e/r")
    potential = result.read_dataset("fields/Phi")

    assert isinstance(radius, np.ndarray)
    assert radius.flags.owndata
    np.testing.assert_array_equal(radius, [8.0, 9.0, 10.0, 11.0])
    np.testing.assert_array_equal(potential, [1 + 2j, 2 + 3j, 3 + 4j, 4 + 5j])
    assert potential.dtype == np.complex128


def test_dataset_metadata_decodes_attributes_and_does_not_expose_h5py(tmp_path: Path) -> None:
    metadata = Result(periodic_file(tmp_path / "result.h5")).dataset_metadata("backs/e/r")

    assert metadata.path == "backs/e/r"
    assert metadata.shape == (4,)
    assert metadata.dtype == "float64"
    assert metadata.attributes == {"comment": "Effective radius", "unit": "cm"}


@pytest.mark.parametrize("path", ["", "/fields/Phi", "../secret", "fields/../Phi"])
def test_dataset_paths_are_bounded_relative_paths(tmp_path: Path, path: str) -> None:
    result = Result(periodic_file(tmp_path / "result.h5"))

    with pytest.raises(ResultError, match="dataset path"):
        result.read_dataset(path)


def test_missing_dataset_and_group_reads_have_explicit_errors(tmp_path: Path) -> None:
    result = Result(periodic_file(tmp_path / "result.h5"))

    with pytest.raises(ResultError, match="dataset not found.*fields/missing"):
        result.read_dataset("fields/missing")
    with pytest.raises(ResultError, match="not a dataset.*fields"):
        result.read_dataset("fields")


def test_unreadable_result_file_has_context(tmp_path: Path) -> None:
    path = tmp_path / "broken.h5"
    path.write_text("not HDF5")

    with pytest.raises(ResultError, match=r"could not read KIM result.*broken\.h5"):
        Result(path).list_datasets()


def test_periodic_view_reads_required_and_optional_fields(tmp_path: Path) -> None:
    periodic = Result(periodic_file(tmp_path / "result.h5")).periodic

    np.testing.assert_array_equal(periodic.radius, [8.0, 9.0, 10.0, 11.0])
    np.testing.assert_array_equal(periodic.potential, [1 + 2j, 2 + 3j, 3 + 4j, 4 + 5j])
    np.testing.assert_array_equal(periodic.parallel_current_density, [1 + 1j] * 4)
    assert periodic.as_is_half_width == pytest.approx(0.6)
    assert periodic.resonance_radius == pytest.approx(10.0)
    assert periodic.available_parallel_current_densities == ("jpar", "jpar_e")
    np.testing.assert_array_equal(
        periodic.read_parallel_current_density("jpar_e"), [0.25 + 0.25j] * 4
    )


def test_periodic_integrals_match_campaign_trapezoid_for_both_regions(
    tmp_path: Path,
) -> None:
    radius = np.arange(7.0, 13.0, 0.5)
    current = radius + 1j * (2.0 * radius)
    periodic = Result(
        periodic_file(tmp_path / "result.h5", radius=radius, current=current)
    ).periodic

    as_is = np.abs(radius - 10.0) <= 0.6

    def campaign_integral(r: np.ndarray, jpar: np.ndarray) -> complex:
        integrand = r * jpar
        return complex(2.0 * np.pi * np.sum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(r)))

    assert periodic.integrated_parallel_current("as_is") == pytest.approx(
        campaign_integral(radius[as_is], current[as_is])
    )
    assert periodic.integrated_parallel_current("full_window") == pytest.approx(
        campaign_integral(radius, current)
    )


@pytest.mark.parametrize(
    ("radius", "potential", "message"),
    [
        ([1.0, 2.0, 3.0], [1 + 0j, 2 + 0j], "compatible one-dimensional shapes"),
        ([1.0, 1.0, 2.0], None, "strictly increasing"),
        ([1.0, 2.0, 4.0], None, "equidistant periodic grid"),
    ],
)
def test_periodic_view_rejects_incompatible_shapes_and_grid(
    tmp_path: Path,
    radius: list[float],
    potential: list[complex] | None,
    message: str,
) -> None:
    path = periodic_file(
        tmp_path / "result.h5",
        radius=np.asarray(radius),
        potential=np.asarray(potential) if potential is not None else None,
    )

    with pytest.raises(ResultError, match=message):
        Result(path).periodic


def test_periodic_view_rejects_partial_files(tmp_path: Path) -> None:
    path = tmp_path / "partial.h5"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("fields/Phi", data=complex_compound(np.ones(3)))

    with pytest.raises(ResultError, match="required periodic dataset.*backs/e/r"):
        Result(path).periodic


def test_periodic_current_name_is_constrained_to_discovered_current_fields(
    tmp_path: Path,
) -> None:
    periodic = Result(periodic_file(tmp_path / "result.h5")).periodic

    with pytest.raises(ResultError, match="parallel-current dataset not found"):
        periodic.read_parallel_current_density("Phi")
