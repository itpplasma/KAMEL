"""Synthetic, read-only contract tests for the QL-Balance HDF5 oracle.

The fixture is deliberately created at test time: a real QL-Balance file is an
oracle only and must never become an input fixture or be modified by the
reader.  Profile arrays are required to share ``r_out``'s length, while the
three equilibrium arrays are required to share ``equilibrium_r``'s length.
"""

from __future__ import annotations

import hashlib
from dataclasses import FrozenInstanceError
from pathlib import Path

import h5py
import numpy as np
import pytest
from kim import ExperimentalInputError
from kim.ql_balance_oracle import QLBalanceOracle, read_ql_balance_oracle

_DATASETS = (
    "preprocprof/r_out",
    "preprocprof/n",
    "preprocprof/Te",
    "preprocprof/Ti",
    "preprocprof/Vz",
    "preprocprof/q",
    "preprocprof/equil/r",
    "preprocprof/equil/psi_pol_norm",
    "preprocprof/equil/q",
)


def _valid_data() -> dict[str, np.ndarray]:
    """Return a three-point fixture with signed q values and bounded psi."""

    return {
        "preprocprof/r_out": np.array([10.0, 20.0, 30.0], dtype=np.float64),
        "preprocprof/n": np.array([1.0e13, 2.0e13, 3.0e13], dtype=np.float64),
        "preprocprof/Te": np.array([100.0, 80.0, 60.0], dtype=np.float64),
        "preprocprof/Ti": np.array([90.0, 70.0, 50.0], dtype=np.float64),
        "preprocprof/Vz": np.array([-2.0e5, 0.0, 2.0e5], dtype=np.float64),
        "preprocprof/q": np.array([-1.25, 0.5, 2.0], dtype=np.float64),
        "preprocprof/equil/r": np.array([5.0, 25.0, 45.0], dtype=np.float64),
        "preprocprof/equil/psi_pol_norm": np.array([0.0, 0.25, 1.0], dtype=np.float64),
        "preprocprof/equil/q": np.array([1.5, -0.75, 0.25], dtype=np.float64),
    }


def _write_oracle(
    path: Path,
    *,
    overrides: dict[str, np.ndarray] | None = None,
) -> None:
    values = _valid_data()
    if overrides:
        values.update(overrides)
    with h5py.File(path, "w") as handle:
        for dataset_path, data in values.items():
            handle.create_dataset(dataset_path, data=data)


def test_reads_exact_ql_balance_dataset_mapping_and_preserves_signed_q(
    tmp_path: Path,
) -> None:
    path = tmp_path / "ql-balance-oracle.h5"
    _write_oracle(path)

    result = read_ql_balance_oracle(path)

    assert isinstance(result, QLBalanceOracle)
    np.testing.assert_array_equal(result.r_out, [10.0, 20.0, 30.0])
    np.testing.assert_array_equal(result.n, [1.0e13, 2.0e13, 3.0e13])
    np.testing.assert_array_equal(result.Te, [100.0, 80.0, 60.0])
    np.testing.assert_array_equal(result.Ti, [90.0, 70.0, 50.0])
    np.testing.assert_array_equal(result.Vz, [-2.0e5, 0.0, 2.0e5])
    np.testing.assert_array_equal(result.q, [-1.25, 0.5, 2.0])
    np.testing.assert_array_equal(result.equilibrium_r, [5.0, 25.0, 45.0])
    np.testing.assert_array_equal(result.equilibrium_psi_pol_norm, [0.0, 0.25, 1.0])
    np.testing.assert_array_equal(result.equilibrium_q, [1.5, -0.75, 0.25])


def test_rejects_missing_oracle_file(tmp_path: Path) -> None:
    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(tmp_path / "missing.h5")


@pytest.mark.parametrize("missing_dataset", _DATASETS)
def test_rejects_missing_required_oracle_dataset(
    tmp_path: Path,
    missing_dataset: str,
) -> None:
    path = tmp_path / "missing-dataset.h5"
    _write_oracle(path)
    with h5py.File(path, "a") as handle:
        del handle[missing_dataset]

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


@pytest.mark.parametrize(
    ("dataset_path", "replacement"),
    [
        ("preprocprof/r_out", np.array([[10.0, 20.0, 30.0]], dtype=np.float64)),
        ("preprocprof/n", np.array([[1.0e13], [2.0e13], [3.0e13]], dtype=np.float64)),
        (
            "preprocprof/equil/r",
            np.array([[5.0, 25.0, 45.0]], dtype=np.float64),
        ),
        (
            "preprocprof/equil/psi_pol_norm",
            np.array([[0.0], [0.25], [1.0]], dtype=np.float64),
        ),
    ],
)
def test_rejects_non_one_dimensional_oracle_dataset(
    tmp_path: Path,
    dataset_path: str,
    replacement: np.ndarray,
) -> None:
    path = tmp_path / "non-one-dimensional.h5"
    _write_oracle(path, overrides={dataset_path: replacement})

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


@pytest.mark.parametrize(
    ("dataset_path", "replacement"),
    [
        ("preprocprof/n", np.array([1.0e13, 2.0e13], dtype=np.float64)),
        ("preprocprof/Te", np.array([100.0, 80.0], dtype=np.float64)),
        ("preprocprof/Ti", np.array([90.0, 70.0], dtype=np.float64)),
        ("preprocprof/Vz", np.array([-2.0e5, 0.0], dtype=np.float64)),
        ("preprocprof/q", np.array([-1.25, 0.5], dtype=np.float64)),
    ],
)
def test_requires_each_profile_to_share_r_out_length(
    tmp_path: Path,
    dataset_path: str,
    replacement: np.ndarray,
) -> None:
    path = tmp_path / "profile-shape-mismatch.h5"
    _write_oracle(path, overrides={dataset_path: replacement})

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


@pytest.mark.parametrize(
    ("dataset_path", "replacement"),
    [
        ("preprocprof/equil/psi_pol_norm", np.array([0.0, 1.0], dtype=np.float64)),
        ("preprocprof/equil/q", np.array([1.5, -0.75], dtype=np.float64)),
    ],
)
def test_requires_equilibrium_arrays_to_share_equilibrium_r_length(
    tmp_path: Path,
    dataset_path: str,
    replacement: np.ndarray,
) -> None:
    path = tmp_path / "equilibrium-shape-mismatch.h5"
    _write_oracle(path, overrides={dataset_path: replacement})

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


@pytest.mark.parametrize(
    ("dataset_path", "replacement"),
    [
        ("preprocprof/n", np.array([1.0e13, np.nan, 3.0e13], dtype=np.float64)),
        ("preprocprof/q", np.array([-1.25, np.inf, 2.0], dtype=np.float64)),
        (
            "preprocprof/equil/psi_pol_norm",
            np.array([0.0, np.nan, 1.0], dtype=np.float64),
        ),
    ],
)
def test_rejects_non_finite_oracle_values(
    tmp_path: Path,
    dataset_path: str,
    replacement: np.ndarray,
) -> None:
    path = tmp_path / "non-finite.h5"
    _write_oracle(path, overrides={dataset_path: replacement})

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


@pytest.mark.parametrize("dataset_path", ["preprocprof/r_out", "preprocprof/equil/r"])
def test_requires_strictly_increasing_radial_coordinates(
    tmp_path: Path,
    dataset_path: str,
) -> None:
    path = tmp_path / "non-increasing-radius.h5"
    _write_oracle(
        path,
        overrides={dataset_path: np.array([10.0, 10.0, 30.0], dtype=np.float64)},
    )

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


@pytest.mark.parametrize(
    "psi_pol_norm",
    [
        np.array([0.0, 0.0, 1.0], dtype=np.float64),
        np.array([0.0, 1.1, 1.2], dtype=np.float64),
        np.array([-0.1, 0.25, 1.0], dtype=np.float64),
    ],
)
def test_requires_strictly_increasing_bounded_normalized_poloidal_flux(
    tmp_path: Path,
    psi_pol_norm: np.ndarray,
) -> None:
    path = tmp_path / "invalid-psi-pol-norm.h5"
    _write_oracle(
        path,
        overrides={"preprocprof/equil/psi_pol_norm": psi_pol_norm},
    )

    with pytest.raises(ExperimentalInputError):
        read_ql_balance_oracle(path)


def test_oracle_is_detached_and_deeply_read_only_after_hdf5_close(tmp_path: Path) -> None:
    path = tmp_path / "detached.h5"
    _write_oracle(path)

    result = read_ql_balance_oracle(path)
    arrays = (
        result.r_out,
        result.n,
        result.Te,
        result.Ti,
        result.Vz,
        result.q,
        result.equilibrium_r,
        result.equilibrium_psi_pol_norm,
        result.equilibrium_q,
    )

    for array in arrays:
        assert isinstance(array, np.ndarray)
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array[0] = array[0]
        with pytest.raises(ValueError):
            array.setflags(write=True)

    with pytest.raises((FrozenInstanceError, AttributeError, TypeError)):
        result.r_out = np.array([1.0])  # type: ignore[misc]


def test_reader_opens_source_read_only_and_preserves_bytes_hash_and_mtime(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "read-only-source.h5"
    _write_oracle(path)
    original_bytes = path.read_bytes()
    original_hash = hashlib.sha256(original_bytes).hexdigest()
    original_mtime_ns = path.stat().st_mtime_ns

    real_file = h5py.File
    opened_modes: list[str] = []

    def spy_file(name: object, *args: object, **kwargs: object) -> h5py.File:
        mode = kwargs.get("mode")
        if mode is None and len(args) >= 1:
            mode = args[0]
        opened_modes.append("r" if mode is None else str(mode))
        return real_file(name, *args, **kwargs)

    monkeypatch.setattr(h5py, "File", spy_file)
    read_ql_balance_oracle(path)

    assert opened_modes == ["r"]
    assert path.read_bytes() == original_bytes
    assert hashlib.sha256(path.read_bytes()).hexdigest() == original_hash
    assert path.stat().st_mtime_ns == original_mtime_ns
