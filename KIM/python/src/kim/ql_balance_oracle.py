"""Read-only access to the selected QL-Balance input-HDF oracle schema."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from kim.errors import ExperimentalInputError
from numpy.typing import NDArray

_FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class QLBalanceOracle:
    """Detached, validated arrays from a QL-Balance input-HDF oracle.

    Each required source dataset may use either the canonical one-dimensional
    ``(N,)`` representation or a single-row ``(1, N)`` representation.  The
    public arrays are always detached one-dimensional ``(N,)`` arrays.
    The arrays are copied from the source file and made read-only before this
    object is constructed.  Keeping only these arrays makes the result
    independent of both the source file and its HDF5 handle.
    """

    r_out: _FloatArray
    n: _FloatArray
    Te: _FloatArray
    Ti: _FloatArray
    Vz: _FloatArray
    q: _FloatArray
    equilibrium_r: _FloatArray
    equilibrium_psi_pol_norm: _FloatArray
    equilibrium_q: _FloatArray


_DATASET_FIELDS: tuple[tuple[str, str], ...] = (
    ("preprocprof/r_out", "r_out"),
    ("preprocprof/n", "n"),
    ("preprocprof/Te", "Te"),
    ("preprocprof/Ti", "Ti"),
    ("preprocprof/Vz", "Vz"),
    ("preprocprof/q", "q"),
    ("preprocprof/equil/r", "equilibrium_r"),
    ("preprocprof/equil/psi_pol_norm", "equilibrium_psi_pol_norm"),
    ("preprocprof/equil/q", "equilibrium_q"),
)

_PROFILE_FIELDS = ("n", "Te", "Ti", "Vz", "q")
_EQUILIBRIUM_FIELDS = ("equilibrium_psi_pol_norm", "equilibrium_q")


def read_ql_balance_oracle(path: Path | str) -> QLBalanceOracle:
    """Read and validate the fixed QL-Balance input-HDF oracle schema.

    Only the nine required real numeric datasets are read.  Each dataset must
    have either the canonical one-dimensional shape ``(N,)`` or a single-row
    shape ``(1, N)``; the returned arrays always have detached shape ``(N,)``.
    The source is opened explicitly read-only and no HDF5 objects escape this
    function.  All path, HDF5, type, shape, and validation failures are
    reported as :class:`ExperimentalInputError`.
    """

    try:
        candidate = Path(path)
        arrays: dict[str, _FloatArray] = {}
        with h5py.File(candidate, mode="r") as handle:
            for dataset_path, field_name in _DATASET_FIELDS:
                arrays[field_name] = _read_dataset(handle, dataset_path)
        _validate_arrays(arrays)
        return QLBalanceOracle(**arrays)
    except ExperimentalInputError:
        raise
    except Exception as error:
        raise ExperimentalInputError(
            f"{path}: unable to read QL-Balance input-HDF oracle"
        ) from error


def _read_dataset(handle: h5py.File, dataset_path: str) -> _FloatArray:
    """Read one required dataset from the approved ``(N,)``/``(1, N)`` forms."""

    try:
        node = handle[dataset_path]
    except Exception as error:
        raise ExperimentalInputError(f"missing required oracle dataset: {dataset_path}") from error

    if not isinstance(node, h5py.Dataset):
        raise ExperimentalInputError(f"required oracle path is not a dataset: {dataset_path}")
    if len(node.shape) != 1 and not (len(node.shape) == 2 and node.shape[0] == 1):
        raise ExperimentalInputError(
            f"required oracle dataset must have shape (N,) or (1, N): {dataset_path}"
        )
    if not _is_real_numeric(node.dtype):
        raise ExperimentalInputError(
            f"required oracle dataset must contain real numeric values: {dataset_path}"
        )

    try:
        raw = np.asarray(node[...], dtype=np.float64)
        values = raw[0] if raw.ndim == 2 and raw.shape[0] == 1 else raw
        copied_values = values.copy()
    except Exception as error:
        raise ExperimentalInputError(f"unable to read oracle dataset: {dataset_path}") from error
    # A normal copied ndarray owns a mutable buffer and can therefore be made
    # writable again with ``setflags``.  Back the public view with immutable
    # bytes so read-only status cannot be reversed by a caller.
    values = np.frombuffer(copied_values.tobytes(), dtype=np.float64)
    if values.size < 2:
        raise ExperimentalInputError(
            f"required oracle dataset must contain at least two samples: {dataset_path}"
        )
    if not np.all(np.isfinite(values)):
        raise ExperimentalInputError(
            f"required oracle dataset must contain only finite values: {dataset_path}"
        )
    values.setflags(write=False)
    return values


def _is_real_numeric(dtype: np.dtype[object]) -> bool:
    return np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.floating)


def _validate_arrays(arrays: dict[str, _FloatArray]) -> None:
    r_out = arrays["r_out"]
    if not np.all(np.diff(r_out) > 0.0):
        raise ExperimentalInputError("oracle r_out must be strictly increasing")
    profile_length = r_out.size
    for field_name in _PROFILE_FIELDS:
        if arrays[field_name].size != profile_length:
            raise ExperimentalInputError(
                f"oracle profile {field_name} must have the same length as r_out"
            )

    equilibrium_r = arrays["equilibrium_r"]
    if not np.all(np.diff(equilibrium_r) > 0.0):
        raise ExperimentalInputError("oracle equilibrium_r must be strictly increasing")
    equilibrium_length = equilibrium_r.size
    for field_name in _EQUILIBRIUM_FIELDS:
        if arrays[field_name].size != equilibrium_length:
            raise ExperimentalInputError(
                f"oracle equilibrium {field_name} must have the same length as equilibrium_r"
            )

    psi_pol_norm = arrays["equilibrium_psi_pol_norm"]
    if not np.all(np.diff(psi_pol_norm) > 0.0):
        raise ExperimentalInputError("oracle equilibrium_psi_pol_norm must be strictly increasing")
    if not np.all((psi_pol_norm >= 0.0) & (psi_pol_norm <= 1.0)):
        raise ExperimentalInputError("oracle equilibrium_psi_pol_norm must be within [0, 1]")


__all__ = ["QLBalanceOracle", "read_ql_balance_oracle"]
