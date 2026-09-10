"""Safe generic and typed access to KIM HDF5 results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Literal

import h5py
import numpy as np
from kim.errors import ResultError
from numpy.typing import NDArray

IntegrationRegion = Literal["as_is", "full_window"]


@dataclass(frozen=True)
class DatasetMetadata:
    """Serializable structural metadata for one HDF5 dataset."""

    path: str
    shape: tuple[int, ...]
    dtype: str
    attributes: dict[str, Any]


class Result:
    """Bounded, file-lifetime-safe access to one KIM HDF5 output file."""

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path).expanduser().absolute()

    def list_datasets(self) -> tuple[str, ...]:
        """Return every dataset path in deterministic order."""

        datasets: list[str] = []
        with self._open() as handle:

            def collect(name: str, item: h5py.Group | h5py.Dataset) -> None:
                if isinstance(item, h5py.Dataset):
                    datasets.append(name)

            handle.visititems(collect)
        return tuple(sorted(datasets))

    def read_dataset(self, path: str) -> Any:
        """Read and detach a dataset, decoding KIM compound complex values."""

        normalized = _dataset_path(path)
        with self._open() as handle:
            item = handle.get(normalized)
            if item is None:
                raise ResultError(f"dataset not found in {self.path}: {normalized}")
            if not isinstance(item, h5py.Dataset):
                raise ResultError(f"HDF5 object is not a dataset: {normalized}")
            value = _decode_value(item[()])
        if isinstance(value, np.ndarray):
            value.setflags(write=False)
        return value

    def dataset_metadata(self, path: str) -> DatasetMetadata:
        """Return detached shape, type, and attribute metadata."""

        normalized = _dataset_path(path)
        with self._open() as handle:
            item = handle.get(normalized)
            if item is None:
                raise ResultError(f"dataset not found in {self.path}: {normalized}")
            if not isinstance(item, h5py.Dataset):
                raise ResultError(f"HDF5 object is not a dataset: {normalized}")
            attributes = {str(name): _decode_attribute(value) for name, value in item.attrs.items()}
            return DatasetMetadata(
                path=normalized,
                shape=tuple(item.shape),
                dtype=str(item.dtype),
                attributes=attributes,
            )

    @property
    def periodic(self) -> PeriodicResult:
        """Return a validated typed view of a forced-periodicity result."""

        return PeriodicResult(self)

    def _open(self) -> h5py.File:
        try:
            return h5py.File(self.path, "r")
        except (OSError, ValueError) as error:
            raise ResultError(f"could not read KIM result {self.path}: {error}") from error


class PeriodicResult:
    """Validated fields and campaign-compatible operations for periodic KIM output."""

    def __init__(self, result: Result) -> None:
        self.result = result
        self._radius = self._required("backs/e/r")
        self._potential = self._required("fields/Phi")
        self._current = self._required("fields/jpar")
        width = self._required("setup/periodic_scale/dx_asis")
        self._validate_arrays(width)
        self._as_is_half_width = float(np.asarray(width))
        spacing = float(self._radius[1] - self._radius[0])
        self._resonance_radius = float(self._radius[0] + 0.5 * self._radius.size * spacing)

    @property
    def radius(self) -> NDArray[np.float64]:
        return self._radius

    @property
    def potential(self) -> NDArray[np.complex128]:
        return self._potential

    @property
    def parallel_current_density(self) -> NDArray[np.complex128]:
        return self._current

    @property
    def as_is_half_width(self) -> float:
        return self._as_is_half_width

    @property
    def resonance_radius(self) -> float:
        return self._resonance_radius

    @property
    def available_parallel_current_densities(self) -> tuple[str, ...]:
        prefix = "fields/"
        return tuple(
            path.removeprefix(prefix)
            for path in self.result.list_datasets()
            if path.startswith("fields/jpar") and "/" not in path.removeprefix(prefix)
        )

    def read_parallel_current_density(self, name: str) -> NDArray[np.complex128]:
        """Read one discovered total or species-resolved current field."""

        if name not in self.available_parallel_current_densities:
            raise ResultError(f"parallel-current dataset not found: {name}")
        current = np.asarray(self.result.read_dataset(f"fields/{name}"))
        if current.shape != self._radius.shape or not np.issubdtype(
            current.dtype, np.complexfloating
        ):
            raise ResultError(
                f"fields/{name} must be a complex array with shape {self._radius.shape}"
            )
        return current

    def integrated_parallel_current(self, region: IntegrationRegion = "as_is") -> complex:
        """Return `2*pi*integral(r*jpar dr)` over a supported periodic region."""

        if region == "as_is":
            selected = np.abs(self._radius - self._resonance_radius) <= self._as_is_half_width
            radius = self._radius[selected]
            current = self._current[selected]
        elif region == "full_window":
            radius = self._radius
            current = self._current
        else:
            raise ResultError(f"unsupported periodic integration region: {region}")
        if radius.size < 2:
            raise ResultError(
                f"periodic integration region {region!r} contains fewer than two grid points"
            )
        integrand = radius * current
        return complex(
            2.0 * np.pi * np.sum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(radius))
        )

    def _required(self, path: str) -> Any:
        try:
            return self.result.read_dataset(path)
        except ResultError as error:
            raise ResultError(
                f"required periodic dataset {path} is unavailable: {error}"
            ) from error

    def _validate_arrays(self, width: Any) -> None:
        arrays = (self._radius, self._potential, self._current)
        if not all(isinstance(value, np.ndarray) and value.ndim == 1 for value in arrays):
            raise ResultError("periodic fields must have compatible one-dimensional shapes")
        if self._radius.size < 2 or any(value.shape != self._radius.shape for value in arrays[1:]):
            raise ResultError("periodic fields must have compatible one-dimensional shapes")
        if not np.issubdtype(self._radius.dtype, np.floating) or not np.all(
            np.isfinite(self._radius)
        ):
            raise ResultError("periodic radius must contain finite floating-point values")
        if not all(
            np.issubdtype(value.dtype, np.complexfloating)
            and np.all(np.isfinite(value.real))
            and np.all(np.isfinite(value.imag))
            for value in arrays[1:]
        ):
            raise ResultError("periodic potential and current must contain finite complex values")
        spacing = np.diff(self._radius)
        if not np.all(spacing > 0.0):
            raise ResultError("periodic radius must be strictly increasing")
        if not np.allclose(spacing, spacing[0], rtol=1.0e-10, atol=1.0e-12):
            raise ResultError("periodic radius must be an equidistant periodic grid")
        width_array = np.asarray(width)
        if width_array.shape != () or not np.issubdtype(width_array.dtype, np.number):
            raise ResultError("periodic dx_asis must be a scalar")
        if not np.isfinite(width_array) or float(width_array) <= 0.0:
            raise ResultError("periodic dx_asis must be finite and positive")


def _dataset_path(path: str) -> str:
    if not isinstance(path, str) or not path:
        raise ResultError("dataset path must be a non-empty relative HDF5 path")
    candidate = PurePosixPath(path)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        raise ResultError(f"dataset path must be bounded and relative: {path!r}")
    normalized = candidate.as_posix()
    if normalized != path:
        raise ResultError(f"dataset path must be canonical: {path!r}")
    return normalized


def _decode_value(value: Any) -> Any:
    array = np.array(value, copy=True)
    fields = array.dtype.fields
    if fields is not None and set(fields) == {"real", "imag"}:
        decoded = np.asarray(array["real"]) + 1j * np.asarray(array["imag"])
        if decoded.shape == ():
            return np.complex128(decoded[()])
        return np.array(decoded, dtype=np.complex128, copy=True)
    if array.shape == ():
        scalar = array[()]
        return scalar.decode("utf-8") if isinstance(scalar, bytes) else scalar
    return array


def _decode_attribute(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "S":
            return [item.decode("utf-8") for item in value.tolist()]
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value
