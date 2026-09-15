"""Explicit source conventions and bounded conversions for KIM inputs.

The KIM executable consumes CGS profile data.  This module is the small,
deliberately strict boundary for callers that have data in one of the
approved source unit systems.  It does not infer units from values or file
names and it does not apply a Fourier conjugation or sign convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Literal, Mapping

import numpy as np
from kim.errors import KimError
from numpy.typing import ArrayLike, NDArray
from pydantic import BaseModel, ConfigDict, Field, field_validator


class ConventionError(KimError):
    """Raised when source conventions are missing, unsupported, or incompatible."""


class ConventionModel(BaseModel):
    """Strict immutable models used in conversion provenance."""

    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
        validate_default=True,
    )


_DENSITY_UNITS = {"m^-3", "1/m^3", "cm^-3", "1/cm^3"}
_RADIUS_UNITS = {"m", "cm"}
_MAGNETIC_UNITS = {"T", "G"}
_ELECTRIC_UNITS = {"V/m", "statV/cm"}
_VELOCITY_UNITS = {"m/s", "cm/s"}
_TEMPERATURE_UNITS = {"eV"}
_FREQUENCY_UNITS = {"rad/s"}


class SourceMetadata(ConventionModel):
    """Required declaration of the source representation.

    ``r_eff`` is the only radial coordinate accepted by this bounded
    conversion boundary.  The remaining fields state the source unit for a
    quantity that may be present in the converted payload.  Keeping all unit
    declarations together makes a report complete even when a particular
    optional profile is absent.
    """

    schema_version: Literal[1] = 1
    source: str = Field(min_length=1)
    coordinate: Literal["r_eff"]
    radius_unit: str
    density_unit: str
    temperature_unit: str
    magnetic_field_unit: str
    magnetic_field_sign: Literal["preserve", "signed"]
    electric_field_unit: str
    velocity_unit: str
    frequency_unit: str
    frequency_convention: Literal["signed_omega"]
    time_phase: Literal["exp(-i omega t)"]
    perturbation_phase: Literal["preserve"]
    equilibrium_provenance: str = Field(min_length=1)
    q_unit: Literal["1", "dimensionless"]
    mode_convention: Literal["signed"]
    fourier_phase: Literal["exp(+i k r)"]
    complex_phase: Literal["preserve"]
    resonance_convention: Literal["q=-m/n"]

    @field_validator("radius_unit")
    @classmethod
    def validate_radius_unit(cls, value: str) -> str:
        return _unit(value, _RADIUS_UNITS, "radius")

    @field_validator("density_unit")
    @classmethod
    def validate_density_unit(cls, value: str) -> str:
        return _unit(value, _DENSITY_UNITS, "density")

    @field_validator("temperature_unit")
    @classmethod
    def validate_temperature_unit(cls, value: str) -> str:
        return _unit(value, _TEMPERATURE_UNITS, "temperature")

    @field_validator("magnetic_field_unit")
    @classmethod
    def validate_magnetic_field_unit(cls, value: str) -> str:
        return _unit(value, _MAGNETIC_UNITS, "magnetic field")

    @field_validator("electric_field_unit")
    @classmethod
    def validate_electric_field_unit(cls, value: str) -> str:
        return _unit(value, _ELECTRIC_UNITS, "electric field")

    @field_validator("velocity_unit")
    @classmethod
    def validate_velocity_unit(cls, value: str) -> str:
        return _unit(value, _VELOCITY_UNITS, "toroidal velocity")

    @field_validator("frequency_unit")
    @classmethod
    def validate_frequency_unit(cls, value: str) -> str:
        return _unit(value, _FREQUENCY_UNITS, "frequency")


class ConversionOperation(ConventionModel):
    """One explicit operation recorded in a :class:`ConversionReport`."""

    quantity: str = Field(min_length=1)
    source_unit: str = Field(min_length=1)
    target_unit: str = Field(min_length=1)
    factor: float
    operation: str = Field(min_length=1)


class ConversionReport(ConventionModel):
    """Serializable provenance for one conversion request."""

    schema_version: Literal[1] = 1
    created_at: datetime
    source: SourceMetadata
    target: Literal["KIM-CGS"] = "KIM-CGS"
    source_hashes: dict[str, str] = Field(default_factory=dict)
    operations: tuple[ConversionOperation, ...] = ()
    rejected: tuple[str, ...] = ()

    @field_validator("source_hashes")
    @classmethod
    def source_hashes_are_sha256(cls, values: dict[str, str]) -> dict[str, str]:
        for name, digest in values.items():
            if (
                not name
                or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                raise ValueError(f"source hash for {name!r} must be a lowercase SHA-256 digest")
        return values

    @field_validator("created_at")
    @classmethod
    def report_timestamp_has_timezone(cls, value: datetime) -> datetime:
        if value.utcoffset() is None:
            raise ValueError("conversion report timestamp must include a timezone")
        return value


@dataclass(frozen=True)
class ConvertedInputs:
    """Converted arrays and their inspectable provenance."""

    values: Mapping[str, NDArray[Any]]
    report: ConversionReport


_Quantity = Literal[
    "radius",
    "density",
    "temperature",
    "magnetic_field",
    "electric_field",
    "toroidal_velocity",
    "frequency",
]


_QUANTITY_UNITS: dict[str, tuple[str, str, float]] = {
    "radius": ("radius_unit", "cm", 100.0),
    "density": ("density_unit", "1/cm^3", 1.0e-6),
    "temperature": ("temperature_unit", "eV", 1.0),
    "magnetic_field": ("magnetic_field_unit", "G", 1.0e4),
    # 1 statV/cm = 29979.2458 V/m.
    "electric_field": ("electric_field_unit", "statV/cm", 1.0 / 29979.2458),
    "toroidal_velocity": ("velocity_unit", "cm/s", 100.0),
    "frequency": ("frequency_unit", "rad/s", 1.0),
}

_KEY_QUANTITIES: dict[str, str] = {
    "r": "radius",
    "radius": "radius",
    "r_eff": "radius",
    "n_e": "density",
    "electron_density": "density",
    "density": "density",
    "Te": "temperature",
    "Ti": "temperature",
    "temperature": "temperature",
    "B": "magnetic_field",
    "btor": "magnetic_field",
    "br": "magnetic_field",
    "Br": "magnetic_field",
    "Er": "electric_field",
    "electric_field": "electric_field",
    "Vz": "toroidal_velocity",
    "toroidal_velocity": "toroidal_velocity",
    "omega": "frequency",
    "frequency": "frequency",
}


def convert_quantity(
    values: ArrayLike,
    quantity: _Quantity | str,
    metadata: SourceMetadata,
) -> NDArray[Any]:
    """Convert one finite source array to KIM's CGS representation.

    The returned array is a detached NumPy array.  Complex values retain
    their real and imaginary parts exactly apart from the approved scalar
    unit factor.
    """

    _source_field, _target_unit, factor = _quantity_spec(quantity, metadata)
    array = np.asarray(values)
    _require_finite(array, str(quantity))
    converted = array * factor
    _require_finite(converted, str(quantity))
    return np.array(converted, copy=True)


def convert_from_kim(
    values: ArrayLike,
    quantity: _Quantity | str,
    metadata: SourceMetadata,
) -> NDArray[Any]:
    """Convert KIM CGS values back to the declared source units.

    This inverse is provided for inspection and round-trip verification.  It
    never broadens the set of accepted source units or changes conventions.
    """

    source_field, target_unit, factor = _quantity_spec(quantity, metadata)
    source_unit = getattr(metadata, source_field)
    array = np.asarray(values)
    _require_finite(array, str(quantity))
    inverse = 1.0 / factor
    converted = array * inverse
    _require_finite(converted, str(quantity))
    return np.array(converted, copy=True)


def convert_inputs(
    values: Mapping[str, ArrayLike],
    metadata: SourceMetadata,
    *,
    created_at: datetime | None = None,
    source_hashes: Mapping[str, str] | None = None,
) -> ConvertedInputs:
    """Convert named source values and return an ordered conversion report.

    Names are intentionally limited to the profile/setup quantities used by
    KIM.  ``m`` and ``n`` mode numbers are copied as signed integers and do
    not pass through a unit conversion.
    """

    converted: dict[str, NDArray[Any]] = {}
    operations: list[ConversionOperation] = []
    for name, raw in values.items():
        if name in {"m", "n", "m_mode", "n_mode"}:
            array = np.asarray(raw)
            _require_mode_values(array, name)
            converted[name] = np.array(array, copy=True)
            continue
        quantity = _KEY_QUANTITIES.get(name)
        if quantity is None:
            raise ConventionError(
                f"unsupported input quantity {name!r}; declare one of {sorted(_KEY_QUANTITIES)}"
            )
        source_field, target_unit, factor = _quantity_spec(quantity, metadata)
        source_unit = getattr(metadata, source_field)
        converted[name] = convert_quantity(raw, quantity, metadata)
        operations.append(
            ConversionOperation(
                quantity=name,
                source_unit=source_unit,
                target_unit=target_unit,
                factor=factor,
                operation=("preserve" if factor == 1.0 else f"multiply by {factor:.16g}"),
            )
        )

    timestamp = created_at or datetime.now(timezone.utc)
    report = ConversionReport(
        created_at=timestamp,
        source=metadata,
        source_hashes=dict(source_hashes or {}),
        operations=tuple(operations),
    )
    return ConvertedInputs(values=converted, report=report)


def resonance_target(m_mode: int, n_mode: int, metadata: SourceMetadata) -> float:
    """Return KIM's signed resonance target ``q=-m/n``."""

    _require_conventions(metadata)
    if not isinstance(m_mode, (int, np.integer)) or not isinstance(n_mode, (int, np.integer)):
        raise ConventionError("m_mode and n_mode must be signed integers")
    if m_mode == 0 or n_mode == 0:
        raise ConventionError("m_mode and n_mode must be nonzero for q=-m/n")
    return -int(m_mode) / int(n_mode)


def reconstruct_fourier(
    radius: ArrayLike,
    modes: ArrayLike,
    coefficients: ArrayLike,
    metadata: SourceMetadata,
    *,
    period: float | None = None,
) -> NDArray[np.complex128]:
    """Reconstruct ``sum Phi_l exp(+i*k_l*r)`` without changing phase.

    ``modes`` are the signed integer Fourier indices and ``period`` is one
    full period in the source radius unit.  The physical wavenumber is
    ``k_l = 2*pi*l/period`` after both radius and period are converted to cm.
    """

    _require_conventions(metadata)
    if period is None:
        raise ConventionError("period is required for Fourier reconstruction")
    period_array = convert_quantity(period, "radius", metadata)
    period_cm = float(period_array)
    if not np.isfinite(period_cm) or period_cm <= 0.0:
        raise ConventionError("period must be a finite positive radius")
    radial = convert_quantity(radius, "radius", metadata).astype(np.float64, copy=False)
    indices = np.asarray(modes)
    amplitudes = np.asarray(coefficients, dtype=np.complex128)
    if indices.ndim != 1 or amplitudes.ndim != 1 or indices.size != amplitudes.size:
        raise ConventionError(
            "modes and coefficients must be one-dimensional arrays of equal length"
        )
    _require_mode_values(indices, "modes", allow_zero=True)
    _require_finite(amplitudes, "coefficients")
    wavenumbers = 2.0 * np.pi * indices / period_cm
    return np.exp(1j * np.outer(radial, wavenumbers)) @ amplitudes


def reconstruct_fourier_derivative(
    radius: ArrayLike,
    modes: ArrayLike,
    coefficients: ArrayLike,
    metadata: SourceMetadata,
    *,
    period: float | None = None,
) -> NDArray[np.complex128]:
    """Evaluate the radial derivative of :func:`reconstruct_fourier`."""

    _require_conventions(metadata)
    if period is None:
        raise ConventionError("period is required for Fourier reconstruction")
    period_cm = float(convert_quantity(period, "radius", metadata))
    if not np.isfinite(period_cm) or period_cm <= 0.0:
        raise ConventionError("period must be a finite positive radius")
    radial = convert_quantity(radius, "radius", metadata).astype(np.float64, copy=False)
    indices = np.asarray(modes)
    amplitudes = np.asarray(coefficients, dtype=np.complex128)
    if indices.ndim != 1 or amplitudes.ndim != 1 or indices.size != amplitudes.size:
        raise ConventionError(
            "modes and coefficients must be one-dimensional arrays of equal length"
        )
    _require_mode_values(indices, "modes", allow_zero=True)
    _require_finite(amplitudes, "coefficients")
    wavenumbers = 2.0 * np.pi * indices / period_cm
    return np.exp(1j * np.outer(radial, wavenumbers)) @ (1j * wavenumbers * amplitudes)


def _quantity_spec(quantity: str, metadata: SourceMetadata) -> tuple[str, str, float]:
    try:
        source_field, target_unit, factor = _QUANTITY_UNITS[quantity]
    except KeyError as error:
        raise ConventionError(f"unsupported quantity {quantity!r}") from error
    _require_conventions(metadata)
    source_unit = getattr(metadata, source_field)
    expected = _expected_units(quantity)
    if source_unit not in expected:
        raise ConventionError(
            f"dimension mismatch for {quantity}: source unit {source_unit!r} is not supported; "
            f"expected one of {sorted(expected)}"
        )
    if quantity == "frequency" and source_unit != "rad/s":
        raise ConventionError("frequency must be declared in rad/s; Hz is not supported")
    if quantity == "temperature" and source_unit != "eV":
        raise ConventionError("temperature must be declared in eV; Kelvin is not supported")
    # ``cm^-3`` and ``1/cm^3`` are equivalent declarations.  Keep the
    # original spelling in provenance while applying the dimensional factor
    # only when the source is genuinely in SI units.
    if quantity == "density" and source_unit in {"cm^-3", "1/cm^3"}:
        factor = 1.0
    elif source_unit == target_unit:
        factor = 1.0
    return source_field, target_unit, factor


def _expected_units(quantity: str) -> set[str]:
    if quantity == "radius":
        return _RADIUS_UNITS
    if quantity == "density":
        return _DENSITY_UNITS
    if quantity == "temperature":
        return _TEMPERATURE_UNITS
    if quantity == "magnetic_field":
        return _MAGNETIC_UNITS
    if quantity == "electric_field":
        return _ELECTRIC_UNITS
    if quantity == "toroidal_velocity":
        return _VELOCITY_UNITS
    if quantity == "frequency":
        return _FREQUENCY_UNITS
    raise ConventionError(f"unsupported quantity {quantity!r}")


def _require_conventions(metadata: SourceMetadata) -> None:
    if metadata.coordinate != "r_eff":
        raise ConventionError("only the r_eff radial coordinate is supported")
    if metadata.mode_convention != "signed":
        raise ConventionError("mode numbers must use the signed convention")
    if metadata.fourier_phase != "exp(+i k r)":
        raise ConventionError("Fourier phase must be declared as exp(+i k r)")
    if metadata.complex_phase != "preserve":
        raise ConventionError("complex phase must be preserved")
    if metadata.resonance_convention != "q=-m/n":
        raise ConventionError("resonance must use signed q=-m/n")
    if metadata.magnetic_field_sign not in {"preserve", "signed"}:
        raise ConventionError("magnetic field sign must be preserved")
    if metadata.frequency_convention != "signed_omega":
        raise ConventionError("frequency must use signed omega")
    if metadata.time_phase != "exp(-i omega t)":
        raise ConventionError("time phase must be exp(-i omega t)")
    if metadata.perturbation_phase != "preserve":
        raise ConventionError("perturbation phase must be preserved")
    if metadata.q_unit not in {"1", "dimensionless"}:
        raise ConventionError("q/profile unit must be dimensionless")


def _unit(value: str, allowed: set[str], quantity: str) -> str:
    if value not in allowed:
        raise ValueError(
            f"unsupported {quantity} unit {value!r}; expected one of {sorted(allowed)}"
        )
    return value


def _require_finite(values: NDArray[Any], quantity: str) -> None:
    try:
        finite = np.isfinite(values)
    except TypeError as error:
        raise ConventionError(f"{quantity} values must be numeric and finite") from error
    if not np.all(finite):
        raise ConventionError(f"{quantity} values must be finite")


def _require_mode_values(values: NDArray[Any], quantity: str, *, allow_zero: bool = False) -> None:
    _require_finite(values, quantity)
    if not np.issubdtype(values.dtype, np.signedinteger):
        raise ConventionError(f"{quantity} must contain signed integer mode numbers")
    if not allow_zero and np.any(values == 0):
        raise ConventionError(f"{quantity} must be nonzero")


__all__ = [
    "ConversionOperation",
    "ConversionReport",
    "ConvertedInputs",
    "ConventionError",
    "SourceMetadata",
    "convert_inputs",
    "convert_from_kim",
    "convert_quantity",
    "reconstruct_fourier",
    "reconstruct_fourier_derivative",
    "resonance_target",
]
