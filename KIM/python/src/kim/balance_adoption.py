"""Pure, explicit transformations for BALANCE profile adoption.

This module is intentionally narrower than :mod:`kim.conventions`.  It accepts
only the conventions used by the approved BALANCE source files and returns
detached arrays together with the operation needed to explain each conversion.
It does not read or write files, map radial coordinates, or infer scientific
metadata from values.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from kim.conventions import ConventionError, ConversionOperation
from numpy.typing import ArrayLike, NDArray

_BALANCE_COORDINATES = {"rho_pol", "sqrt_psiN"}
_BALANCE_COORDINATE_UNIT = "1"


@dataclass(frozen=True)
class BalanceTransformation:
    """Detached BALANCE arrays and their explicit conversion operation."""

    coordinate: NDArray[np.float64]
    values: NDArray[np.float64]
    operation: ConversionOperation


def convert_balance_density(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
) -> BalanceTransformation:
    """Convert an explicitly declared BALANCE density to ``1/cm^3``."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if source_unit in {"1/m^3", "m^-3"}:
        factor = 1.0e-6
        operation = "density_si_to_cgs"
    elif source_unit in {"1/cm^3", "cm^-3"}:
        factor = 1.0
        operation = "preserve"
    else:
        raise ConventionError(
            f"density unit {source_unit!r} is unsupported; declare SI or CGS density explicitly"
        )

    converted = _scaled_values(value_array, factor, "density")
    return BalanceTransformation(
        coordinate=coordinate_array,
        values=converted,
        operation=ConversionOperation(
            quantity="density",
            source_unit=source_unit,
            target_unit="1/cm^3",
            factor=factor,
            operation=operation,
        ),
    )


def convert_balance_temperature(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    species: Literal["electron", "ion"] | str,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
) -> BalanceTransformation:
    """Validate an explicitly declared BALANCE temperature in eV."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if species not in {"electron", "ion"}:
        raise ConventionError("temperature species must be explicitly declared as electron or ion")
    if source_unit != "eV":
        raise ConventionError("temperature unit must be declared in eV; conversion is not inferred")

    return BalanceTransformation(
        coordinate=coordinate_array,
        values=_readonly_copy(value_array),
        operation=ConversionOperation(
            quantity=f"{species}_temperature",
            source_unit="eV",
            target_unit="eV",
            factor=1.0,
            operation="preserve",
        ),
    )


def convert_balance_rotation(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
    r0_cm: float | None,
) -> BalanceTransformation:
    """Convert signed angular rotation to signed toroidal velocity in cm/s."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if source_unit != "rad/s":
        raise ConventionError(
            "rotation unit must be declared in rad/s; velocity conversion is not inferred"
        )
    try:
        radius = float(r0_cm)
    except (OverflowError, TypeError, ValueError) as error:
        raise ConventionError("r0_cm must be a finite positive radius") from error
    if not np.isfinite(radius) or radius <= 0.0:
        raise ConventionError("r0_cm must be a finite positive radius")

    converted = _scaled_values(value_array, radius, "rotation")
    return BalanceTransformation(
        coordinate=coordinate_array,
        values=converted,
        operation=ConversionOperation(
            quantity="toroidal_velocity",
            source_unit="rad/s",
            target_unit="cm/s",
            factor=radius,
            operation="omega_to_v_phi",
            parameters={"major_radius_cm": radius},
        ),
    )


def convert_balance_q(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
    operation: Literal["preserve", "negate"] | str,
) -> BalanceTransformation:
    """Apply an explicitly selected signed-q operation."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if source_unit != "1":
        raise ConventionError("q unit must be declared as dimensionless unit 1")
    if operation not in {"preserve", "negate"}:
        raise ConventionError("q operation must be explicitly selected as preserve or negate")

    factor = 1.0 if operation == "preserve" else -1.0
    converted = _scaled_values(value_array, factor, "q")
    return BalanceTransformation(
        coordinate=coordinate_array,
        values=converted,
        operation=ConversionOperation(
            quantity="q",
            source_unit="1",
            target_unit="1",
            factor=factor,
            operation=operation,
        ),
    )


def _validate_profile(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if coordinate not in _BALANCE_COORDINATES:
        raise ConventionError("coordinate must be explicitly declared as rho_pol or sqrt_psiN")
    if coordinate_unit != _BALANCE_COORDINATE_UNIT:
        raise ConventionError("coordinate_unit must be 1 for rho_pol or sqrt_psiN")

    # Strict BALANCE readers enforce strictly increasing coordinates upstream.
    coordinate_array = _as_finite_vector(coordinate_values, "coordinate")
    value_array = _as_finite_vector(values, "profile")
    if coordinate_array.size != value_array.size:
        raise ConventionError("coordinate and values must have equal lengths")
    return coordinate_array, value_array


def _as_finite_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as error:
        raise ConventionError(f"{name} values must be numeric and finite") from error
    if np.iscomplexobj(raw) or (
        raw.dtype.kind == "O" and any(np.iscomplexobj(item) for item in raw.flat)
    ):
        raise ConventionError(f"{name} values must be real; complex values are unsupported")
    try:
        array = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ConventionError(f"{name} values must be numeric and finite") from error
    if array.ndim != 1:
        raise ConventionError(f"{name} values must be one-dimensional")
    if not np.all(np.isfinite(array)):
        raise ConventionError(f"{name} values must be finite")
    result = np.array(array, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _scaled_values(values: NDArray[np.float64], factor: float, name: str) -> NDArray[np.float64]:
    converted = values * factor
    if not np.all(np.isfinite(converted)):
        raise ConventionError(f"converted {name} values must be finite")
    result = np.array(converted, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _readonly_copy(values: NDArray[np.float64]) -> NDArray[np.float64]:
    result = np.array(values, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


__all__ = [
    "BalanceTransformation",
    "convert_balance_density",
    "convert_balance_q",
    "convert_balance_rotation",
    "convert_balance_temperature",
]
