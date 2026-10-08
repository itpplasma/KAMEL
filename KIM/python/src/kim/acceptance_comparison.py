"""Pure, in-memory comparison of an oracle and a prepared radial profile.

The oracle is always the reference profile.  This module intentionally has no
file, HDF5, plotting, or solver dependencies: it reports measurements and
coverage facts for callers that already have validated arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np
from kim.errors import KimError
from numpy.typing import ArrayLike, NDArray


class ComparisonError(KimError):
    """Raised when a profile comparison request is invalid."""


@dataclass(frozen=True)
class ComparisonMeasurements:
    """Absolute and oracle-relative error measurements on the comparison grid."""

    absolute_rms: float
    absolute_max: float
    relative_rms: float | None
    relative_max: float | None


@dataclass(frozen=True)
class ExclusionInterval:
    """One typed interval excluded from a profile comparison.

    Endpoint flags are explicit because domain and shared-overlap exclusions
    meet at boundaries without sharing the boundary point.
    """

    lower_cm: float
    upper_cm: float
    lower_inclusive: bool
    upper_inclusive: bool


@dataclass(frozen=True)
class ProfileExclusions:
    """Points and continuous intervals excluded for one source profile."""

    outside_domain_points_cm: tuple[float, ...]
    outside_domain_intervals_cm: tuple[ExclusionInterval, ...]
    outside_shared_overlap_points_cm: tuple[float, ...]
    outside_shared_overlap_intervals_cm: tuple[ExclusionInterval, ...]


@dataclass(frozen=True)
class ComparisonExclusions:
    """Exclusions separated by source profile."""

    oracle: ProfileExclusions
    prepared: ProfileExclusions


@dataclass(frozen=True)
class ResonanceComparison:
    """Independent piecewise-linear resonance information for both profiles.

    ``*_crossing_covered`` means at least one crossing lies in the requested
    domain. It does not assert shared profile support at that crossing;
    consult the exclusions and evaluated radii for comparison coverage.
    """

    target_q: float
    reference_crossing_radii_cm: tuple[float, ...]
    candidate_crossing_radii_cm: tuple[float, ...]
    reference_crossing_covered: bool
    candidate_crossing_covered: bool


@dataclass(frozen=True)
class ComparisonResult:
    """Complete comparison report, with no implicit acceptance policy."""

    comparison_radius_cm: NDArray[np.float64]
    measurements: ComparisonMeasurements
    exclusions: ComparisonExclusions
    warnings: tuple[str, ...]
    threshold_decisions: Mapping[str, bool] | None
    overall_pass: bool | None
    resonance: ResonanceComparison | None = None


_METRIC_NAMES = ("absolute_rms", "absolute_max", "relative_rms", "relative_max")
_DIRECTIONS = ("prepared_to_oracle", "oracle_to_prepared")


def compare_profiles(
    *,
    oracle_radius_cm: ArrayLike,
    oracle_values: ArrayLike,
    prepared_radius_cm: ArrayLike,
    prepared_values: ArrayLike,
    domain_cm: tuple[float, float],
    interpolation_direction: str,
    method: str,
    relative_floor: float | None,
    tolerances: Mapping[str, float] | None = None,
    resonance: tuple[int, int] | float | None = None,
) -> ComparisonResult:
    """Compare two finite, strictly increasing piecewise-linear profiles.

    ``interpolation_direction`` names the source-to-target interpolation
    direction (``prepared_to_oracle`` evaluates on oracle nodes and
    ``oracle_to_prepared`` evaluates on prepared nodes).  Only existing target
    nodes within the requested domain and mutual source support are evaluated;
    endpoints are never synthesized and no extrapolation is performed.
    """

    oracle_radius, oracle_profile = _validate_profile(oracle_radius_cm, oracle_values, "oracle")
    prepared_radius, prepared_profile = _validate_profile(
        prepared_radius_cm, prepared_values, "prepared"
    )
    domain_lower, domain_upper = _validate_domain(domain_cm)
    direction = _validate_direction(interpolation_direction)
    _validate_method(method)
    floor = None if relative_floor is None else _validate_floor(relative_floor)
    validated_tolerances = _validate_tolerances(tolerances)
    if (
        floor is None
        and validated_tolerances is not None
        and any(metric.startswith("relative_") for metric in validated_tolerances)
    ):
        raise ComparisonError("relative tolerances require a relative floor")
    target_q = _validate_resonance(resonance)

    shared_lower = max(float(oracle_radius[0]), float(prepared_radius[0]))
    shared_upper = min(float(oracle_radius[-1]), float(prepared_radius[-1]))
    if shared_lower > shared_upper:
        raise ComparisonError("requested domain has no continuous profile overlap")
    if domain_upper < shared_lower or domain_lower > shared_upper:
        raise ComparisonError("requested domain has no continuous profile overlap")

    if direction == "prepared_to_oracle":
        target_radius = oracle_radius
        oracle_on_grid = oracle_profile
        prepared_on_grid = None
        reference_radius = prepared_radius
        reference_values = prepared_profile
    else:
        target_radius = prepared_radius
        prepared_on_grid = prepared_profile
        oracle_on_grid = None
        reference_radius = oracle_radius
        reference_values = oracle_profile

    target_mask = (
        (target_radius >= domain_lower)
        & (target_radius <= domain_upper)
        & (target_radius >= shared_lower)
        & (target_radius <= shared_upper)
    )
    comparison_radius = target_radius[target_mask]
    if comparison_radius.size == 0:
        raise ComparisonError("continuous overlap contains no comparison point")

    # Source support has already been applied to the target mask.  The bounded
    # interpolator therefore evaluates only inside its source grid.
    interpolated_values = _linear_interpolate_bounded(
        comparison_radius, reference_radius, reference_values
    )
    if direction == "prepared_to_oracle":
        assert oracle_on_grid is not None
        oracle_values_on_grid = oracle_on_grid[target_mask]
        prepared_values_on_grid = interpolated_values
    else:
        assert prepared_on_grid is not None
        oracle_values_on_grid = interpolated_values
        prepared_values_on_grid = prepared_on_grid[target_mask]
    # Keep the arithmetic inside the floating-point domain long enough to
    # report a useful domain error rather than leaking inf/nan metrics.  The
    # scaled RMS below avoids squaring the original values.
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        error = prepared_values_on_grid - oracle_values_on_grid
    _validate_computed_array(error, "comparison error")
    relative_rms = None
    relative_max = None
    if floor is not None:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
            relative_error = error / np.maximum(np.abs(oracle_values_on_grid), floor)
        _validate_computed_array(relative_error, "relative comparison error")
        relative_rms = _stable_rms(relative_error, "relative RMS")
        relative_max = float(np.max(np.abs(relative_error)))
    measurements = ComparisonMeasurements(
        absolute_rms=_stable_rms(error, "absolute RMS"),
        absolute_max=float(np.max(np.abs(error))),
        relative_rms=relative_rms,
        relative_max=relative_max,
    )
    _validate_measurements(measurements)

    oracle_exclusions = _profile_exclusions(
        oracle_radius,
        domain_lower,
        domain_upper,
        shared_lower,
        shared_upper,
    )
    prepared_exclusions = _profile_exclusions(
        prepared_radius,
        domain_lower,
        domain_upper,
        shared_lower,
        shared_upper,
    )
    exclusions = ComparisonExclusions(oracle_exclusions, prepared_exclusions)

    resonance_report = None
    if target_q is not None:
        reference_crossings = _piecewise_linear_roots(oracle_radius, oracle_profile, target_q)
        candidate_crossings = _piecewise_linear_roots(prepared_radius, prepared_profile, target_q)
        resonance_report = ResonanceComparison(
            target_q=target_q,
            reference_crossing_radii_cm=reference_crossings,
            candidate_crossing_radii_cm=candidate_crossings,
            reference_crossing_covered=_covered(reference_crossings, domain_lower, domain_upper),
            candidate_crossing_covered=_covered(candidate_crossings, domain_lower, domain_upper),
        )

    warnings = _warnings(exclusions, resonance_report)
    if floor is None:
        warnings = (*warnings, "relative metrics unavailable: no relative floor was supplied")
    decisions = None
    overall_pass = None
    if validated_tolerances is not None:
        decisions_dict = {
            metric: getattr(measurements, metric) <= limit
            for metric, limit in validated_tolerances.items()
        }
        decisions = MappingProxyType(decisions_dict)
        overall_pass = all(decisions_dict.values())

    return ComparisonResult(
        comparison_radius_cm=_immutable_array(comparison_radius),
        measurements=measurements,
        exclusions=exclusions,
        warnings=warnings,
        threshold_decisions=decisions,
        overall_pass=overall_pass,
        resonance=resonance_report,
    )


def _validate_profile(
    radius: ArrayLike, values: ArrayLike, name: str
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    radius_array = _as_float_array(radius, f"{name} radius")
    values_array = _as_float_array(values, f"{name} values")
    if radius_array.ndim != 1 or values_array.ndim != 1:
        raise ComparisonError(f"{name} radius and values must be one-dimensional")
    if radius_array.size < 2 or values_array.size < 2:
        raise ComparisonError(f"{name} profile must contain at least two samples")
    if radius_array.size != values_array.size:
        raise ComparisonError(f"{name} radius and values must have the same length")
    if not np.all(np.isfinite(radius_array)) or not np.all(np.isfinite(values_array)):
        raise ComparisonError(f"{name} radius and values must contain only finite values")
    if not np.all(radius_array[1:] > radius_array[:-1]):
        raise ComparisonError(f"{name} radius must be strictly increasing without duplicates")
    return radius_array, values_array


def _as_float_array(value: ArrayLike, name: str) -> NDArray[np.float64]:
    try:
        candidate = np.asarray(value)
    except Exception as error:
        raise ComparisonError(f"{name} must be a numeric array") from error
    if candidate.ndim != 1:
        raise ComparisonError(f"{name} must be one-dimensional")
    if not np.issubdtype(candidate.dtype, np.number) or np.issubdtype(
        candidate.dtype, np.complexfloating
    ):
        raise ComparisonError(f"{name} must be a real numeric array")
    try:
        return np.asarray(candidate, dtype=np.float64).copy()
    except (TypeError, ValueError, OverflowError) as error:
        raise ComparisonError(f"{name} must be a real numeric array") from error


def _validate_domain(domain: tuple[float, float]) -> tuple[float, float]:
    try:
        if len(domain) != 2:
            raise ValueError
    except (TypeError, ValueError, IndexError, OverflowError) as error:
        raise ComparisonError("domain must contain two finite endpoints") from error
    lower = _real_scalar(domain[0], "domain endpoint")
    upper = _real_scalar(domain[1], "domain endpoint")
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise ComparisonError("domain must be finite and non-empty with lower < upper")
    return lower, upper


def _validate_direction(direction: str) -> str:
    if not isinstance(direction, str):
        raise ComparisonError("interpolation direction must be a scalar string")
    if direction not in _DIRECTIONS:
        raise ComparisonError(
            "interpolation direction must be explicitly 'prepared_to_oracle' or "
            "'oracle_to_prepared'"
        )
    return direction


def _validate_method(method: str) -> None:
    if not isinstance(method, str):
        raise ComparisonError("interpolation method must be a scalar string")
    if method != "linear":
        raise ComparisonError("only explicit linear interpolation method is supported")


def _validate_floor(value: float) -> float:
    try:
        floor = _real_scalar(value, "relative floor")
    except ComparisonError as error:
        raise ComparisonError("relative floor must be finite and positive") from error
    if not np.isfinite(floor) or floor <= 0.0:
        raise ComparisonError("relative floor must be finite and positive")
    return floor


def _validate_tolerances(
    tolerances: Mapping[str, float] | None,
) -> dict[str, float] | None:
    if tolerances is None:
        return None
    if not isinstance(tolerances, Mapping) or not tolerances:
        raise ComparisonError("tolerance mapping must be non-empty")
    try:
        keys = tuple(tolerances.keys())
    except Exception as error:
        raise ComparisonError("tolerance mapping keys must be strings") from error
    if any(not isinstance(key, str) for key in keys):
        raise ComparisonError("tolerance mapping keys must be strings")
    unknown = tuple(key for key in keys if key not in _METRIC_NAMES)
    if unknown:
        raise ComparisonError(f"tolerance mapping contains unknown metric: {unknown!r}")
    validated: dict[str, float] = {}
    for metric in _METRIC_NAMES:
        if metric not in tolerances:
            continue
        value = tolerances[metric]
        try:
            limit = _real_scalar(value, f"tolerance for {metric}")
        except ComparisonError as error:
            raise ComparisonError(
                f"tolerance for {metric} must be finite and non-negative"
            ) from error
        if not np.isfinite(limit) or limit < 0.0:
            raise ComparisonError(f"tolerance for {metric} must be finite and non-negative")
        validated[metric] = limit
    return validated


def _validate_resonance(resonance: tuple[int, int] | float | None) -> float | None:
    if resonance is None:
        return None
    if isinstance(resonance, tuple):
        if len(resonance) != 2:
            raise ComparisonError("resonance mode must contain m and n")
        m_mode, n_mode = resonance
        if (
            isinstance(m_mode, (bool, np.bool_))
            or isinstance(n_mode, (bool, np.bool_))
            or not isinstance(m_mode, (int, np.integer))
            or not isinstance(n_mode, (int, np.integer))
            or m_mode == 0
            or n_mode == 0
        ):
            raise ComparisonError("resonance mode numbers must be nonzero integers")
        try:
            target = -float(m_mode) / float(n_mode)
        except (OverflowError, ZeroDivisionError) as error:
            raise ComparisonError("resonance target must be a finite signed number") from error
        if not np.isfinite(target):
            raise ComparisonError("resonance target must be a finite signed number")
        return target
    try:
        target = _real_scalar(resonance, "resonance target")
    except ComparisonError as error:
        raise ComparisonError("resonance target must be a finite signed number") from error
    if not np.isfinite(target):
        raise ComparisonError("resonance target must be a finite signed number")
    return target


def _real_scalar(value: object, name: str) -> float:
    """Convert one finite-checkable real scalar without accepting text or arrays."""

    if isinstance(value, (bool, np.bool_, str, bytes)):
        raise ComparisonError(f"{name} must be a real scalar")
    try:
        candidate = np.asarray(value)
    except Exception as error:
        raise ComparisonError(f"{name} must be a real scalar") from error
    if candidate.ndim != 0 or not np.issubdtype(candidate.dtype, np.number):
        raise ComparisonError(f"{name} must be a real scalar")
    if np.issubdtype(candidate.dtype, np.complexfloating):
        raise ComparisonError(f"{name} must be a real scalar")
    try:
        return float(candidate)
    except (TypeError, ValueError, OverflowError) as error:
        raise ComparisonError(f"{name} must be a real scalar") from error


def _validate_computed_array(values: NDArray[np.float64], name: str) -> None:
    if not np.all(np.isfinite(values)):
        raise ComparisonError(f"{name} is not finite and cannot be represented")


def _stable_rms(values: NDArray[np.float64], name: str) -> float:
    """Compute RMS without squaring values at their original scale."""

    scale = float(np.max(np.abs(values)))
    if not np.isfinite(scale):
        raise ComparisonError(f"{name} is not finite and cannot be represented")
    if scale == 0.0:
        return 0.0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
        scaled = np.abs(values) / scale
        scaled_mean_square = float(np.mean(scaled * scaled))
        result = scale * float(np.sqrt(scaled_mean_square))
    if not np.isfinite(result):
        raise ComparisonError(f"{name} is not finite and cannot be represented")
    return result


def _linear_interpolate_bounded(
    target_radius: NDArray[np.float64],
    source_radius: NDArray[np.float64],
    source_values: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Linearly interpolate bounded targets without overflow-prone spans."""

    if target_radius.size == 0:
        return np.empty(0, dtype=np.float64)
    source_count = source_radius.size
    insertion_points = np.searchsorted(source_radius, target_radius, side="left")
    interpolated = np.empty(target_radius.size, dtype=np.float64)
    for position, target in enumerate(target_radius):
        insertion = int(insertion_points[position])
        if insertion < source_count and target == source_radius[insertion]:
            interpolated[position] = source_values[insertion]
            continue
        lower_index = insertion - 1
        if lower_index < 0 or lower_index >= source_count - 1:
            raise ComparisonError("interpolation target lies outside source support")
        lower_radius = float(source_radius[lower_index])
        upper_radius = float(source_radius[lower_index + 1])
        fraction = _bounded_fraction(float(target), lower_radius, upper_radius)
        interpolated[position] = _stable_convex_value(
            float(source_values[lower_index]),
            float(source_values[lower_index + 1]),
            fraction,
        )
    _validate_computed_array(interpolated, "interpolated profile")
    return interpolated


def _bounded_fraction(target: float, lower: float, upper: float) -> float:
    """Return the bounded interval fraction without extreme-coordinate overflow."""

    if lower < 0.0 < upper:
        scale = max(abs(lower), abs(upper), abs(target))
        lower_scaled = lower / scale
        upper_scaled = upper / scale
        target_scaled = target / scale
        fraction = (target_scaled - lower_scaled) / (upper_scaled - lower_scaled)
    else:
        fraction = (target - lower) / (upper - lower)
    if not np.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
        raise ComparisonError("interpolation fraction cannot be represented")
    return float(fraction)


def _stable_convex_value(lower: float, upper: float, fraction: float) -> float:
    """Blend finite values as a normalized convex combination."""

    scale = max(abs(lower), abs(upper))
    if scale == 0.0:
        return 0.0
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        lower_scaled = lower / scale
        upper_scaled = upper / scale
        value = scale * (lower_scaled * (1.0 - fraction) + upper_scaled * fraction)
    if not np.isfinite(value):
        raise ComparisonError("interpolated value cannot be represented")
    return float(value)


def _validate_measurements(measurements: ComparisonMeasurements) -> None:
    absolute_values = np.asarray(
        (measurements.absolute_rms, measurements.absolute_max), dtype=np.float64
    )
    _validate_computed_array(absolute_values, "absolute comparison measurements")
    relative_values = tuple(
        value
        for value in (measurements.relative_rms, measurements.relative_max)
        if value is not None
    )
    if relative_values:
        _validate_computed_array(np.asarray(relative_values, dtype=np.float64), "relative metrics")


def _profile_exclusions(
    radius: NDArray[np.float64],
    domain_lower: float,
    domain_upper: float,
    shared_lower: float,
    shared_upper: float,
) -> ProfileExclusions:
    outside_domain_points = tuple(
        float(point) for point in radius if point < domain_lower or point > domain_upper
    )
    outside_domain_intervals = _domain_intervals(
        float(radius[0]), float(radius[-1]), domain_lower, domain_upper
    )
    outside_shared_points = tuple(
        float(point)
        for point in radius
        if domain_lower <= point <= domain_upper and (point < shared_lower or point > shared_upper)
    )
    outside_shared_intervals = _shared_intervals(
        float(radius[0]),
        float(radius[-1]),
        domain_lower,
        domain_upper,
        shared_lower,
        shared_upper,
    )
    return ProfileExclusions(
        outside_domain_points_cm=outside_domain_points,
        outside_domain_intervals_cm=outside_domain_intervals,
        outside_shared_overlap_points_cm=outside_shared_points,
        outside_shared_overlap_intervals_cm=outside_shared_intervals,
    )


def _domain_intervals(
    radius_lower: float, radius_upper: float, domain_lower: float, domain_upper: float
) -> tuple[ExclusionInterval, ...]:
    intervals: list[ExclusionInterval] = []
    left_upper = min(radius_upper, domain_lower)
    if radius_lower < left_upper:
        intervals.append(
            ExclusionInterval(
                radius_lower,
                left_upper,
                True,
                left_upper < domain_lower,
            )
        )
    right_lower = max(radius_lower, domain_upper)
    if right_lower < radius_upper:
        intervals.append(
            ExclusionInterval(
                right_lower,
                radius_upper,
                right_lower > domain_upper,
                True,
            )
        )
    return tuple(intervals)


def _shared_intervals(
    radius_lower: float,
    radius_upper: float,
    domain_lower: float,
    domain_upper: float,
    shared_lower: float,
    shared_upper: float,
) -> tuple[ExclusionInterval, ...]:
    intervals: list[ExclusionInterval] = []
    left_lower = max(radius_lower, domain_lower)
    left_upper = min(radius_upper, domain_upper, shared_lower)
    if left_lower < left_upper:
        intervals.append(
            ExclusionInterval(
                left_lower,
                left_upper,
                True,
                left_upper < shared_lower,
            )
        )
    right_lower = max(radius_lower, domain_lower, shared_upper)
    right_upper = min(radius_upper, domain_upper)
    if right_lower < right_upper:
        intervals.append(
            ExclusionInterval(
                right_lower,
                right_upper,
                right_lower > shared_upper,
                True,
            )
        )
    return tuple(intervals)


def _piecewise_linear_roots(
    radius: NDArray[np.float64], values: NDArray[np.float64], target: float
) -> tuple[float, ...]:
    """Return all isolated piecewise-linear roots of ``values == target``.

    A segment whose two adjacent nodes equal the target is rejected: its
    continuum of crossings cannot be represented by a finite root tuple.
    """

    roots: list[float] = []
    for index, value in enumerate(values):
        if index + 1 < values.size and value == target and values[index + 1] == target:
            raise ComparisonError(
                "resonance target is flat across adjacent profile nodes; " "crossing is ambiguous"
            )
        if value == target:
            roots.append(float(radius[index]))
        if index + 1 >= values.size:
            continue
        next_value = values[index + 1]
        crosses_target = (value < target < next_value) or (next_value < target < value)
        if not crosses_target:
            continue

        # Normalize before subtraction so ±tiny and ±max q values retain the
        # correct crossing fraction without underflow or overflow.
        scale = max(abs(float(value)), abs(float(next_value)), abs(target))
        if scale == 0.0 or not np.isfinite(scale):
            raise ComparisonError("resonance crossing cannot be represented")
        value_scaled = float(value) / scale
        next_value_scaled = float(next_value) / scale
        target_scaled = target / scale
        denominator = next_value_scaled - value_scaled
        fraction = (target_scaled - value_scaled) / denominator
        if not np.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
            raise ComparisonError("resonance crossing cannot be represented")
        radius_lower = float(radius[index])
        radius_upper = float(radius[index + 1])
        root = radius_lower * (1.0 - fraction) + radius_upper * fraction
        if not np.isfinite(root):
            raise ComparisonError("resonance crossing cannot be represented")
        roots.append(float(root))
    roots.sort()
    unique: list[float] = []
    for root in roots:
        # Segment roots never duplicate an exact node (segments touching a
        # node root are skipped), so exact equality is sufficient here and
        # preserves genuinely distinct, very closely spaced crossings.
        if not unique or root != unique[-1]:
            unique.append(root)
    return tuple(unique)


def _covered(roots: tuple[float, ...], lower: float, upper: float) -> bool:
    return any(lower <= root <= upper for root in roots)


def _warnings(
    exclusions: ComparisonExclusions, resonance: ResonanceComparison | None
) -> tuple[str, ...]:
    warnings: list[str] = []
    profiles = (exclusions.oracle, exclusions.prepared)
    if any(
        profile.outside_domain_points_cm or profile.outside_domain_intervals_cm
        for profile in profiles
    ):
        warnings.append("domain-exclusions")
    if any(
        profile.outside_shared_overlap_points_cm or profile.outside_shared_overlap_intervals_cm
        for profile in profiles
    ):
        warnings.append("shared-overlap-exclusions")
    if resonance is not None and (
        not resonance.reference_crossing_covered or not resonance.candidate_crossing_covered
    ):
        warnings.append("resonance-not-covered")
    return tuple(warnings)


def _immutable_array(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return a read-only array whose write flag cannot be re-enabled."""

    copied = np.asarray(values, dtype=np.float64).copy()
    immutable = np.frombuffer(copied.tobytes(), dtype=np.float64)
    immutable.setflags(write=False)
    return immutable


__all__ = [
    "ComparisonError",
    "ComparisonExclusions",
    "ComparisonMeasurements",
    "ComparisonResult",
    "ExclusionInterval",
    "ProfileExclusions",
    "ResonanceComparison",
    "compare_profiles",
]
