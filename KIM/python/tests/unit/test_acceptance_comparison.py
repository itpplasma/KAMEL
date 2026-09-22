"""RED contract tests for the pure KIM scientific profile comparator.

The comparator is deliberately exercised with analytic, in-memory arrays.  The
oracle is always the reference and the prepared profile is always the
candidate; swapping those meanings would make acceptance results physically
misleading.  This suite specifies the small API expected from
``kim.acceptance_comparison`` without involving HDF5, plotting, or I/O.
"""

from __future__ import annotations

import numpy as np
import pytest
from kim.acceptance_comparison import ComparisonError, compare_profiles


def _compare(
    oracle_radius_cm: np.ndarray | list[float],
    oracle_values: np.ndarray | list[float],
    prepared_radius_cm: np.ndarray | list[float],
    prepared_values: np.ndarray | list[float],
    *,
    domain_cm: tuple[float, float] = (1.0, 7.0),
    interpolation_direction: str = "prepared_to_oracle",
    method: str = "linear",
    relative_floor: float = 1.0,
    tolerances: dict[str, float] | None = None,
    resonance: tuple[int, int] | float | None = None,
):
    return compare_profiles(
        oracle_radius_cm=oracle_radius_cm,
        oracle_values=oracle_values,
        prepared_radius_cm=prepared_radius_cm,
        prepared_values=prepared_values,
        domain_cm=domain_cm,
        interpolation_direction=interpolation_direction,
        method=method,
        relative_floor=relative_floor,
        tolerances=tolerances,
        resonance=resonance,
    )


def _raw_call(**kwargs: object):
    values: dict[str, object] = {
        "oracle_radius_cm": np.array([0.0, 2.0, 5.0, 9.0]),
        "oracle_values": np.array([1.0, 5.0, 11.0, 19.0]),
        "prepared_radius_cm": np.array([1.0, 3.0, 6.0, 7.0]),
        "prepared_values": np.array([3.0, 7.0, 13.0, 15.0]),
        "domain_cm": (2.0, 7.0),
        "interpolation_direction": "prepared_to_oracle",
        "method": "linear",
        "relative_floor": 1.0,
    }
    values.update(kwargs)
    return compare_profiles(**values)


def test_reports_exact_absolute_and_relative_measurements_without_decisions() -> None:
    oracle_radius = np.array([1.0, 3.0, 7.0])
    oracle_values = np.array([2.0, 4.0, 10.0])
    prepared_values = np.array([3.0, 8.0, 7.0])

    result = _compare(
        oracle_radius,
        oracle_values,
        oracle_radius.copy(),
        prepared_values,
        domain_cm=(1.0, 7.0),
        relative_floor=5.0,
    )

    assert result.measurements.absolute_rms == pytest.approx(np.sqrt(26.0 / 3.0), rel=1e-13)
    assert result.measurements.absolute_max == pytest.approx(4.0, rel=1e-13)
    assert result.measurements.relative_rms == pytest.approx(np.sqrt(0.77 / 3.0), rel=1e-13)
    assert result.measurements.relative_max == pytest.approx(0.8, rel=1e-13)
    assert result.threshold_decisions is None
    assert result.overall_pass is None
    assert isinstance(result.warnings, tuple)


def test_relative_denominator_uses_absolute_oracle_reference_with_a_floor() -> None:
    result = _compare(
        np.array([1.0, 2.0, 3.0]),
        np.array([-2.0, 0.0, 4.0]),
        np.array([1.0, 2.0, 3.0]),
        np.array([-1.0, 2.0, 2.0]),
        domain_cm=(1.0, 3.0),
        relative_floor=3.0,
    )

    # Signed errors (prepared minus oracle) are [+1, +2, -2].  The denominators are [3, 3, 4],
    # i.e. max(abs(oracle), floor), including negative and zero oracle values.
    assert result.measurements.relative_rms == pytest.approx(
        np.sqrt(((1.0 / 3.0) ** 2 + (2.0 / 3.0) ** 2 + 0.5**2) / 3.0), rel=1e-13
    )
    assert result.measurements.relative_max == pytest.approx(2.0 / 3.0, rel=1e-13)


@pytest.mark.parametrize(
    ("direction", "expected_radius", "expected_absolute_rms", "expected_relative_rms"),
    [
        (
            "prepared_to_oracle",
            (2.0, 5.0),
            np.sqrt(2.5),
            np.sqrt(((1.0 / 4.0) ** 2 + (2.0 / 25.0) ** 2) / 2.0),
        ),
        (
            "oracle_to_prepared",
            (3.0, 6.0, 7.0),
            2.0,
            np.sqrt(((2.0 / 11.0) ** 2 + (2.0 / 38.0) ** 2 + (2.0 / 51.0) ** 2) / 3.0),
        ),
    ],
)
def test_nonlinear_profiles_keep_oracle_as_reference_in_both_directions(
    direction: str,
    expected_radius: tuple[float, ...],
    expected_absolute_rms: float,
    expected_relative_rms: float,
) -> None:
    oracle_radius = np.array([0.0, 2.0, 5.0, 8.0])
    prepared_radius = np.array([1.0, 3.0, 6.0, 7.0])
    # A quadratic profile is intentionally not reproduced exactly by linear
    # interpolation.  This makes direction and the oracle denominator visible.
    oracle_values = oracle_radius**2
    prepared_values = prepared_radius**2

    result = _compare(
        oracle_radius,
        oracle_values,
        prepared_radius,
        prepared_values,
        domain_cm=(2.0, 7.0),
        interpolation_direction=direction,
        relative_floor=1.0,
    )

    np.testing.assert_array_equal(result.comparison_radius_cm, expected_radius)
    assert result.measurements.absolute_rms == pytest.approx(expected_absolute_rms, rel=1e-13)
    assert result.measurements.relative_rms == pytest.approx(expected_relative_rms, rel=1e-13)


def test_oracle_is_reference_and_signed_q_resonance_uses_q_equals_negative_m_over_n() -> None:
    oracle_radius = np.array([0.0, 2.0, 4.0, 8.0])
    prepared_radius = np.array([1.0, 3.0, 6.0, 7.0])
    # Monotonic analytic q(r) = 0.5*r - 3 crosses the signed target -3/2 at r=3.
    oracle_q = 0.5 * oracle_radius - 3.0
    prepared_q = 0.5 * prepared_radius - 3.0

    result = _compare(
        oracle_radius,
        oracle_q,
        prepared_radius,
        prepared_q,
        domain_cm=(2.0, 4.0),
        resonance=(3, 2),
    )

    np.testing.assert_array_equal(result.comparison_radius_cm, [2.0, 4.0])
    assert result.resonance.target_q == pytest.approx(-1.5)
    assert result.resonance.reference_crossing_covered is True
    assert result.resonance.candidate_crossing_covered is True
    assert result.resonance.reference_crossing_radii_cm == pytest.approx((3.0,))
    assert result.resonance.candidate_crossing_radii_cm == pytest.approx((3.0,))


def test_accepts_an_explicit_signed_resonance_target_and_reports_uncovered_crossing() -> None:
    oracle_radius = np.array([0.0, 2.0, 4.0, 8.0])
    prepared_radius = np.array([1.0, 3.0, 6.0, 7.0])
    oracle_q = 0.5 * oracle_radius - 3.0
    prepared_q = 0.5 * prepared_radius - 3.0

    result = _compare(
        oracle_radius,
        oracle_q,
        prepared_radius,
        prepared_q,
        domain_cm=(4.0, 7.0),
        resonance=-1.5,
    )

    assert result.resonance.target_q == pytest.approx(-1.5)
    assert result.resonance.reference_crossing_covered is False
    assert result.resonance.candidate_crossing_covered is False
    # The crossing is reported even when it falls outside the requested domain.
    assert result.resonance.reference_crossing_radii_cm == pytest.approx((3.0,))
    assert result.resonance.candidate_crossing_radii_cm == pytest.approx((3.0,))


def test_resonance_includes_profile_and_domain_boundary_crossings_and_all_multiple_crossings() -> (
    None
):
    radius = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
    # Signed target -1.5 occurs exactly at r=0, and then at r=3, 5, and 7.
    q = np.array([-1.5, -0.5, -2.5, -0.5, -2.5])

    result = _compare(
        radius,
        q,
        radius.copy(),
        q.copy(),
        domain_cm=(0.0, 5.0),
        resonance=-1.5,
    )

    assert result.resonance.reference_crossing_radii_cm == pytest.approx((0.0, 3.0, 5.0, 7.0))
    assert result.resonance.candidate_crossing_radii_cm == pytest.approx((0.0, 3.0, 5.0, 7.0))
    # At least one crossing in the inclusive requested domain and source
    # support means covered; crossings outside remain reported above.
    assert result.resonance.reference_crossing_covered is True
    assert result.resonance.candidate_crossing_covered is True


@pytest.mark.parametrize(
    ("direction", "expected_radius"),
    [
        ("prepared_to_oracle", (2.0, 5.0)),
        # Both endpoints of the requested domain are included when present in
        # the target grid; no target point is synthesized at an endpoint.
        ("oracle_to_prepared", (3.0, 6.0, 7.0)),
    ],
)
def test_interpolates_only_in_requested_direction_on_irregular_grids_and_includes_endpoints(
    direction: str, expected_radius: tuple[float, ...]
) -> None:
    oracle_radius = np.array([0.0, 2.0, 5.0, 9.0])
    prepared_radius = np.array([1.0, 3.0, 6.0, 7.0])
    oracle_values = 1.0 + 2.0 * oracle_radius
    prepared_values = 1.0 + 2.0 * prepared_radius

    result = _compare(
        oracle_radius,
        oracle_values,
        prepared_radius,
        prepared_values,
        domain_cm=(2.0, 7.0),
        interpolation_direction=direction,
    )

    np.testing.assert_array_equal(result.comparison_radius_cm, expected_radius)
    assert result.measurements.absolute_max == pytest.approx(0.0, abs=1e-14)
    assert result.measurements.relative_max == pytest.approx(0.0, abs=1e-14)


@pytest.mark.parametrize(
    ("direction", "expected_radius"),
    [
        ("prepared_to_oracle", (2.0, 5.0)),
        ("oracle_to_prepared", (1.0, 3.0, 7.0)),
    ],
)
def test_never_extrapolates_in_either_direction_and_reports_shared_overlap_edges(
    direction: str, expected_radius: tuple[float, ...]
) -> None:
    oracle_radius = np.array([0.0, 2.0, 5.0, 9.0])
    prepared_radius = np.array([1.0, 3.0, 7.0])
    oracle_values = 1.0 + 2.0 * oracle_radius
    prepared_values = 1.0 + 2.0 * prepared_radius

    result = _compare(
        oracle_radius,
        oracle_values,
        prepared_radius,
        prepared_values,
        domain_cm=(0.0, 9.0),
        interpolation_direction=direction,
    )

    # Target nodes outside the mutual continuous support are not extrapolated
    # in either direction and remain visible as exclusions.
    np.testing.assert_array_equal(result.comparison_radius_cm, expected_radius)
    assert result.exclusions.oracle.outside_domain_points_cm == ()
    assert result.exclusions.prepared.outside_domain_points_cm == ()
    assert result.exclusions.oracle.outside_shared_overlap_points_cm == pytest.approx((0.0, 9.0))
    assert result.exclusions.prepared.outside_shared_overlap_points_cm == ()
    assert result.exclusions.oracle.outside_shared_overlap_intervals_cm == pytest.approx(
        ((0.0, 1.0), (7.0, 9.0))
    )
    assert result.exclusions.prepared.outside_shared_overlap_intervals_cm == ()


def test_reports_requested_domain_exclusions_separately_from_overlap_exclusions() -> None:
    oracle_radius = np.array([0.0, 2.0, 5.0, 8.0, 9.0])
    prepared_radius = np.array([1.0, 3.0, 6.0, 7.0])

    result = _compare(
        oracle_radius,
        1.0 + 2.0 * oracle_radius,
        prepared_radius,
        1.0 + 2.0 * prepared_radius,
        domain_cm=(2.0, 6.0),
    )

    assert result.exclusions.oracle.outside_domain_points_cm == pytest.approx((0.0, 8.0, 9.0))
    assert result.exclusions.prepared.outside_domain_points_cm == pytest.approx((1.0, 7.0))
    assert result.exclusions.oracle.outside_domain_intervals_cm == pytest.approx(
        ((0.0, 2.0), (6.0, 9.0))
    )
    assert result.exclusions.prepared.outside_domain_intervals_cm == pytest.approx(
        ((1.0, 2.0), (6.0, 7.0))
    )
    # The requested-domain and shared-overlap categories are disjoint.  Every
    # point inside [2, 6] is supported by both continuous input grids, so no
    # point or interval is additionally excluded by overlap.
    assert result.exclusions.oracle.outside_shared_overlap_points_cm == ()
    assert result.exclusions.prepared.outside_shared_overlap_points_cm == ()
    assert result.exclusions.oracle.outside_shared_overlap_intervals_cm == ()
    assert result.exclusions.prepared.outside_shared_overlap_intervals_cm == ()


def test_supplied_tolerances_produce_per_metric_decisions_and_an_overall_result() -> None:
    tolerances = {
        "absolute_rms": 3.0,
        "absolute_max": 4.0,
        "relative_rms": 0.6,
        "relative_max": 0.8,
    }
    result = _compare(
        np.array([1.0, 3.0, 7.0]),
        np.array([2.0, 4.0, 10.0]),
        np.array([1.0, 3.0, 7.0]),
        np.array([3.0, 8.0, 7.0]),
        domain_cm=(1.0, 7.0),
        relative_floor=5.0,
        tolerances=tolerances,
    )

    assert result.threshold_decisions == {
        "absolute_rms": True,
        "absolute_max": True,
        "relative_rms": True,
        "relative_max": True,
    }
    assert result.overall_pass is True

    failing = _compare(
        np.array([1.0, 3.0, 7.0]),
        np.array([2.0, 4.0, 10.0]),
        np.array([1.0, 3.0, 7.0]),
        np.array([3.0, 8.0, 7.0]),
        domain_cm=(1.0, 7.0),
        relative_floor=5.0,
        tolerances={**tolerances, "absolute_max": 3.9},
    )
    assert failing.threshold_decisions["absolute_max"] is False
    assert failing.overall_pass is False


def test_nonempty_partial_tolerances_decide_only_supplied_metrics() -> None:
    result = _compare(
        np.array([1.0, 3.0, 7.0]),
        np.array([2.0, 4.0, 10.0]),
        np.array([1.0, 3.0, 7.0]),
        np.array([3.0, 8.0, 7.0]),
        domain_cm=(1.0, 7.0),
        relative_floor=5.0,
        tolerances={"absolute_max": 4.0},
    )

    assert result.threshold_decisions == {"absolute_max": True}
    assert result.overall_pass is True

    failing = _compare(
        np.array([1.0, 3.0, 7.0]),
        np.array([2.0, 4.0, 10.0]),
        np.array([1.0, 3.0, 7.0]),
        np.array([3.0, 8.0, 7.0]),
        domain_cm=(1.0, 7.0),
        relative_floor=5.0,
        tolerances={"relative_rms": 0.5},
    )
    assert failing.threshold_decisions == {"relative_rms": False}
    assert failing.overall_pass is False


@pytest.mark.parametrize(
    "bad_tolerances",
    [
        {},
        {"unknown": 1.0},
        {"absolute_rms": -1.0, "absolute_max": 1.0, "relative_rms": 1.0, "relative_max": 1.0},
        {"absolute_rms": np.nan, "absolute_max": 1.0, "relative_rms": 1.0, "relative_max": 1.0},
        {"absolute_rms": 1.0, "absolute_max": np.inf, "relative_rms": 1.0, "relative_max": 1.0},
    ],
)
def test_rejects_empty_unknown_negative_or_nonfinite_tolerances(
    bad_tolerances: dict[str, float],
) -> None:
    with pytest.raises(ComparisonError, match="tolerance"):
        _compare(
            np.array([1.0, 3.0, 7.0]),
            np.array([2.0, 4.0, 10.0]),
            np.array([1.0, 3.0, 7.0]),
            np.array([3.0, 8.0, 7.0]),
            domain_cm=(1.0, 7.0),
            relative_floor=5.0,
            tolerances=bad_tolerances,
        )


@pytest.mark.parametrize("bad_floor", [0.0, -1.0, np.nan, np.inf, -np.inf])
def test_rejects_non_positive_or_nonfinite_relative_floor(bad_floor: float) -> None:
    with pytest.raises(ComparisonError, match="floor"):
        _compare(
            np.array([0.0, 1.0, 2.0]),
            np.array([1.0, 2.0, 3.0]),
            np.array([0.0, 1.0, 2.0]),
            np.array([1.0, 2.0, 3.0]),
            domain_cm=(0.0, 2.0),
            relative_floor=bad_floor,
        )


@pytest.mark.parametrize(
    ("bad_domain", "message"),
    [
        ((2.0, 1.0), "domain"),
        ((1.0, 1.0), "domain"),
        ((np.nan, 2.0), "domain"),
        ((1.0, np.inf), "domain"),
        ((10.0, 20.0), "overlap"),
    ],
)
def test_rejects_reversed_empty_nonfinite_or_non_overlapping_domains(
    bad_domain: tuple[float, float], message: str
) -> None:
    with pytest.raises(ComparisonError, match=message):
        _compare(
            np.array([0.0, 2.0, 5.0, 9.0]),
            np.array([1.0, 5.0, 11.0, 19.0]),
            np.array([1.0, 3.0, 6.0, 7.0]),
            np.array([3.0, 7.0, 13.0, 15.0]),
            domain_cm=bad_domain,
        )


def test_rejects_omitted_comparison_domain() -> None:
    with pytest.raises((ComparisonError, TypeError), match="domain"):
        compare_profiles(
            oracle_radius_cm=np.array([0.0, 1.0, 2.0]),
            oracle_values=np.array([1.0, 2.0, 3.0]),
            prepared_radius_cm=np.array([0.0, 1.0, 2.0]),
            prepared_values=np.array([1.0, 2.0, 3.0]),
            interpolation_direction="prepared_to_oracle",
            method="linear",
            relative_floor=1.0,
        )


@pytest.mark.parametrize("bad_resonance", [(0, 2), (2, 0), np.nan, np.inf, -np.inf])
def test_rejects_zero_mode_numbers_and_nonfinite_signed_resonance_targets(
    bad_resonance: tuple[int, int] | float,
) -> None:
    with pytest.raises(ComparisonError, match="resonance|mode|target"):
        _compare(
            np.array([0.0, 1.0, 2.0]),
            np.array([-2.0, -1.0, 0.0]),
            np.array([0.0, 1.0, 2.0]),
            np.array([-2.0, -1.0, 0.0]),
            domain_cm=(0.0, 2.0),
            resonance=bad_resonance,
        )


@pytest.mark.parametrize(
    ("bad_arrays", "message"),
    [
        (
            {
                "oracle_radius_cm": np.array([0.0, 1.0, 1.0]),
                "oracle_values": np.array([1.0, 2.0, 3.0]),
            },
            "increasing|duplicate",
        ),
        (
            {
                "prepared_radius_cm": np.array([0.0, 1.0, 1.0]),
                "prepared_values": np.array([1.0, 2.0, 3.0]),
            },
            "increasing|duplicate",
        ),
        (
            {"oracle_radius_cm": np.array([0.0]), "oracle_values": np.array([1.0])},
            "at least|length|sample",
        ),
        (
            {"prepared_radius_cm": np.array([0.0]), "prepared_values": np.array([1.0])},
            "at least|length|sample",
        ),
    ],
)
def test_rejects_duplicate_or_too_short_input_arrays(
    bad_arrays: dict[str, np.ndarray], message: str
) -> None:
    with pytest.raises(ComparisonError, match=message):
        _raw_call(**bad_arrays)


@pytest.mark.parametrize("direction", ["prepared_to_oracle", "oracle_to_prepared"])
def test_rejects_continuous_overlap_without_a_target_evaluation_point(direction: str) -> None:
    with pytest.raises(ComparisonError, match="point|evaluation|comparison"):
        _compare(
            np.array([0.0, 2.0, 4.0]),
            np.array([0.0, 4.0, 16.0]),
            np.array([0.0, 1.0, 3.0, 4.0]),
            np.array([0.0, 1.0, 9.0, 16.0]),
            domain_cm=(1.2, 1.8),
            interpolation_direction=direction,
        )


@pytest.mark.parametrize(
    ("direction", "method"),
    [
        ("implicit", "linear"),
        ("prepared_to_oracle", "implicit"),
        ("prepared_to_candidate", "linear"),
        ("prepared_to_oracle", "cubic"),
    ],
)
def test_rejects_implicit_or_unsupported_interpolation_direction_and_method(
    direction: str, method: str
) -> None:
    with pytest.raises((ComparisonError, TypeError), match="(direction|method|interpol)"):
        if direction == "implicit":
            # Omitted direction is intentionally invalid; the explicit method
            # remains present so this probes only direction handling.
            compare_profiles(
                oracle_radius_cm=np.array([0.0, 1.0, 2.0]),
                oracle_values=np.array([1.0, 2.0, 3.0]),
                prepared_radius_cm=np.array([0.0, 1.0, 2.0]),
                prepared_values=np.array([1.0, 2.0, 3.0]),
                domain_cm=(0.0, 2.0),
                method=method,
                relative_floor=1.0,
            )
        else:
            _compare(
                np.array([0.0, 1.0, 2.0]),
                np.array([1.0, 2.0, 3.0]),
                np.array([0.0, 1.0, 2.0]),
                np.array([1.0, 2.0, 3.0]),
                domain_cm=(0.0, 2.0),
                interpolation_direction=direction,
                method=method,
            )


def test_rejects_implicit_method_even_when_direction_is_explicit() -> None:
    with pytest.raises((ComparisonError, TypeError), match="(method|interpol)"):
        compare_profiles(
            oracle_radius_cm=np.array([0.0, 1.0, 2.0]),
            oracle_values=np.array([1.0, 2.0, 3.0]),
            prepared_radius_cm=np.array([0.0, 1.0, 2.0]),
            prepared_values=np.array([1.0, 2.0, 3.0]),
            domain_cm=(0.0, 2.0),
            interpolation_direction="prepared_to_oracle",
            relative_floor=1.0,
        )


@pytest.mark.parametrize(
    "bad_arrays",
    [
        {
            "oracle_radius_cm": np.array([0.0, np.nan, 2.0]),
        },
        {
            "oracle_values": np.array([1.0, np.inf, 3.0]),
        },
        {
            "prepared_radius_cm": np.array([0.0, 2.0, 1.0]),
        },
        {
            "oracle_values": np.array([1.0, 2.0]),
        },
        {
            "prepared_values": np.array([[1.0, 2.0, 3.0]]),
        },
    ],
)
def test_rejects_nonfinite_unsorted_nonmatching_or_non_1d_arrays(
    bad_arrays: dict[str, np.ndarray],
) -> None:
    with pytest.raises(ComparisonError, match="(finite|increasing|length|one-dimensional|shape)"):
        _raw_call(**bad_arrays)


def test_does_not_mutate_any_input_array() -> None:
    oracle_radius = np.array([0.0, 2.0, 5.0, 9.0])
    oracle_values = np.array([1.0, 5.0, 11.0, 19.0])
    prepared_radius = np.array([1.0, 3.0, 6.0, 7.0])
    prepared_values = np.array([3.0, 7.0, 13.0, 15.0])
    original = tuple(
        array.copy() for array in (oracle_radius, oracle_values, prepared_radius, prepared_values)
    )

    _compare(oracle_radius, oracle_values, prepared_radius, prepared_values, domain_cm=(2.0, 7.0))

    for actual, expected in zip(
        (oracle_radius, oracle_values, prepared_radius, prepared_values), original, strict=True
    ):
        np.testing.assert_array_equal(actual, expected)
