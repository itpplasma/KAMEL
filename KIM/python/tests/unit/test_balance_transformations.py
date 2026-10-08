"""TDD contract for the pure BALANCE-to-KIM transformations.

These tests deliberately target the narrow ``kim.balance_adoption`` boundary.
The adapter owns only the approved BALANCE conventions; it must not become a
generic arbitrary-unit conversion layer or change the existing MARS-F path.
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pytest
from kim import ConventionError
from kim.balance_adoption import (
    BalanceTransformation,
    convert_balance_density,
    convert_balance_q,
    convert_balance_rotation,
    convert_balance_temperature,
)
from kim.conventions import ConversionOperation

R0_CM = 165.0


def _assert_operation(
    result: BalanceTransformation,
    *,
    quantity: str,
    source_unit: str,
    target_unit: str,
    factor: float,
    operation: str,
) -> None:
    record = result.operation
    assert type(record) is ConversionOperation
    assert record.quantity == quantity
    assert record.source_unit == source_unit
    assert record.target_unit == target_unit
    assert record.factor == pytest.approx(factor)
    assert record.operation == operation


@pytest.mark.parametrize(
    ("function", "values", "kwargs", "missing"),
    [
        (
            convert_balance_density,
            [1.0e19, 2.0e19],
            {"coordinate_unit": "1", "source_unit": "1/m^3"},
            "coordinate",
        ),
        (
            convert_balance_density,
            [1.0e19, 2.0e19],
            {"coordinate": "rho_pol", "source_unit": "1/m^3"},
            "coordinate_unit",
        ),
        (
            convert_balance_density,
            [1.0e19, 2.0e19],
            {"coordinate": "rho_pol", "coordinate_unit": "1"},
            "source_unit",
        ),
        (
            convert_balance_temperature,
            [100.0, 80.0],
            {"coordinate": "rho_pol", "coordinate_unit": "1", "source_unit": "eV"},
            "species",
        ),
        (
            convert_balance_q,
            [1.0, 2.0],
            {
                "coordinate": "rho_pol",
                "coordinate_unit": "1",
                "operation": "preserve",
            },
            "source_unit",
        ),
    ],
)
def test_transformations_require_explicit_metadata(
    function: object, values: list[float], kwargs: dict[str, object], missing: str
) -> None:
    """The API must not fill scientific metadata from quantity names or values."""

    with pytest.raises(TypeError, match=missing):
        function([0.0, 1.0], values, **kwargs)  # type: ignore[operator]


def test_density_conversion_uses_explicit_si_to_cgs_factor_and_preserves_grid() -> None:
    rho_pol = np.array([0.0, 0.5, 1.0], dtype=np.float64)
    density = np.array([1.0e19, 2.0e19, 3.0e19], dtype=np.float64)
    rho_original = rho_pol.copy()
    density_original = density.copy()

    result = convert_balance_density(
        rho_pol,
        density,
        coordinate="rho_pol",
        coordinate_unit="1",
        source_unit="1/m^3",
    )

    np.testing.assert_allclose(result.coordinate, rho_pol)
    np.testing.assert_array_equal(result.values, [1.0e13, 2.0e13, 3.0e13])
    assert len(result.coordinate) == len(rho_pol)
    assert len(result.values) == len(density)
    assert not np.shares_memory(result.coordinate, rho_pol)
    assert not np.shares_memory(result.values, density)
    np.testing.assert_array_equal(rho_pol, rho_original)
    np.testing.assert_array_equal(density, density_original)
    _assert_operation(
        result,
        quantity="density",
        source_unit="1/m^3",
        target_unit="1/cm^3",
        factor=1.0e-6,
        operation="density_si_to_cgs",
    )


@pytest.mark.parametrize(
    ("species", "quantity"),
    [
        ("electron", "electron_temperature"),
        ("ion", "ion_temperature"),
    ],
)
def test_temperature_conversion_is_explicit_eV_identity_with_detached_outputs(
    species: str, quantity: str
) -> None:
    rho_pol = np.array([0.0, 0.4, 1.0], dtype=np.float64)
    temperature = np.array([100.0, 80.0, 40.0], dtype=np.float64)
    rho_original = rho_pol.copy()
    temperature_original = temperature.copy()

    result = convert_balance_temperature(
        rho_pol,
        temperature,
        species=species,
        coordinate="rho_pol",
        coordinate_unit="1",
        source_unit="eV",
    )

    np.testing.assert_array_equal(result.coordinate, rho_pol)
    np.testing.assert_array_equal(result.values, temperature)
    assert result.coordinate is not rho_pol
    assert result.values is not temperature
    assert not np.shares_memory(result.coordinate, rho_pol)
    assert not np.shares_memory(result.values, temperature)
    np.testing.assert_array_equal(rho_pol, rho_original)
    np.testing.assert_array_equal(temperature, temperature_original)
    _assert_operation(
        result,
        quantity=quantity,
        source_unit="eV",
        target_unit="eV",
        factor=1.0,
        operation="preserve",
    )


def test_angular_rotation_becomes_signed_toroidal_velocity_using_r0() -> None:
    rho_pol = np.array([0.0, 0.5, 1.0], dtype=np.float64)
    omega = np.array([-2.0, 0.0, 3.0], dtype=np.float64)
    omega_original = omega.copy()

    result = convert_balance_rotation(
        rho_pol,
        omega,
        coordinate="rho_pol",
        coordinate_unit="1",
        source_unit="rad/s",
        r0_cm=R0_CM,
    )

    np.testing.assert_array_equal(result.coordinate, rho_pol)
    np.testing.assert_array_equal(result.values, [-330.0, 0.0, 495.0])
    assert len(result.coordinate) == len(rho_pol)
    assert len(result.values) == len(omega)
    assert not np.shares_memory(result.coordinate, rho_pol)
    assert not np.shares_memory(result.values, omega)
    np.testing.assert_array_equal(omega, omega_original)
    _assert_operation(
        result,
        quantity="toroidal_velocity",
        source_unit="rad/s",
        target_unit="cm/s",
        factor=R0_CM,
        operation="omega_to_v_phi",
    )
    assert result.operation.model_dump(mode="json")["parameters"]["major_radius_cm"] == R0_CM


@pytest.mark.parametrize(
    "r0_cm",
    [None, np.nan, 0.0, -1.0, 10**1000],
    ids=["missing", "nonfinite", "zero", "negative", "overflow"],
)
def test_rotation_rejects_missing_nonfinite_or_nonpositive_r0(r0_cm: float | None) -> None:
    with pytest.raises(ConventionError, match="(?i)(radius|r0)"):
        convert_balance_rotation(
            [0.0, 1.0],
            [1.0, 2.0],
            coordinate="rho_pol",
            coordinate_unit="1",
            source_unit="rad/s",
            r0_cm=r0_cm,
        )


@pytest.mark.parametrize(
    ("operation", "expected", "factor"),
    [
        ("preserve", [1.0, -3.5, 2.0], 1.0),
        ("negate", [-1.0, 3.5, -2.0], -1.0),
    ],
)
def test_q_operation_is_explicit_and_preserves_signed_values(
    operation: str, expected: list[float], factor: float
) -> None:
    rho_pol = np.array([0.0, 0.5, 1.0], dtype=np.float64)
    q = np.array([1.0, -3.5, 2.0], dtype=np.float64)
    q_original = q.copy()

    result = convert_balance_q(
        rho_pol,
        q,
        coordinate="rho_pol",
        coordinate_unit="1",
        source_unit="1",
        operation=operation,
    )

    np.testing.assert_array_equal(result.coordinate, rho_pol)
    np.testing.assert_array_equal(result.values, expected)
    assert not np.shares_memory(result.coordinate, rho_pol)
    assert not np.shares_memory(result.values, q)
    np.testing.assert_array_equal(q, q_original)
    _assert_operation(
        result,
        quantity="q",
        source_unit="1",
        target_unit="1",
        factor=factor,
        operation=operation,
    )


def test_q_rejects_an_operation_that_would_require_inference() -> None:
    with pytest.raises(ConventionError, match="(?i)q"):
        convert_balance_q(
            [0.0, 1.0],
            [1.0, 2.0],
            coordinate="rho_pol",
            coordinate_unit="1",
            source_unit="1",
            operation="infer",
        )


@pytest.mark.parametrize(
    ("function", "values", "kwargs", "message"),
    [
        (
            convert_balance_density,
            [1.0e19, 2.0e19],
            {"source_unit": "kg/m^3"},
            "unit",
        ),
        (
            convert_balance_temperature,
            [100.0, 80.0],
            {"species": "electron", "source_unit": "K"},
            "unit",
        ),
        (
            convert_balance_rotation,
            [1.0, 2.0],
            {"source_unit": "m/s", "r0_cm": R0_CM},
            "unit",
        ),
    ],
)
def test_rejects_unsupported_units_instead_of_inferring_a_conversion(
    function: object, values: list[float], kwargs: dict[str, object], message: str
) -> None:
    with pytest.raises(ConventionError, match=message):
        function(  # type: ignore[operator]
            [0.0, 1.0],
            values,
            coordinate="rho_pol",
            coordinate_unit="1",
            **kwargs,
        )


@pytest.mark.parametrize(
    "coordinate_metadata",
    [
        {"coordinate": "rho_pol", "coordinate_unit": "cm"},
        {"coordinate": "unknown_coordinate", "coordinate_unit": "1"},
    ],
)
def test_rejects_invalid_coordinate_metadata_instead_of_inference(
    coordinate_metadata: dict[str, str],
) -> None:
    with pytest.raises(ConventionError, match="coordinate"):
        convert_balance_density(
            [0.0, 1.0],
            [1.0e19, 2.0e19],
            source_unit="1/m^3",
            **coordinate_metadata,
        )


@pytest.mark.parametrize("coordinate", ["rho_pol", "sqrt_psiN"])
def test_accepts_explicit_balance_coordinate_spellings(coordinate: str) -> None:
    result = convert_balance_density(
        [0.0, 1.0],
        [1.0e19, 2.0e19],
        coordinate=coordinate,
        coordinate_unit="1",
        source_unit="1/m^3",
    )

    np.testing.assert_array_equal(result.coordinate, [0.0, 1.0])
    np.testing.assert_array_equal(result.values, [1.0e13, 2.0e13])


def test_rejects_coordinate_and_value_length_mismatch() -> None:
    with pytest.raises(ConventionError, match="(?i)(length|shape|size)"):
        convert_balance_density(
            [0.0, 0.5, 1.0],
            [1.0e19, 2.0e19],
            coordinate="rho_pol",
            coordinate_unit="1",
            source_unit="1/m^3",
        )


@pytest.mark.parametrize("coordinate_values", [[0.0, np.nan], [0.0, np.inf]])
def test_rejects_nonfinite_coordinate_values(coordinate_values: list[float]) -> None:
    with pytest.raises(ConventionError, match="(?i)(finite|coordinate)"):
        convert_balance_density(
            coordinate_values,
            [1.0e19, 2.0e19],
            coordinate="rho_pol",
            coordinate_unit="1",
            source_unit="1/m^3",
        )


def test_rejects_nonfinite_profile_values_before_conversion() -> None:
    with pytest.raises(ConventionError, match="finite"):
        convert_balance_density(
            [0.0, 1.0],
            [1.0e19, np.nan],
            coordinate="rho_pol",
            coordinate_unit="1",
            source_unit="1/m^3",
        )


@pytest.mark.parametrize("complex_position", ["coordinate", "values"])
def test_rejects_complex_profile_inputs_without_discarding_imaginary_parts(
    complex_position: str,
) -> None:
    coordinate = [0.0, 1.0]
    values = [1.0e19, 2.0e19]
    if complex_position == "coordinate":
        coordinate = [0.0 + 1.0j, 1.0 + 0.0j]  # type: ignore[assignment]
    else:
        values = [1.0e19 + 1.0j, 2.0e19 + 0.0j]  # type: ignore[assignment]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(ConventionError, match="(?i)complex"):
            convert_balance_density(
                coordinate,
                values,
                coordinate="rho_pol",
                coordinate_unit="1",
                source_unit="1/m^3",
            )
    assert not any("complex" in str(warning.message).lower() for warning in caught)


def test_transformation_output_arrays_are_read_only() -> None:
    result = convert_balance_density(
        [0.0, 1.0],
        [1.0e19, 2.0e19],
        coordinate="rho_pol",
        coordinate_unit="1",
        source_unit="1/m^3",
    )

    with pytest.raises(ValueError, match="read-only"):
        result.coordinate[0] = 0.25
    with pytest.raises(ValueError, match="read-only"):
        result.values[0] = 0.25


def test_conversion_operation_parameters_are_deeply_immutable_and_serializable() -> None:
    operation = ConversionOperation(
        quantity="example",
        source_unit="1",
        target_unit="1",
        factor=1.0,
        operation="preserve",
        parameters={"nested": {"items": [1, 2]}},
    )

    assert operation.model_dump(mode="json")["parameters"] == {"nested": {"items": [1, 2]}}
    assert '"parameters":{"nested":{"items":[1,2]}}' in operation.model_dump_json()
    with pytest.raises(TypeError, match="immutable"):
        operation.parameters["new"] = 3
    with pytest.raises(TypeError, match="immutable"):
        operation.parameters["nested"]["new"] = 3
    with pytest.raises(TypeError):
        operation.parameters["nested"]["items"][0] = 9


def test_conversion_operation_without_parameters_keeps_empty_default() -> None:
    operation = ConversionOperation(
        quantity="example",
        source_unit="1",
        target_unit="1",
        factor=1.0,
        operation="preserve",
    )

    assert operation.parameters == {}
    assert operation.model_dump(mode="json")["parameters"] == {}


def test_conversion_operation_copies_preserve_immutable_parameters() -> None:
    operation = ConversionOperation(
        quantity="example",
        source_unit="1",
        target_unit="1",
        factor=1.0,
        operation="preserve",
        parameters={"nested": {"items": [1, 2]}},
    )

    copies = [operation.model_copy(deep=True), copy.deepcopy(operation)]
    for copied in copies:
        assert copied.parameters is operation.parameters
        assert copied.model_dump(mode="json")["parameters"] == {"nested": {"items": [1, 2]}}
        with pytest.raises(TypeError, match="immutable"):
            copied.parameters["new"] = 3
