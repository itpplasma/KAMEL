from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest
from kim import (
    ConventionError,
    SourceMetadata,
    convert_from_kim,
    convert_inputs,
    convert_quantity,
    reconstruct_fourier,
    reconstruct_fourier_derivative,
    resonance_target,
)
from pydantic import ValidationError


def metadata(**updates: object) -> SourceMetadata:
    values: dict[str, object] = {
        "source": "approved-test-source",
        "coordinate": "r_eff",
        "radius_unit": "m",
        "density_unit": "m^-3",
        "temperature_unit": "eV",
        "magnetic_field_unit": "T",
        "magnetic_field_sign": "preserve",
        "electric_field_unit": "V/m",
        "velocity_unit": "m/s",
        "frequency_unit": "rad/s",
        "frequency_convention": "signed_omega",
        "time_phase": "exp(-i omega t)",
        "perturbation_phase": "preserve",
        "equilibrium_provenance": "approved-test-equilibrium",
        "q_unit": "1",
        "mode_convention": "signed",
        "fourier_phase": "exp(+i k r)",
        "complex_phase": "preserve",
        "resonance_convention": "q=-m/n",
    }
    values.update(updates)
    return SourceMetadata(**values)


def test_source_metadata_is_strict_and_serializable() -> None:
    source = metadata()

    assert source.model_dump(mode="json")["frequency_unit"] == "rad/s"
    missing = source.model_dump()
    missing.pop("fourier_phase")
    with pytest.raises(ValidationError, match="fourier_phase"):
        SourceMetadata.model_validate(missing)
    with pytest.raises(ValidationError, match="frequency_unit"):
        metadata(frequency_unit="Hz")
    with pytest.raises(ValidationError, match="temperature_unit"):
        metadata(temperature_unit="K")
    with pytest.raises(ValidationError, match="time_phase"):
        metadata(time_phase="exp(+i omega t)")
    with pytest.raises(ValidationError, match="extra"):
        SourceMetadata.model_validate(metadata().model_dump() | {"unknown": "value"})


@pytest.mark.parametrize(
    ("quantity", "source", "expected"),
    [
        ("radius", 1.25, 125.0),
        ("density", 2.5e19, 2.5e13),
        ("magnetic_field", -2.0, -2.0e4),
        ("electric_field", -29979.2458, -1.0),
        ("toroidal_velocity", -4.0, -400.0),
        ("temperature", 12.0, 12.0),
        ("frequency", -7.0, -7.0),
    ],
)
def test_approved_unit_conversions_have_independent_reference_values(
    quantity: str, source: float, expected: float
) -> None:
    assert convert_quantity(source, quantity, metadata()) == pytest.approx(expected)


def test_conversion_report_records_operations_and_signed_complex_values() -> None:
    report_time = datetime(2026, 9, 14, 12, tzinfo=timezone.utc)
    result = convert_inputs(
        {
            "radius": np.array([0.0, 0.5]),
            "n_e": np.array([1.0e19, 2.0e19]),
            "Er": np.array([-29979.2458 + 1.0j]),
            "m_mode": np.array([-6]),
            "n_mode": np.array([2]),
        },
        metadata(),
        created_at=report_time,
    )

    np.testing.assert_allclose(result.values["radius"], [0.0, 50.0])
    np.testing.assert_allclose(result.values["n_e"], [1.0e13, 2.0e13])
    assert result.values["Er"][0] == pytest.approx(-1.0 + 1.0j / 29979.2458)
    assert result.values["m_mode"].tolist() == [-6]
    assert [operation.quantity for operation in result.report.operations] == [
        "radius",
        "n_e",
        "Er",
    ]
    assert result.report.created_at == report_time
    assert result.report.model_dump(mode="json")["source"]["source"] == "approved-test-source"
    mode_result = convert_inputs({"n": -2}, metadata())
    assert mode_result.values["n"].item() == -2


def test_conversions_round_trip_for_invertible_linear_units() -> None:
    source = metadata()
    converted = {
        "radius": convert_quantity([0.0, 1.2], "radius", source),
        "n_e": convert_quantity([1.0e19, 2.0e19], "density", source),
        "B": convert_quantity([-1.5], "magnetic_field", source),
        "Er": convert_quantity([-2.0], "electric_field", source),
        "Vz": convert_quantity([3.0], "toroidal_velocity", source),
    }
    target = metadata(
        radius_unit="cm",
        density_unit="1/cm^3",
        magnetic_field_unit="G",
        electric_field_unit="statV/cm",
        velocity_unit="cm/s",
    )
    np.testing.assert_allclose(convert_from_kim(converted["radius"], "radius", source), [0.0, 1.2])
    np.testing.assert_allclose(convert_from_kim(converted["n_e"], "density", source), [1e19, 2e19])
    np.testing.assert_allclose(convert_from_kim(converted["B"], "magnetic_field", source), [-1.5])
    np.testing.assert_allclose(convert_from_kim(converted["Er"], "electric_field", source), [-2.0])
    np.testing.assert_allclose(
        convert_from_kim(converted["Vz"], "toroidal_velocity", source), [3.0]
    )


def test_nonfinite_values_and_dimension_mismatch_fail_closed() -> None:
    with pytest.raises(ConventionError, match="density.*finite"):
        convert_quantity([1.0, np.nan], "density", metadata())
    with pytest.raises(ConventionError, match="unsupported quantity"):
        convert_quantity([1.0], "pressure", metadata())
    with pytest.raises(ConventionError, match="unsupported input quantity"):
        convert_inputs({"pressure": [1.0]}, metadata())
    with pytest.raises(ConventionError, match="signed integer"):
        convert_inputs({"n": 2.0}, metadata())


def test_phase_vector_preserves_signed_modes_and_complex_phase() -> None:
    values = reconstruct_fourier(
        [0.0, np.pi],
        [-1, 0, 1],
        [1.0 + 2.0j, 3.0 - 1.0j, -2.0 + 0.5j],
        metadata(radius_unit="cm"),
        period=2.0 * np.pi,
    )

    np.testing.assert_allclose(values, [2.0 + 1.5j, 4.0 - 3.5j])
    derivative = reconstruct_fourier_derivative(
        [0.0],
        [-1, 0, 1],
        [1.0 + 2.0j, 3.0 - 1.0j, -2.0 + 0.5j],
        metadata(radius_unit="cm"),
        period=2.0 * np.pi,
    )
    np.testing.assert_allclose(derivative, [1.5 - 3.0j])


def test_zero_mode_is_constant_for_any_period() -> None:
    values = reconstruct_fourier(
        [0.0, 10.0, 25.0],
        [0],
        [2.0 + 3.0j],
        metadata(radius_unit="cm"),
        period=7.0,
    )
    np.testing.assert_allclose(values, [2.0 + 3.0j] * 3)


def test_fourier_rejects_float_mode_arrays_even_when_integral_valued() -> None:
    with pytest.raises(ConventionError, match="signed integer"):
        reconstruct_fourier(
            [0.0],
            [0.0],
            [2.0 + 3.0j],
            metadata(radius_unit="cm"),
            period=7.0,
        )


@pytest.mark.parametrize(
    ("m_mode", "n_mode", "expected"),
    [(-7, 2, 3.5), (7, -2, 3.5), (-7, -2, -3.5)],
)
def test_resonance_target_preserves_signed_fortran_convention(
    m_mode: int, n_mode: int, expected: float
) -> None:
    assert resonance_target(m_mode, n_mode, metadata()) == expected
