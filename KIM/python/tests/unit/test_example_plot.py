from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pytest
from kim import ResultError

_PLOT_MODULE_PATH = Path(__file__).parents[2] / "examples" / "plot_periodic.py"
_PLOT_MODULE_SPEC = importlib.util.spec_from_file_location(
    "plot_periodic_example", _PLOT_MODULE_PATH
)
assert _PLOT_MODULE_SPEC is not None and _PLOT_MODULE_SPEC.loader is not None
_PLOT_MODULE = importlib.util.module_from_spec(_PLOT_MODULE_SPEC)
_PLOT_MODULE_SPEC.loader.exec_module(_PLOT_MODULE)
plot_periodic_result = _PLOT_MODULE.plot_periodic_result


def _compound(values: np.ndarray) -> np.ndarray:
    encoded = np.empty(values.shape, dtype=[("real", "<f8"), ("imag", "<f8")])
    encoded["real"] = values.real
    encoded["imag"] = values.imag
    return encoded


def _periodic_file(path: Path) -> Path:
    radius = np.array([8.0, 9.0, 10.0, 11.0])
    potential = np.array([1.0 + 2.0j, 2.0 + 3.0j, 3.0 + 4.0j, 4.0 + 5.0j])
    current = np.array([0.5 + 1.0j, 1.0 + 1.5j, 1.5 + 2.0j, 2.0 + 2.5j])
    with h5py.File(path, "w") as handle:
        handle.create_dataset("backs/e/r", data=radius)
        handle.create_dataset("fields/Phi", data=_compound(potential))
        handle.create_dataset("fields/jpar", data=_compound(current))
        handle.create_dataset("setup/periodic_scale/dx_asis", data=0.6)
    return path


def test_plot_periodic_result_saves_and_labels_stored_fields(tmp_path: Path) -> None:
    result_path = _periodic_file(tmp_path / "result.h5")
    figure_path = tmp_path / "periodic.png"
    before = result_path.read_bytes()

    figure = plot_periodic_result(result_path, figure_path)

    assert figure_path.is_file()
    assert figure_path.stat().st_size > 0
    assert figure.axes[0].get_xlabel() == "r_eff [cm]"
    assert figure.axes[0].get_ylabel() == "Phi [statV]"
    assert figure.axes[1].get_xlabel() == "r_eff [cm]"
    assert figure.axes[1].get_ylabel() == "jpar [statA/cm^2]"

    potential_lines = {line.get_label(): line for line in figure.axes[0].get_lines()}
    current_lines = {line.get_label(): line for line in figure.axes[1].get_lines()}
    np.testing.assert_array_equal(potential_lines["Re(Phi)"].get_xdata(), [8.0, 9.0, 10.0, 11.0])
    np.testing.assert_array_equal(potential_lines["Re(Phi)"].get_ydata(), [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(potential_lines["Im(Phi)"].get_ydata(), [2.0, 3.0, 4.0, 5.0])
    np.testing.assert_array_equal(current_lines["Re(jpar)"].get_ydata(), [0.5, 1.0, 1.5, 2.0])
    np.testing.assert_array_equal(current_lines["Im(jpar)"].get_ydata(), [1.0, 1.5, 2.0, 2.5])
    assert potential_lines["resonance"].get_xdata()[0] == pytest.approx(10.0)
    assert current_lines["resonance"].get_xdata()[0] == pytest.approx(10.0)
    assert any(patch.get_label() == "as-is interval" for patch in figure.axes[0].patches)
    assert any(patch.get_label() == "as-is interval" for patch in figure.axes[1].patches)
    assert plt.get_fignums() == []
    assert result_path.read_bytes() == before


def test_plot_periodic_result_reports_missing_periodic_data(tmp_path: Path) -> None:
    result_path = tmp_path / "missing.h5"
    with h5py.File(result_path, "w") as handle:
        handle.create_dataset("backs/e/r", data=np.array([8.0, 9.0]))

    with pytest.raises(ResultError, match="required periodic dataset.*fields/Phi"):
        plot_periodic_result(result_path, tmp_path / "missing.png")
    assert not (tmp_path / "missing.png").exists()
    assert plt.get_fignums() == []


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_plot_periodic_result_rejects_alias_of_input_without_modifying_it(
    tmp_path: Path, alias_kind: str
) -> None:
    result_path = _periodic_file(tmp_path / "result.h5")
    figure_path = tmp_path / "figure.png"
    if alias_kind == "direct":
        figure_path = result_path
    elif alias_kind == "symlink":
        figure_path.symlink_to(result_path)
    else:
        os.link(result_path, figure_path)
    before = result_path.read_bytes()

    with pytest.raises(ValueError, match="different files"):
        plot_periodic_result(result_path, figure_path)

    assert result_path.read_bytes() == before
    assert plt.get_fignums() == []


def test_plot_periodic_result_closes_figure_when_output_parent_cannot_be_created(
    tmp_path: Path,
) -> None:
    result_path = _periodic_file(tmp_path / "result.h5")
    blocked_parent = tmp_path / "existing-file"
    blocked_parent.write_text("not a directory", encoding="utf-8")

    with pytest.raises(OSError):
        plot_periodic_result(result_path, blocked_parent / "figure.png")

    assert plt.get_fignums() == []
