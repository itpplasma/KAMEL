"""KiLCA antenna current I0: machine defaults and explicit overrides."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from KiLCA_interface.KiLCA_interface import KiLCA_interface

AUG_DEFAULT_I0 = 5.1e12
MASTU_DEFAULT_I0 = 6.0e12


def make_interface(tmp_path, machine):
    return KiLCA_interface(1, 1, str(tmp_path) + "/", "vacuum", machine)


@pytest.mark.filterwarnings("ignore:No profile directory")
@pytest.mark.parametrize(
    ("machine", "default_I0"), [("AUG", AUG_DEFAULT_I0), ("MASTU", MASTU_DEFAULT_I0)]
)
def test_I0_none_keeps_machine_default(tmp_path, machine, default_I0):
    interface = make_interface(tmp_path, machine)
    interface.set_machine(delta_r_antenna=3.0, a_minor=None, I0=None)
    assert interface.antenna.data["I0"] == default_I0


@pytest.mark.filterwarnings("ignore:No profile directory")
@pytest.mark.parametrize("machine", ["AUG", "MASTU"])
def test_explicit_I0_overrides_machine_default(tmp_path, machine):
    interface = make_interface(tmp_path, machine)
    interface.set_machine(delta_r_antenna=3.0, a_minor=None, I0=1.2e12)
    assert interface.antenna.data["I0"] == 1.2e12


@pytest.mark.filterwarnings("ignore:No profile directory")
def test_I0_without_machine_is_rejected(tmp_path):
    interface = make_interface(tmp_path, None)
    with pytest.raises(ValueError, match="I0 needs a machine"):
        interface.set_machine(I0=1.2e12)
