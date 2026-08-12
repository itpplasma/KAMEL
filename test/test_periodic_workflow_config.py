"""Regression tests for the explicit periodic KIM production contract."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from balance_interface.balance_conf import balance_conf


def test_periodic_kim_configuration_is_explicit():
    conf = balance_conf()
    conf.configure_periodic_kim([(6, 2), (-6, 2)], target_current=4.0)
    nml = conf.conf["balancenml"]
    assert nml["wave_code"] == "KIM"
    assert nml["kim_run_type"] == "electrostatic_periodic"
    assert nml["kim_electron_transport_model"] == "drift_kinetic"
    assert nml["kim_ion_transport_model"] == "integral"
    assert nml["kim_bparallel_source"] == "periodic"
    assert nml["kim_m_list"] == [6, -6]
    assert nml["I_par_toroidal"] == 4.0


@pytest.mark.parametrize("modes", [[], [(0, 2)], [(1, 0)]])
def test_periodic_kim_rejects_invalid_modes(modes):
    conf = balance_conf()
    with pytest.raises(ValueError):
        conf.configure_periodic_kim(modes)
