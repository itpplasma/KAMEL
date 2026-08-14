"""Regression tests for the explicit periodic KIM production contract."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from balance_interface.balance_conf import balance_conf
from balance_interface.balance_interface import QL_Balance_interface


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


def test_run_periodic_kim_derives_preparation_mode(tmp_path, monkeypatch):
    interface = QL_Balance_interface(tmp_path, 1, 2.0, "periodic", debug=False)
    interface.conf = balance_conf()
    observed = {}

    def fake_prepare(Btor, a_minor):
        observed["mode"] = (interface.m_mode, interface.n_mode)

    monkeypatch.setattr(interface, "prepare_balance_kim", fake_prepare)
    monkeypatch.setattr(interface, "set_config_nml", lambda: None)
    monkeypatch.setattr(interface, "write_config_nml", lambda path: None)
    monkeypatch.setattr(interface, "run_balance", lambda **kwargs: 0)

    result = interface.run_periodic_kim(-2.0, 50.0, [(-6, 2), (-7, 2)])
    assert result == 0
    assert observed["mode"] == (-6, 2)
