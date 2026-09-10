"""Regression tests for the explicit periodic KIM workflow contract."""

import hashlib
import math
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
    assert nml["kim_ion_transport_model"] == "finite_larmor_radius"
    assert nml["kim_bparallel_source"] == "disabled"
    assert nml["kim_transport_benchmark"] is False
    assert nml["kim_m_list"] == [6, -6]
    assert nml["I_par_toroidal"] == 4.0


@pytest.mark.parametrize(
    "modes",
    [
        [],
        [(0, 2)],
        [(1, 0)],
        [(1.5, 2)],
        [(True, 2)],
        [(1, 2), (1, 2)],
        [1],
        [(1, 2, 3)],
        ["12"],
    ],
)
def test_periodic_kim_rejects_invalid_modes(modes):
    conf = balance_conf()
    with pytest.raises(ValueError):
        conf.configure_periodic_kim(modes)


@pytest.mark.parametrize("target", [-1.0, math.inf, -math.inf, math.nan, True])
def test_periodic_kim_rejects_invalid_target_current(target):
    conf = balance_conf()
    with pytest.raises(ValueError):
        conf.configure_periodic_kim([(1, 2)], target_current=target)


def test_periodic_kim_benchmark_mode_controls_solver_flag():
    conf = balance_conf()
    conf.configure_periodic_kim([(6, 2)], benchmark_mode="drift_kinetic_limit")
    assert conf.conf["balancenml"]["kim_transport_benchmark"] is True


def test_run_periodic_kim_prepares_all_modes(tmp_path, monkeypatch):
    interface = QL_Balance_interface(tmp_path, 1, 2.0, "periodic", debug=False)
    interface.conf = balance_conf()
    observed = {}

    def fake_prepare(Btor, a_minor):
        observed["modes"] = (interface.m_mode, interface.n_mode)

    monkeypatch.setattr(interface, "prepare_balance_kim", fake_prepare)
    monkeypatch.setattr(interface, "set_config_nml", lambda: None)
    monkeypatch.setattr(interface, "write_config_nml", lambda path: None)
    monkeypatch.setattr(interface, "run_balance", lambda **kwargs: 0)

    kim_config = tmp_path / "source_KIM_config.nml"
    kim_config.write_text(
        "&KIM_CONFIG\n turn_off_ions=.false., turn_off_electrons=.false.\n/\n"
        "&KIM_PERIODIC\n periodic_Bparallel_ratio=(0.0,0.0)\n/\n"
    )
    result = interface.run_periodic_kim(-2.0, 50.0, [(-6, 2), (-7, 2)], kim_config_file=kim_config)
    assert result == 0
    assert observed["modes"] == ([-6, -7], [2, 2])
    expected_digest = hashlib.sha256(kim_config.read_bytes()).hexdigest()
    assert interface.conf.conf["balancenml"]["kim_config_sha256"] == expected_digest


def test_run_periodic_kim_requires_kim_config(tmp_path):
    interface = QL_Balance_interface(tmp_path, 1, 2.0, "periodic", debug=False)
    interface.conf = balance_conf()
    with pytest.raises(ValueError, match="kim_config_file is required"):
        interface.run_periodic_kim(-2.0, 50.0, [(-6, 2)])


@pytest.mark.parametrize(
    "kim_config",
    [
        "&KIM_CONFIG\n turn_off_ions=.true., turn_off_electrons=.false.\n/\n",
        "&KIM_CONFIG\n turn_off_ions=.false., turn_off_electrons=.true.\n/\n",
    ],
)
def test_run_periodic_kim_rejects_disabled_species(tmp_path, kim_config):
    config_path = tmp_path / "source_KIM_config.nml"
    config_path.write_text(kim_config)
    interface = QL_Balance_interface(tmp_path, 1, 2.0, "periodic", debug=False)
    interface.conf = balance_conf()
    with pytest.raises(ValueError):
        interface.run_periodic_kim(-2.0, 50.0, [(-6, 2)], kim_config_file=config_path)


def test_run_periodic_kim_accepts_prescribed_bparallel(tmp_path, monkeypatch):
    config_path = tmp_path / "source_KIM_config.nml"
    config_path.write_text(
        "&KIM_CONFIG\n"
        " turn_off_ions=.false., turn_off_electrons=.false.,\n"
        " collision_model='FokkerPlanck', ion_collision_model='FokkerPlanck',\n"
        " artificial_debye_case=0\n/\n"
        "&KIM_SETUP\n mphi_max=0\n/\n"
        "&KIM_PERIODIC\n"
        " periodic_Bparallel_ratio=(0.25,-0.1),\n"
        " periodic_match_global_kernel_approximations=.false.\n/\n"
    )
    interface = QL_Balance_interface(tmp_path, 1, 2.0, "periodic", debug=False)
    interface.conf = balance_conf()
    monkeypatch.setattr(interface, "prepare_balance_kim", lambda *_: None)
    monkeypatch.setattr(interface, "set_config_nml", lambda: None)
    monkeypatch.setattr(interface, "write_config_nml", lambda path: None)
    monkeypatch.setattr(interface, "run_balance", lambda **kwargs: 0)
    interface.run_periodic_kim(-2.0, 50.0, [(-6, 2)], kim_config_file=config_path)
    assert interface.conf.conf["balancenml"]["kim_bparallel_source"] == "prescribed_zero_mode"


@pytest.mark.parametrize(
    ("physics_settings", "setup_settings", "periodic_settings", "message"),
    [
        ("collision_model='Krook'", "mphi_max=0", "", "FokkerPlanck collisions"),
        (
            "collision_model='FokkerPlanck', ion_collision_model='collisionless'",
            "mphi_max=0",
            "",
            "FokkerPlanck ions",
        ),
        (
            "collision_model='FokkerPlanck', artificial_debye_case=1",
            "mphi_max=0",
            "",
            "artificial_debye_case",
        ),
        (
            "collision_model='FokkerPlanck'",
            "mphi_max=1",
            "",
            "mphi_max=0",
        ),
        (
            "collision_model='FokkerPlanck'",
            "mphi_max=0",
            "periodic_match_global_kernel_approximations=.true.",
            "full periodic gyrogeometry",
        ),
    ],
)
def test_run_periodic_kim_rejects_unsupported_prescribed_bparallel(
    tmp_path, physics_settings, setup_settings, periodic_settings, message
):
    config_path = tmp_path / "source_KIM_config.nml"
    config_path.write_text(
        "&KIM_CONFIG\n"
        " turn_off_ions=.false., turn_off_electrons=.false.,\n"
        f" {physics_settings}\n/\n"
        f"&KIM_SETUP\n {setup_settings}\n/\n"
        "&KIM_PERIODIC\n periodic_Bparallel_ratio=(0.25,-0.1),\n"
        f" {periodic_settings}\n/\n"
    )
    interface = QL_Balance_interface(tmp_path, 1, 2.0, "periodic", debug=False)
    interface.conf = balance_conf()
    with pytest.raises(ValueError, match=message):
        interface.run_periodic_kim(-2.0, 50.0, [(-6, 2)], kim_config_file=config_path)
