"""Prepare a QL-Balance tutorial run directory with KAMELpy.

The script creates parabolic plasma profiles, runs KiLCA for the plasma and
the vacuum response, writes ``balance_conf.nml`` and links ``ql-balance.x``
into the run directory. Run ``./ql-balance`` there afterwards.

Usage: python setup_run.py RUN_DIRECTORY
"""

import argparse
from pathlib import Path

from balance_interface import QL_Balance_interface
from KiLCA_interface.KiLCA_zone_flre import KiLCA_zone_flre
from utility import create_parabolic_profiles_from_res_surf

# Resonant surface and plasma values at the resonance (CGS units).
M_MODE = 6  # poloidal mode number
N_MODE = 2  # toroidal mode number
Q_RES = M_MODE / N_MODE  # safety factor at the resonance
DENSITY = 2e13  # 1/cm^3
ELECTRON_TEMPERATURE = 1000.0  # eV
ION_TEMPERATURE = 1000.0  # eV
TOROIDAL_VELOCITY = 1e6  # cm/s
RADIAL_ELECTRIC_FIELD = 0.2  # statV/cm
POLOIDAL_VELOCITY = 1e5  # cm/s

# Geometry and profile grid.
R_MIN = 3.0  # cm, first profile point
R_MAX = 80.0  # cm, last profile point
N_POINTS = 300  # number of profile points
PLASMA_RADIUS = 67.0  # cm
TOROIDAL_FIELD = -17000.0  # G, on axis

# A coarser KiLCA conductivity grid than the default keeps the run short.
KILCA_CONDUCTIVITY_GRID = 3001


def setup_run(run_path: Path) -> None:
    run_path.mkdir(parents=True, exist_ok=True)
    KiLCA_zone_flre.data["cond_maxdim"] = KILCA_CONDUCTIVITY_GRID

    balance = QL_Balance_interface(run_path=run_path, shot=0, time=0, name="tutorial")
    balance.read_config_nml()  # start from QL-Balance/namelists/balance_conf.nml
    balance.set_type_of_run(run_type="TimeEvolution")
    balance.set_modes(m_mode=M_MODE, n_mode=N_MODE)

    create_parabolic_profiles_from_res_surf(
        str(run_path / "profiles"),
        Q_RES,
        DENSITY,
        ELECTRON_TEMPERATURE,
        ION_TEMPERATURE,
        TOROIDAL_VELOCITY,
        RADIAL_ELECTRIC_FIELD,
        POLOIDAL_VELOCITY,
        M_MODE,
        N_MODE,
        R_MIN,
        R_MAX,
        N_POINTS,
        PLASMA_RADIUS,
    )

    # Runs KiLCA for the plasma (flre/) and vacuum (vacuum/) response and
    # links the ql-balance executable into the run directory.
    balance.prepare_balance(Btor=TOROIDAL_FIELD, a_minor=PLASMA_RADIUS)

    balance.set_config_nml()
    settings = balance.conf.conf["balancenml"]
    settings["ramp_up_mode"] = 3  # switch the full coil current on at once
    settings["t_max_ramp_up"] = 1.2
    settings["suppression_mode"] = False  # save the full profile history
    settings["save_prof_time_step"] = 1  # save profiles at every time step
    balance.write_config_nml(path=str(run_path / "balance_conf.nml"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_path", type=Path, help="directory to create the run in")
    setup_run(parser.parse_args().run_path.resolve())
