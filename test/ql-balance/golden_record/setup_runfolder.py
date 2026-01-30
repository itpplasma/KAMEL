#!/usr/bin/env python3
"""
Setup script for the QL-Balance golden record test runfolder.

==============================================================================
USER TODO: Implement this file!

This script must create a complete runfolder with all necessary input files
for running ql-balance.x. The runfolder structure should match the example at:
    ../runfolders/runfolder_33353_2900/

Required files typically include:
    - balance_conf.nml (QL-Balance configuration namelist)
    - background.in (background plasma parameters)
    - eigmode.in (eigenmode configuration)
    - output.in (output configuration)
    - profiles/ directory with:
        - n.dat (electron density)
        - Te.dat (electron temperature)
        - Ti.dat (ion temperature)
        - q.dat (safety factor)
        - Er.dat (radial electric field)
        - Vz.dat (toroidal rotation) - optional
        - Vth.dat (thermal velocity) - optional
        - Da.dat (anomalous diffusion) - optional
    - Equilibrium data file (e.g., *_equil_r_q_psi.dat)
    - Input HDF5 file (e.g., *_inp.hdf5)

You can either:
1. Copy files from a reference runfolder
2. Generate files programmatically
3. Use a combination of both

==============================================================================
"""

import os
import sys
from pathlib import Path

# Add the python package directory to sys.path for standalone execution
_python_dir = Path(__file__).resolve().parents[3] / "python"
if str(_python_dir) not in sys.path:
    sys.path.insert(0, str(_python_dir))

from balance_interface import QL_Balance_interface
from utility import create_parabolic_profiles_from_res_surf


def setup_runfolder(run_path: str) -> None:
    """
    Set up a runfolder in the target directory with all necessary input files.

    This function is called by ensure_golden.py and test_golden_record.py to
    create the runfolder before running ql-balance.x.

    Args:
        run_path: The directory where the runfolder should be created.
                  This directory already exists when this function is called.

    Raises:
        RuntimeError: If setup fails for any reason.
    """
    mpol = 6  # poloidal mode number
    ntor = 2  # toroidal mode number
    profile_path = os.path.join(run_path, "profiles")

    # all values at the resonant surface:
    q0 = mpol / ntor  # 1, safety factor at resonance
    n0 = 2e13  # cm^-3, particle density
    Te0 = 1000  # eV, electron temperature
    Ti0 = 1000  # eV, ion temperature
    Vz0 = 1e6  # cm/s, physical toroidal rotation velocity
    Er0 = 0.2  # statV/cm, radial electric field
    Vth0 = 1e5  # cm/s, physical poloidal rotation velocity
    rmin = 3.0  # cm, minimum effective radius of the generated grid
    rmax = 80.0  # cm, maximum effective radius
    num = 300  # 1, number of radial grid points
    a = 67.0  # cm, plasma radius (outside profiles have small constant values)

    Btor = -17000  # toroidal magnetic field on axis in Gauss

    bi = QL_Balance_interface(run_path=run_path, shot=0, time=0, name="test")
    bi.read_config_nml()
    bi.set_type_of_run(run_type="TimeEvolution")
    bi.set_modes(m_mode=mpol, n_mode=ntor)

    create_parabolic_profiles_from_res_surf(
        profile_path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, mpol, ntor, rmin, rmax, num, a
    )

    bi.prepare_balance(Btor=Btor, a_minor=a)
    bi.set_config_nml()
    bi.conf.conf["balancenml"]["ramp_up_mode"] = 3  # instant max RMP coil current
    bi.conf.conf["balancenml"]["t_max_ramp_up"] = 1.2
    bi.conf.conf["balancenml"]["diagnostics_output"] = True  # write diagnostic data
    bi.conf.conf["balancenml"]["suppression_mode"] = False  # write full LinearProfiles/KinProfiles
    bi.conf.conf["balancenml"]["save_prof_time_step"] = 1  # save profiles every time step
    bi.write_config_nml(path=os.path.join(run_path, "balance_conf.nml"))


def get_output_hdf5_path(runfolder: Path) -> Path:
    """
    Return the path to the output HDF5 file produced by ql-balance.x.

    This function is called after ql-balance.x has finished running to locate
    the output file for comparison.

    Args:
        runfolder: The runfolder directory where ql-balance.x was run.

    Returns:
        Path to the output HDF5 file.

    Raises:
        FileNotFoundError: If the output file doesn't exist.
    """
    # Output file matches the pattern from setup_ql_balance():
    # shot=0, time=0, name="test" -> out/0_0_test.hdf5
    output_file = Path(runfolder) / "out" / "0_0_test.hdf5"
    if not output_file.exists():
        raise FileNotFoundError(f"Output file not found: {output_file}")
    return output_file
