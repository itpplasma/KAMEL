"""KIMpy collision frequencies against KIM's Fortran rates (species_mod.f90)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from KIMpy.collision_frequency import collision_frequencies

# Deuterium, n_e = n_D = 1e13 cm^-3, T_e = 1 keV, T_i = 800 eV. The ion value is
# the one KIM's Fortran test_collision_scale_z0 asserts for calculate_plasma_backs.
FORTRAN_DEUTERIUM_NU = 785.04935762193702


def deuterium_plasma():
    return {
        "e": {"n": 1.0e13, "T": 1.0e3},
        "D": {"n": 1.0e13, "T": 8.0e2, "Zi": 1, "Ai": 2.0},
    }


def two_ion_plasma(order):
    ions = {
        "D": {"n": 7.0e12, "T": 8.0e2, "Zi": 1, "Ai": 2.0},
        "H": {"n": 3.0e12, "T": 6.0e2, "Zi": 1, "Ai": 1.0},
    }
    plasma = {"e": {"n": 1.0e13, "T": 1.0e3}}
    plasma.update({name: ions[name] for name in order})
    return plasma


def test_deuterium_matches_fortran():
    nu = collision_frequencies(deuterium_plasma(), ["D"])
    assert nu["D"] == pytest.approx(FORTRAN_DEUTERIUM_NU, rel=1e-12)


def test_electrons_follow_kilca_rates():
    log_ee = (
        23.5 - np.log(np.sqrt(1.0e13) / 1.0e3**1.25) - np.sqrt(1e-5 + (np.log(1.0e3) - 2) ** 2 / 16)
    )
    log_ei = 24 - np.log(np.sqrt(1.0e13) / 1.0e3)
    expected = (5.8e-6 * log_ee + 7.7e-6 * log_ei) * 1.0e13 / 1.0e3**1.5
    nu = collision_frequencies(deuterium_plasma(), ["D"])
    assert nu["e"] == pytest.approx(expected, rel=1e-12)


def test_ion_electron_collisions_are_mass_suppressed():
    plasma = deuterium_plasma()
    with_electrons = collision_frequencies(plasma, ["D"])["D"]
    plasma["e"]["n"] = 1.0e-30
    without_electrons = collision_frequencies(plasma, ["D"])["D"]
    assert (with_electrons - without_electrons) / without_electrons < 0.05


@pytest.mark.parametrize("species", ["D", "H"])
def test_ion_rates_independent_of_species_order(species):
    forward = collision_frequencies(two_ion_plasma(["D", "H"]), ["D", "H"])
    backward = collision_frequencies(two_ion_plasma(["H", "D"]), ["H", "D"])
    assert backward[species] == pytest.approx(forward[species], rel=1e-14)
