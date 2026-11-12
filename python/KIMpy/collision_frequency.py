"""Collision frequencies of electrons and ions, mirroring KIM's species_mod.f90.

Densities are in cm^-3 and temperatures in eV. ``spec_dat`` maps a species name
to a dict with ``n`` and ``T`` (arrays or scalars) and, for ions, ``Zi`` and
``Ai``; the electrons are ``"e"``.
"""

import numpy as np

NU_LIKE_ELECTRON = 5.8e-6
NU_ELECTRON_ION = 7.7e-6
NU_ION_ION = 1.8e-7
NU_ION_ELECTRON = 1.8e-7


def electron_coulomb_logs(n_e, T_e):
    """Return the electron-electron and electron-ion Coulomb logarithms (NRL)."""
    log_ee = (
        23.5
        - np.log(np.sqrt(n_e) * T_e ** (-5 / 4))
        - np.sqrt(1e-5 + ((np.log(T_e) - 2) ** 2) / 16)
    )
    log_ei = 24 - np.log(np.sqrt(n_e) / T_e)
    return log_ee, log_ei


def ion_ion_coulomb_log(test, field):
    """Return the NRL mixed ion-ion Coulomb logarithm of two ion species."""
    return 23 - np.log(
        test["Zi"]
        * field["Zi"]
        * (test["Ai"] + field["Ai"])
        / (test["T"] * field["Ai"] + field["T"] * test["Ai"])
        * np.sqrt(
            test["n"] * test["Zi"] ** 2 / test["T"] + field["n"] * field["Zi"] ** 2 / field["T"]
        )
    )


def collision_frequencies(spec_dat, ion_species):
    """Return a dict of collision frequencies [1/s] for the electrons and each ion species."""
    electrons = spec_dat["e"]
    log_ee, log_ei = electron_coulomb_logs(electrons["n"], electrons["T"])

    nu = {"e": NU_LIKE_ELECTRON * electrons["n"] * log_ee * electrons["T"] ** (-3 / 2)}
    for test in ion_species:
        ion = spec_dat[test]
        nu["e"] = (
            nu["e"] + NU_ELECTRON_ION * electrons["T"] ** -1.5 * ion["n"] * log_ei * ion["Zi"] ** 2
        )
        nu[test] = (
            NU_ION_ELECTRON
            * ion["Ai"] ** (-1 / 2)
            * ion["T"] ** (-3 / 2)
            * electrons["n"]
            * ion["Zi"] ** 2
            * log_ei
        )
        for field in ion_species:
            partner = spec_dat[field]
            nu[test] = nu[test] + (
                NU_ION_ION
                * ion["Ai"] ** -0.5
                * ion["T"] ** -1.5
                * partner["n"]
                * ion["Zi"] ** 2
                * partner["Zi"] ** 2
                * ion_ion_coulomb_log(ion, partner)
            )
    return nu
