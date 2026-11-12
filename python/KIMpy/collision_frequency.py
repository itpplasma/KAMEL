"""Collision frequencies of electrons and ions, mirroring KIM's species_mod.f90.

Densities are in cm^-3 and temperatures in eV. ``spec_dat`` maps a species name
to a dict with ``n`` and ``T`` (arrays or scalars) and, for ions, ``Zi`` and
``Ai``; the electrons are ``"e"``.
"""

import numpy as np

# KiLCA's convention (KiLCA/background/background_data_m.f90). The ion
# coefficients are the electron ones scaled by the mass ratio: like particles
# 5.8e-6 sqrt(m_e/m_p) = 1.4e-7, slow heavy field particles 7.7e-6 sqrt(m_e/m_p)
# = 1.8e-7, ions on electrons 5.8e-6 m_e/m_p = 3.2e-9. Unlike ion pairs keep
# the heavy-field coefficient, which is not validated for lighter field ions.
NU_LIKE_ELECTRON = 5.8e-6
NU_ELECTRON_ION = 7.7e-6
NU_LIKE_ION = 1.4e-7
NU_HEAVY_FIELD_ION = 1.8e-7
NU_ION_ELECTRON = 3.2e-9


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
        # Ions on electrons: momentum exchange suppressed by m_e/m_i,
        # with Lambda_ie = Lambda_ei.
        nu[test] = (
            NU_ION_ELECTRON
            * electrons["n"]
            * log_ei
            * ion["Zi"] ** 2
            / ion["Ai"]
            / (np.sqrt(electrons["T"]) * ion["T"])
        )
        for field in ion_species:
            partner = spec_dat[field]
            coefficient = NU_LIKE_ION if field == test else NU_HEAVY_FIELD_ION
            nu[test] = nu[test] + (
                coefficient
                * partner["n"]
                * ion["Zi"] ** 2
                * partner["Zi"] ** 2
                * ion_ion_coulomb_log(ion, partner)
                / np.sqrt(ion["Ai"])
                / ion["T"] ** 1.5
            )
    return nu
