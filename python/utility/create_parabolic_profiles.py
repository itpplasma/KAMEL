import os

import numpy as np


def create_parabolic_profiles_from_res_surf(
    path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, m_mode, n_mode, rmin, rmax, num, a
):
    """Create parabolic profiles for fixed density and electron
    temperature values at the rational surface."""

    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

    r = np.linspace(rmin, rmax, num)
    q = -(1.05 + q0 * (r / a) ** 2)
    rres = np.interp(m_mode / n_mode, np.abs(q), r)

    fac_par = 1 - (r / a) ** 2
    n = n0 * fac_par / (1 - (rres / a) ** 2)
    Te = Te0 * fac_par / (1 - (rres / a) ** 2)
    Ti = Ti0 * fac_par
    Vz = Vz0 * fac_par
    Er = Er0 * fac_par
    Vth = Vth0 * fac_par

    # set profiles zero outside of plasma
    n[np.where(n < n0 / 10)] = n0 / 10
    Te[np.where(Te < Te0 / 10)] = Te0 / 10
    Ti[np.where(Ti < Ti0 / 10)] = Ti0 / 10
    Vz[np.where(Vz < Vz0 / 10)] = Vz0 / 10
    Er[np.where(Er < Er0 / 10)] = Er0 / 10
    Vth[np.where(Vth < Vth0 / 10)] = Vth0 / 10

    np.savetxt(os.path.join(path, "q.dat"), np.array((r, q)).transpose())
    np.savetxt(os.path.join(path, "Te.dat"), np.array((r, Te)).transpose())
    np.savetxt(os.path.join(path, "Ti.dat"), np.array((r, Ti)).transpose())
    np.savetxt(os.path.join(path, "n.dat"), np.array((r, n)).transpose())
    np.savetxt(os.path.join(path, "Vz.dat"), np.array((r, Vz)).transpose())
    np.savetxt(os.path.join(path, "Er.dat"), np.array((r, Er)).transpose())
    np.savetxt(os.path.join(path, "Vth.dat"), np.array((r, Vth)).transpose())
