"""
KIMPoissonSolver -- load KIM HDF5 kernels and solve the Poisson equation.

Provides a class that loads kernel matrices, background quantities, and the
aligned potential from a KIM electrostatic HDF5 output file, assembles the
Poisson equation, and solves for the potential phi.

Example
-------
>>> solver = KIMPoissonSolver("path/to/out_ES_D_bench.h5", m_mode=6, n_mode=2)
>>> solver.set_delta_br(offset=-5.0)
>>> solver.solve()
>>> r_fine, phi_fine = solver.spline_phi()
>>> kr = solver.compute_local_kr()
"""

import h5py
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares


def impose_poisson_boundary_conditions(a_mat, b_vec, phi_aligned):
    """Impose zero-misalignment BC (Phi = -phi_aligned at boundaries).

    Parameters
    ----------
    a_mat : ndarray, shape (n, n)
        System matrix (Laplacian + density-response kernels).
    b_vec : ndarray, shape (n,)
        Right-hand side vector.
    phi_aligned : ndarray, shape (n,)
        Aligned potential on the field-line grid.

    Returns
    -------
    a_mat, b_vec : tuple of ndarray
        Modified system with Dirichlet BCs applied.
    """
    a_mat = a_mat.copy()
    b_vec = b_vec.copy()
    npts = len(b_vec)

    phi_left = -phi_aligned[0]
    phi_right = -phi_aligned[-1]

    b_vec[1 : npts - 1] = (
        b_vec[1 : npts - 1]
        - a_mat[1 : npts - 1, 0] * phi_left
        - a_mat[1 : npts - 1, npts - 1] * phi_right
    )

    a_mat[:, 0] = 0.0
    a_mat[:, -1] = 0.0
    a_mat[0, :] = 0.0
    a_mat[-1, :] = 0.0
    a_mat[0, 0] = 1.0
    a_mat[-1, -1] = 1.0
    b_vec[0] = phi_left
    b_vec[-1] = phi_right
    return a_mat, b_vec


class KIMPoissonSolver:
    """Load KIM kernels from HDF5 and solve the Poisson equation.

    Parameters
    ----------
    h5_path : str
        Path to KIM electrostatic HDF5 output file
        (e.g. ``out_ES_D_bench.h5``).
    m_mode : int
        Poloidal mode number (default 6).
    n_mode : int
        Toroidal mode number (default 2).

    Attributes
    ----------
    xl : ndarray, shape (n,)
        Field-line grid.
    laplace : ndarray, shape (n, n)
        Laplacian FEM matrix.
    k_rho_phi : ndarray, shape (n, n)
        Combined (electron + ion) density-response kernel for phi.
    k_rho_b : ndarray, shape (n, n)
        Combined (electron + ion) density-response kernel for B.
    phi_aligned : ndarray, shape (n,)
        Aligned potential ``1j * E0r / (kp * B0)`` interpolated onto *xl*.
    r_res : float
        Resonant-surface position [cm].
    rho_l : ndarray
        Ion Larmor radius profile on the background grid.
    rb : ndarray
        Background radial grid.
    """

    def __init__(self, h5_path, m_mode=6, n_mode=2):
        self.h5_path = h5_path
        self.m_mode = m_mode
        self.n_mode = n_mode

        self._load()
        self.r_res = self._compute_r_res()

    # ------------------------------------------------------------------
    def _load(self):
        """Read HDF5 kernels, backgrounds, and grid data."""
        with h5py.File(self.h5_path, "r") as h5f:
            xl = h5f["grid/xl_xb"][:]
            rg = h5f["grid/rg_xb"][:]
            e0r = h5f["backs/E0r"][:]
            kp = h5f["backs/kp"][:]
            b0 = h5f["backs/B0"][:]
            rho_l = h5f["backs/i/rho_L"][:]
            rb = h5f["backs/i/r"][:]
            q = h5f["backs/q"][:]

            laplace = h5f["kernel/Laplace_in_FEM"][:]
            k_rho_phi_e = h5f["kernel/K_rho_phi_e_re"][:] + 1j * h5f["kernel/K_rho_phi_e_im"][:]
            k_rho_phi_i = h5f["kernel/K_rho_phi_i_re"][:] + 1j * h5f["kernel/K_rho_phi_i_im"][:]
            k_rho_b_e = h5f["kernel/K_rho_B_e_re"][:] + 1j * h5f["kernel/K_rho_B_e_im"][:]
            k_rho_b_i = h5f["kernel/K_rho_B_i_re"][:] + 1j * h5f["kernel/K_rho_B_i_im"][:]

        # B0 may live on a different grid; use mean as reference
        b0_ref = np.mean(b0)

        self.xl = xl
        self.laplace = laplace
        self.k_rho_phi = k_rho_phi_e + k_rho_phi_i
        self.k_rho_b = k_rho_b_e + k_rho_b_i
        self.phi_aligned = 1j * np.interp(xl, rg, e0r) / (np.interp(xl, rg, kp) * b0_ref)
        self.rho_l = rho_l
        self.rb = rb
        self._q = q
        self.br = None
        self.r_probe = None
        self.phi = None

    # ------------------------------------------------------------------
    def _compute_r_res(self):
        """Find the resonant surface where |q(r)| = m/n.

        Uses linear interpolation between grid points to find the precise
        crossing.  Returns the outermost crossing (relevant for edge modes).

        Returns
        -------
        float
            Radial position of the resonant surface [cm].
        """
        q_res = self.m_mode / self.n_mode
        q_abs = np.abs(self._q)
        r = self.rb

        crossings = []
        for i in range(len(q_abs) - 1):
            if (q_abs[i] - q_res) * (q_abs[i + 1] - q_res) < 0:
                r_cross = r[i] + (q_res - q_abs[i]) * (r[i + 1] - r[i]) / (q_abs[i + 1] - q_abs[i])
                crossings.append(r_cross)

        if not crossings:
            raise ValueError(f"No resonant surface found for m={self.m_mode}, " f"n={self.n_mode}")
        return crossings[-1]

    # ------------------------------------------------------------------
    def set_delta_br(self, offset=2.0):
        """Create a delta-function perturbation in B_r at *r_res + offset*.

        Parameters
        ----------
        offset : float
            Radial distance from the resonant surface [cm] (default 2.0).
        """
        self.br = np.zeros_like(self.xl, dtype=np.complex128)
        r_p = self.r_res + offset
        idx = int(np.argmin(np.abs(self.xl - r_p)))
        self.br[idx] = 1.0
        self.r_probe = self.xl[idx]

    # ------------------------------------------------------------------
    def solve(self):
        """Assemble and solve the Poisson equation for the potential *phi*.

        Raises
        ------
        RuntimeError
            If :meth:`set_delta_br` has not been called first.
        """
        if self.br is None:
            raise RuntimeError("self.br is None -- call set_delta_br() before solve()")

        a_mat = self.laplace + 4.0 * np.pi * self.k_rho_phi
        b_vec = -4.0 * np.pi * (self.k_rho_b @ self.br)
        a_mat, b_vec = impose_poisson_boundary_conditions(a_mat, b_vec, self.phi_aligned)
        self.phi = np.linalg.solve(a_mat, b_vec)

    # ------------------------------------------------------------------
    def spline_phi(self, n_points=5000):
        """Spline phi onto a fine equidistant grid.

        The original grid *xl* is non-equidistant (dense near resonance,
        sparse away from it).  Splining gives better resolution for
        plotting and kr extraction away from the resonant surface.

        Parameters
        ----------
        n_points : int
            Number of points in the fine grid (default 5000).

        Returns
        -------
        r_fine : ndarray
            Equidistant radial grid.
        phi_fine : ndarray, complex
            Splined potential on the fine grid.
        """
        if self.phi is None:
            raise RuntimeError("self.phi is None -- call solve() before spline_phi()")

        r_fine = np.linspace(self.xl[0], self.xl[-1], n_points)
        cs_re = CubicSpline(self.xl, self.phi.real)
        cs_im = CubicSpline(self.xl, self.phi.imag)
        phi_fine = cs_re(r_fine) + 1j * cs_im(r_fine)
        return r_fine, phi_fine

    # ------------------------------------------------------------------
    def compute_local_kr(self):
        """Compute the local radial wavenumber from the solved potential.

        Returns
        -------
        kr : ndarray, complex
            Local wavenumber ``kr(r) = -i d(ln phi)/dr``.

        Raises
        ------
        RuntimeError
            If :meth:`solve` has not been called first.
        """
        if self.phi is None:
            raise RuntimeError("self.phi is None -- call solve() before compute_local_kr()")
        return -1j * np.gradient(np.log(self.phi), self.xl)

    # ------------------------------------------------------------------
    def fit_kr(self, r_min, r_max):
        """Fit a complex exponential ``A exp(i kr (r - r_probe))`` to *phi*.

        Parameters
        ----------
        r_min, r_max : float
            Radial window [cm] for the fit.

        Returns
        -------
        dict
            ``{'kr': complex, 'amplitude': complex, 'r_fit': ndarray,
            'phi_fit': ndarray, 'success': bool}``

        Raises
        ------
        RuntimeError
            If :meth:`solve` (or manual assignment) has not set *phi*.
        """
        if self.phi is None:
            raise RuntimeError("self.phi is None -- call solve() before fit_kr()")

        mask = (self.xl >= r_min) & (self.xl <= r_max)
        r_fit = self.xl[mask]
        phi_fit = self.phi[mask]
        r_shifted = r_fit - self.r_probe

        # Initial guess from local wavenumber at window centre
        kr_local = self.compute_local_kr()
        kr_local_fit = kr_local[mask]
        mid_idx = len(kr_local_fit) // 2
        kr0 = kr_local_fit[mid_idx]
        a0 = phi_fit[mid_idx]

        p0 = [kr0.real, kr0.imag, a0.real, a0.imag]

        def _residuals(p):
            kr_c = p[0] + 1j * p[1]
            a_c = p[2] + 1j * p[3]
            model = a_c * np.exp(1j * kr_c * r_shifted)
            diff = model - phi_fit
            return np.concatenate([diff.real, diff.imag])

        result = least_squares(_residuals, p0, method="lm")

        kr_fit = result.x[0] + 1j * result.x[1]
        amp_fit = result.x[2] + 1j * result.x[3]
        return {
            "kr": kr_fit,
            "amplitude": amp_fit,
            "r_fit": r_fit,
            "phi_fit": phi_fit,
            "success": bool(result.success),
        }

    # ------------------------------------------------------------------
    def fit_exponential_components(self, r_left, r_right=None, xl=None, phi=None):
        """Fit Re(phi) and Im(phi) with ``A exp(b (r - r_probe)) + d``.

        The fit is performed on a window to the left of the probe by default
        (i.e. *r < r_probe*), capturing the outward-propagating wave.

        Parameters
        ----------
        r_left : float
            Left edge of the fitting window [cm].
        r_right : float or None
            Right edge of the fitting window [cm].
            Defaults to ``r_probe`` if *None*.
        xl : ndarray or None
            Radial grid to use. Defaults to ``self.xl``.
        phi : ndarray or None
            Potential array to use. Defaults to ``self.phi``.

        Returns
        -------
        dict with keys:
            ``'r_fit'`` : radial grid in window,
            ``'fit_re'`` : dict ``{'A', 'b', 'd', 'model', 'success'}``
            for Re(phi),
            ``'fit_im'`` : dict ``{'A', 'b', 'd', 'model', 'success'}``
            for Im(phi).
        """
        if xl is None:
            xl = self.xl
        if phi is None:
            phi = self.phi
        if phi is None:
            raise RuntimeError("phi is None -- call solve() before " "fit_exponential_components()")

        if r_right is None:
            r_right = self.r_probe

        mask = (xl >= r_left) & (xl <= r_right)
        r_fit = xl[mask]
        phi_fit = phi[mask]
        r_shifted = r_fit - self.r_probe

        def _fit_component(y_data):
            """Fit y = A * exp(b * x) + d to real data."""
            out = {"A": None, "b": None, "d": None, "model": None, "success": False}
            if len(y_data) < 4:
                return out

            valid = np.abs(y_data) > 0
            if np.sum(valid) < 4:
                out["A"] = 0.0
                out["b"] = 0.0
                out["d"] = 0.0
                out["model"] = np.zeros_like(y_data)
                out["success"] = True
                return out

            log_abs = np.log(np.abs(y_data[valid]))
            b0 = np.polyfit(r_shifted[valid], log_abs, 1)[0]
            a0 = y_data[-1]

            def residuals(p):
                return p[0] * np.exp(p[1] * r_shifted) + p[2] - y_data

            try:
                res = least_squares(residuals, [a0, b0, 0.0], method="lm")
                out["A"] = res.x[0]
                out["b"] = res.x[1]
                out["d"] = res.x[2]
                out["model"] = res.x[0] * np.exp(res.x[1] * r_shifted) + res.x[2]
                out["success"] = bool(res.success)
            except Exception:
                pass
            return out

        return {
            "r_fit": r_fit,
            "fit_re": _fit_component(phi_fit.real),
            "fit_im": _fit_component(phi_fit.imag),
        }

    # ------------------------------------------------------------------
    def fit_sinexp_components(self, r_left, r_right=None, xl=None, phi=None):
        """Fit Re/Im(phi) with ``A sin(c (r-r_p) + d) exp(b (r-r_p))``.

        Parameters
        ----------
        r_left : float
            Left edge of the fitting window [cm].
        r_right : float or None
            Right edge. Defaults to ``r_probe``.
        xl : ndarray or None
            Radial grid. Defaults to ``self.xl``.
        phi : ndarray or None
            Potential array. Defaults to ``self.phi``.

        Returns
        -------
        dict with keys:
            ``'r_fit'`` : radial grid in window,
            ``'fit_re'`` : dict ``{'A', 'b', 'c', 'd', 'model', 'success'}``,
            ``'fit_im'`` : dict ``{'A', 'b', 'c', 'd', 'model', 'success'}``.
        """
        if xl is None:
            xl = self.xl
        if phi is None:
            phi = self.phi
        if phi is None:
            raise RuntimeError("phi is None -- call solve() before " "fit_sinexp_components()")

        if r_right is None:
            r_right = self.r_probe

        mask = (xl >= r_left) & (xl <= r_right)
        r_fit = xl[mask]
        phi_fit = phi[mask]
        r_shifted = r_fit - self.r_probe

        def _fit_component(y_data):
            """Fit y = A * sin(c * x + d) * exp(b * x) to real data."""
            out = {"A": None, "b": None, "c": None, "d": None, "model": None, "success": False}
            if len(y_data) < 6:
                return out

            valid = np.abs(y_data) > 0
            if np.sum(valid) < 6:
                return out

            log_abs = np.log(np.abs(y_data[valid]))
            b0 = np.polyfit(r_shifted[valid], log_abs, 1)[0]

            sign_changes = np.diff(np.sign(y_data))
            n_zeros = np.count_nonzero(sign_changes)
            span = r_shifted[-1] - r_shifted[0]
            if n_zeros > 0 and abs(span) > 0:
                c0 = n_zeros * np.pi / abs(span)
            else:
                c0 = 2.0

            a0 = np.max(np.abs(y_data))

            def residuals(p):
                return p[0] * np.sin(p[2] * r_shifted + p[3]) * np.exp(p[1] * r_shifted) - y_data

            best = None
            for a_sign in [1.0, -1.0]:
                for c_mult in [1.0, 0.5, 2.0]:
                    for d0 in [0.0, np.pi / 2]:
                        p0 = [a_sign * a0, b0, c0 * c_mult, d0]
                        try:
                            res = least_squares(residuals, p0, method="trf")
                            if res.success and (best is None or res.cost < best.cost):
                                best = res
                        except Exception:
                            pass

            if best is not None and best.success:
                out["A"] = best.x[0]
                out["b"] = best.x[1]
                out["c"] = best.x[2]
                out["d"] = best.x[3]
                out["model"] = (
                    best.x[0]
                    * np.sin(best.x[2] * r_shifted + best.x[3])
                    * np.exp(best.x[1] * r_shifted)
                )
                out["success"] = True
            return out

        return {
            "r_fit": r_fit,
            "fit_re": _fit_component(phi_fit.real),
            "fit_im": _fit_component(phi_fit.imag),
        }
