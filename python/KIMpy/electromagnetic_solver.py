"""
KIMElectromagneticSolver -- solve the coupled Poisson-Ampere block system.

Loads kernel matrices from a KIM electromagnetic HDF5 output file,
constructs the Ampere FEM matrices, and solves the coupled 2N x 2N
block system for Phi and Br self-consistently. Can also solve the
Poisson-only (electrostatic) problem for comparison.

Example
-------
>>> solver = KIMElectromagneticSolver("path/to/em_output.h5")
>>> solver.solve(Br_boundary=1.0)
>>> solver.solve_poisson(br_prescribed=np.ones(solver.N))
"""

import h5py
import numpy as np
from scipy.interpolate import CubicSpline

from .constants import sol, pi


class KIMElectromagneticSolver:
    """Solve coupled Poisson-Ampere or Poisson-only from KIM HDF5 kernels.

    Parameters
    ----------
    h5_path : str
        Path to KIM electromagnetic HDF5 output file.
    m_mode : int or None
        Poloidal mode number. Read from HDF5 ``setup/m_mode`` if None.
    n_mode : int or None
        Toroidal mode number. Read from HDF5 ``setup/n_mode`` if None.
    R0 : float or None
        Major radius [cm]. Read from HDF5 ``setup/R0`` if None.
    """

    def __init__(self, h5_path, m_mode=None, n_mode=None, R0=None):
        self.h5_path = h5_path
        self._load(m_mode, n_mode, R0)
        self.N = len(self.xl)
        self.r_res = self._compute_r_res()
        self._build_fem_matrices()
        # B0 on xl grid from geometry (avoids HDF5 grid mismatch for backs/B0)
        self._b0_xl = np.abs(self.btor) * np.sqrt(1.0 + (self.xl / (self._q_xl * self.R0))**2)
        # Aligned potential for BCs
        self.phi_aligned = 1j * self._e0r_xl / (self._kp_xl * self._b0_xl)

    def _load(self, m_mode, n_mode, R0):
        """Read HDF5 kernels, backgrounds, grid, and setup data."""
        with h5py.File(self.h5_path, "r") as h5f:
            # Setup parameters
            if m_mode is None:
                m_mode = int(h5f["setup/m_mode"][()])
            if n_mode is None:
                n_mode = int(h5f["setup/n_mode"][()])
            if R0 is None:
                R0 = float(h5f["setup/R0"][()])

            self.m_mode = m_mode
            self.n_mode = n_mode
            self.R0 = R0
            self.kz = n_mode / R0
            self.btor = float(h5f["setup/btor"][()])

            # Grid
            self.xl = h5f["grid/xl_xb"][:]
            rg = h5f["grid/rg_xb"][:]

            # Background quantities on rg grid
            e0r = h5f["backs/E0r"][:]
            kp = h5f["backs/kp"][:]
            q = h5f["backs/q"][:]

            # Interpolate backgrounds onto xl grid
            self._e0r_xl = np.interp(self.xl, rg, e0r)
            self._kp_xl = np.interp(self.xl, rg, kp)
            self._q_xl = np.interp(self.xl, rg, q)
            self._q_rg = q
            self._rg = rg

            # Laplace matrix
            self.laplace = h5f["kernel/Laplace_in_FEM"][:]

            # Density-response kernels (total, pre-summed over all species)
            self.k_rho_phi = h5f["kernel/K_rho_phi_re"][:] + 1j * h5f["kernel/K_rho_phi_im"][:]
            self.k_rho_b = h5f["kernel/K_rho_B_re"][:] + 1j * h5f["kernel/K_rho_B_im"][:]

            # Current-response kernels (total, pre-summed over all species)
            self.k_j_phi = h5f["kernel/K_j_phi_re"][:] + 1j * h5f["kernel/K_j_phi_im"][:]
            self.k_j_b = h5f["kernel/K_j_B_re"][:] + 1j * h5f["kernel/K_j_B_im"][:]

        # Solution storage
        self.phi = None
        self.br = None
        self.apar = None
        self.phi_es = None

    def _compute_r_res(self):
        """Find the resonant surface where |q(r)| = m/n."""
        q_res = self.m_mode / self.n_mode
        q_abs = np.abs(self._q_rg)
        r = self._rg
        crossings = []
        for i in range(len(q_abs) - 1):
            if (q_abs[i] - q_res) * (q_abs[i + 1] - q_res) < 0:
                r_cross = r[i] + (q_res - q_abs[i]) * (r[i + 1] - r[i]) / (
                    q_abs[i + 1] - q_abs[i]
                )
                crossings.append(r_cross)
        if not crossings:
            raise ValueError(
                f"No resonant surface found for m={self.m_mode}, n={self.n_mode}"
            )
        return crossings[-1]

    def _build_fem_matrices(self):
        """Construct mass matrix, potential matrix, and weighted variants."""
        N = self.N
        xl = self.xl

        # --- Mass matrix M (hat-function FEM) ---
        M = np.zeros((N, N))
        for i in range(N - 1):
            h = xl[i + 1] - xl[i]
            M[i, i] += 2.0 * h / 6.0
            M[i, i + 1] += h / 6.0
            M[i + 1, i] += h / 6.0
            M[i + 1, i + 1] += 2.0 * h / 6.0
        # Dirichlet BC on last node (match Fortran)
        M[-1, :] = 0.0
        M[:, -1] = 0.0
        M[-1, -1] = 1.0
        self.mass_matrix = M

        # --- Potential matrix Q = int (1/r) phi_l phi_l' dr ---
        Q = np.zeros((N, N))
        for i in range(1, N - 1):
            r_l = xl[i]
            r_lm1 = xl[i - 1]
            r_lp1 = xl[i + 1]

            # Diagonal
            Q[i, i] = 0.5 * (
                (r_l**2 + 2.0 * r_lm1**2 * (np.log(r_l) - np.log(r_lm1))
                 - 4.0 * r_l * r_lm1 + 3.0 * r_lm1**2) / (r_l - r_lm1)**2
                - (r_l**2 + 2.0 * r_lp1**2 * (np.log(r_l) - np.log(r_lp1))
                   - 4.0 * r_l * r_lp1 + 3.0 * r_lp1**2) / (r_l - r_lp1)**2
            )
            # Off-diagonal: l+1
            if i < N - 1:
                Q[i, i + 1] = (
                    -r_l**2
                    + 2.0 * r_l * r_lp1 * (np.log(r_l) - np.log(r_lp1))
                    + r_lp1**2
                ) / (2.0 * (r_l - r_lp1)**2)
            # Off-diagonal: l-1
            if i > 0:
                Q[i, i - 1] = (
                    r_l**2
                    - 2.0 * r_l * r_lm1 * (np.log(r_l) - np.log(r_lm1))
                    - r_lm1**2
                ) / (2.0 * (r_l - r_lm1)**2)
        self.potential_matrix = Q

        # --- hz, hth from q profile (with correct btor sign) ---
        q_xl = self._q_xl
        r_over_qR0 = xl / (q_xl * self.R0)
        btor_sign = np.sign(self.btor)
        hz_xl = btor_sign / np.sqrt(1.0 + r_over_qR0**2)
        hth_xl = r_over_qR0 * hz_xl
        self._hz_xl = hz_xl
        self._hth_xl = hth_xl

        # --- Weighted mass matrix M^{hth} ---
        M_hth = np.zeros((N, N))
        for i in range(N):
            M_hth[i, i] = hth_xl[i] * M[i, i]
            if i < N - 1:
                hth_mid = 0.5 * (hth_xl[i] + hth_xl[i + 1])
                M_hth[i, i + 1] = hth_mid * M[i, i + 1]
                M_hth[i + 1, i] = hth_mid * M[i + 1, i]
        self.M_hth = M_hth

        # --- Weighted potential matrix Q^{hz} ---
        Q_hz = np.zeros((N, N))
        for i in range(N):
            Q_hz[i, i] = hz_xl[i] * Q[i, i]
            if i < N - 1:
                hz_mid = 0.5 * (hz_xl[i] + hz_xl[i + 1])
                Q_hz[i, i + 1] = hz_mid * Q[i, i + 1]
                Q_hz[i + 1, i] = hz_mid * Q[i + 1, i]
        self.Q_hz = Q_hz

    # ------------------------------------------------------------------
    def solve_apar(self, Br_boundary=1.0 + 0j):
        """Solve coupled Poisson-Ampere using A_∥ formulation.

        Solves for (Φ, A_∥) where -∇²_⊥ A_∥ = (4π/c) j_∥,
        then recovers Br = α(r) · A_∥.

        Parameters
        ----------
        Br_boundary : complex
            Br value at the right boundary (default 1.0).
        """
        N = self.N
        r = self.xl

        # α(r) = i·(m/r · h_z - kz · h_θ) converts A_∥ → Br
        alpha = 1j * (self.m_mode / r * self._hz_xl - self.kz * self._hth_xl)
        D_alpha = np.diag(alpha)

        # Build ∇²_⊥ FEM matrix (negative semi-definite, like self.laplace).
        # self.laplace has zeroed boundary couplings (L[1,0]=0, L[N-2,N-1]=0)
        # so we build a fresh matrix with full tridiagonal coupling.
        laplace_perp = np.zeros((N, N))
        for i in range(1, N - 1):
            hi_m = r[i] - r[i - 1]
            hi_p = r[i + 1] - r[i]
            laplace_perp[i, i - 1] = 1.0 / hi_m
            laplace_perp[i, i] = -(1.0 / hi_m + 1.0 / hi_p)
            laplace_perp[i, i + 1] = 1.0 / hi_p

        # Assemble 2N×2N block system for (Φ, A_∥)
        A = np.zeros((2 * N, 2 * N), dtype=complex)
        b = np.zeros(2 * N, dtype=complex)

        # Top-left: Δ + 4π K^{ρΦ}  (same as Br formulation)
        A[:N, :N] = self.laplace + 4.0 * pi * self.k_rho_phi

        # Top-right: 4π K^{ρB} · D_α  (since Br = D_α · A_∥)
        A[:N, N:] = 4.0 * pi * self.k_rho_b @ D_alpha

        # Bottom-left: (4π/c) K^{jΦ}
        A[N:, :N] = (4.0 * pi / sol) * self.k_j_phi

        # Bottom-right: ∇²_⊥ + (4π/c) K^{jB} · D_α
        A[N:, N:] = laplace_perp + (4.0 * pi / sol) * self.k_j_b @ D_alpha

        # --- Boundary conditions ---
        # Φ(left) = 0
        A[0, :] = 0.0
        A[0, 0] = 1.0
        b[0] = 0.0

        # Φ(right) = zero-misalignment BC
        A[N - 1, :] = 0.0
        A[N - 1, N - 1] = 1.0
        b[N - 1] = -1j * self._e0r_xl[-1] * Br_boundary / (
            self._b0_xl[-1] * self._kp_xl[-1])

        # A_∥(left) = 0
        A[N, :] = 0.0
        A[N, N] = 1.0
        b[N] = 0.0

        # A_∥(right) = Br_boundary / α(right)
        A[2 * N - 1, :] = 0.0
        A[2 * N - 1, 2 * N - 1] = 1.0
        b[2 * N - 1] = Br_boundary / alpha[-1]

        # Solve and recover Br
        x = np.linalg.solve(A, b)
        self.phi = x[:N]
        self.apar = x[N:]
        self.br = D_alpha @ self.apar

    # ------------------------------------------------------------------
    def solve_poisson(self, br_prescribed):
        """Solve the Poisson-only (electrostatic) problem.

        Parameters
        ----------
        br_prescribed : ndarray, shape (N,)
            Prescribed Br field on the xl grid.
        """
        a_mat = self.laplace + 4.0 * pi * self.k_rho_phi
        b_vec = -4.0 * pi * (self.k_rho_b @ br_prescribed)

        # Zero-misalignment BCs (same as KIMPoissonSolver)
        a_mat = a_mat.copy()
        b_vec = b_vec.copy()
        N = self.N
        phi_left = -self.phi_aligned[0]
        phi_right = -self.phi_aligned[-1]

        b_vec[1:N - 1] = (
            b_vec[1:N - 1]
            - a_mat[1:N - 1, 0] * phi_left
            - a_mat[1:N - 1, N - 1] * phi_right
        )
        a_mat[:, 0] = 0.0
        a_mat[:, -1] = 0.0
        a_mat[0, :] = 0.0
        a_mat[-1, :] = 0.0
        a_mat[0, 0] = 1.0
        a_mat[-1, -1] = 1.0
        b_vec[0] = phi_left
        b_vec[-1] = phi_right

        self.phi_es = np.linalg.solve(a_mat, b_vec)

    # ------------------------------------------------------------------
    def compute_jpar(self, phi, br):
        """Compute parallel current density j_par = M^{-1}(K_j_phi*Phi + K_j_B*Br).

        Parameters
        ----------
        phi : ndarray, shape (N,)
        br : ndarray, shape (N,)

        Returns
        -------
        ndarray, shape (N,), complex
        """
        M_inv = np.linalg.inv(self.mass_matrix)
        return M_inv @ (self.k_j_phi @ phi + self.k_j_b @ br)

    # ------------------------------------------------------------------
    def compute_rho(self, phi, br):
        """Compute charge density from kernel response.

        Parameters
        ----------
        phi : ndarray, shape (N,)
        br : ndarray, shape (N,)

        Returns
        -------
        ndarray, shape (N,), complex
        """
        return self.k_rho_phi @ phi + self.k_rho_b @ br

    # ------------------------------------------------------------------
    def compute_E_perp_MA(self, phi, br):
        """Compute misalignment perpendicular electric field.

        E_perp = -i * ks * Phi  (from potential)
        E_perp_psi = Er * Br * ks / (B0 * kp)  (from perturbed flux surfaces)
        E_perp_MA = E_perp + E_perp_psi

        Parameters
        ----------
        phi : ndarray, shape (N,)
        br : ndarray, shape (N,)

        Returns
        -------
        dict with keys 'E_perp', 'E_perp_psi', 'E_perp_MA'
        """
        ks_xl = (self.m_mode * self._hz_xl - self.n_mode * self._hth_xl / self.R0) / self.xl
        E_perp = -1j * ks_xl * phi
        E_perp_psi = self._e0r_xl * br * ks_xl / (self._b0_xl * self._kp_xl)
        E_perp_MA = E_perp + E_perp_psi
        return {"E_perp": E_perp, "E_perp_psi": E_perp_psi, "E_perp_MA": E_perp_MA}

    # ------------------------------------------------------------------
    def spline_field(self, field, n_points=5000):
        """Interpolate a complex field onto a fine equidistant grid.

        Parameters
        ----------
        field : ndarray, shape (N,), complex
        n_points : int

        Returns
        -------
        r_fine, field_fine : tuple of ndarray
        """
        r_fine = np.linspace(self.xl[0], self.xl[-1], n_points)
        cs_re = CubicSpline(self.xl, field.real)
        cs_im = CubicSpline(self.xl, field.imag)
        return r_fine, cs_re(r_fine) + 1j * cs_im(r_fine)

    # ------------------------------------------------------------------
    def load_fortran_solution(self):
        """Load the Fortran EM solution from the same HDF5 file for validation.

        Returns
        -------
        dict with keys 'phi', 'br', 'jpar' (or None if not found)
        """
        result = {}
        with h5py.File(self.h5_path, "r") as h5f:
            for key, name in [
                ("fields/Phi", "phi"),
                ("fields/Br_selfconsistent", "br"),
                ("fields/jpar", "jpar"),
            ]:
                if key not in h5f:
                    result[name] = None
                    continue
                data = h5f[key][:]
                if data.dtype.names and "real" in data.dtype.names:
                    result[name] = data["real"] + 1j * data["imag"]
                else:
                    result[name] = data
        return result
