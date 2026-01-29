import os
from typing import Dict, Optional

import h5py
import numpy as np
from numpy.typing import NDArray


class KIMData:

    data_types = ["dat", "h5"]  # possible data types of the KIM output files

    def __init__(
        self, data_path: str, m=-6, n=2, collision_model="FokkerPlanck", data_type="dat"
    ) -> None:
        self.r_is_set = False
        self.data_path = data_path
        self.mode_string = f"m{m}_n{n}"
        self.mpol = m
        self.ntor = n
        self.collision_model = collision_model
        self.data_type = data_type

        if self.collision_model == "FokkerPlanck":
            self.collision_model_acronym = "fp"

        if self.data_type == "h5":
            if not os.path.isfile(self.data_path):
                raise FileNotFoundError(f"H5 file not found: {self.data_path}")
            if not self.data_path.endswith(".h5"):
                raise ValueError(f"Invalid H5 file: {self.data_path}")

            self.h5f = h5py.File(self.data_path, "r")

        # main data dictionary
        self.data: Dict[str, Optional[NDArray[np.float64]]] = {
            "r": None,
            "Phi_m": None,
            "Phi_m_e": None,
            "Phi_m_i": None,
            "B0": None,
        }

    def __del__(self):
        self.r_is_set = False
        self.data = {}

    def get_Phi_m(self) -> None:
        """
        Read the electrostatic potential perturbation $\Phi_{\bm}$,
        which is the solution to the Poisson equation in KIM.
        """
        if self.data_type == "h5":
            r = self.h5f[f"grid/xl_xb"][:]

            phi = np.array(self.h5f[f"fields/Phi_m"][:])
            self.set_profile(r, phi["real"] + 1j * phi["imag"], "Phi_m")

            phi_e = self.h5f[f"/fields/Phi_m_e"][:]
            self.set_profile(r, phi_e["real"] + 1j * phi_e["imag"], "Phi_m_e")

            phi_i = self.h5f[f"/fields/Phi_m_i"][:]
            self.set_profile(r, phi_i["real"] + 1j * phi_i["imag"], "Phi_m_i")

            return

        path = os.path.join(self.data_path, self.mode_string, "fields", "Phi_m.dat")
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], "Phi_m")
        except FileNotFoundError:
            print("Phi_m file not found in: " + path)

        path = os.path.join(self.data_path, self.mode_string, "fields", "Phi_m_e.dat")
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], "Phi_m_e")
        except FileNotFoundError:
            print("Phi_m file not found in: " + path)

        path = os.path.join(self.data_path, self.mode_string, "fields", "Phi_m_i.dat")
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], "Phi_m_i")
        except FileNotFoundError:
            print("Phi_m file not found in: " + path)

    def get_jpar(self) -> None:
        """
        Read the electrostatic potential perturbation $j_{\parallel}$,
        which is the solution to the Poisson equation in KIM.
        """

        if self.data_type == "h5":
            r = self.h5f[f"grid/xl_xb"][:]

            jpar = self.h5f[f"fields/jpar"][:]
            self.set_profile(r, jpar["real"] + 1j * jpar["imag"], "jpar")

            jpar_e = self.h5f[f"/fields/jpar_e"][:]
            self.set_profile(r, jpar_e["real"] + 1j * jpar_e["imag"], "jpar_e")

            jpar_i = self.h5f[f"/fields/jpar_i"][:]
            self.set_profile(r, jpar_i["real"] + 1j * jpar_i["imag"], "jpar_i")

            return

        path = os.path.join(self.data_path, self.mode_string, "fields", "jpar.dat")
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], "jpar")
        except FileNotFoundError:
            print("jpar file not found in: " + path)

        path = os.path.join(self.data_path, self.mode_string, "fields", "jpar_e.dat")
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], "jpar_e")
        except FileNotFoundError:
            print("jpar file not found in: " + path)

        path = os.path.join(self.data_path, self.mode_string, "fields", "jpar_i.dat")
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], "jpar_i")
        except FileNotFoundError:
            print("jpar file not found in: " + path)

    def get_Phi_aligned(self) -> None:
        """
        Read the aligned electrostatic potential $\Phi_{A}$,
        which is used as the boundary condition in the Poisson equation in KIM.
        """

        if self.data_type == "h5":
            r = self.h5f[f"grid/xl_xb"][:]

            phi_A = self.h5f[f"fields/Phi_aligned"][:]
            self.set_profile(r, phi_A["real"] + 1j * phi_A["imag"], "Phi_aligned")

            return

        path = os.path.join(self.data_path, self.mode_string, "fields", "Phi_aligned.dat")
        try:
            phi_A = np.loadtxt(path)
            self.set_profile(phi_A[:, 0], phi_A[:, 1] + 1j * phi_A[:, 2], "Phi_aligned")
        except FileNotFoundError:
            print("Phi_aligned file not found in: " + path)

    def get_Phi_MA(self) -> None:
        """
        Read the misaligned electrostatic potential $\Phi_{MA}$,
        which is used as the boundary condition in the Poisson equation in KIM.
        """

        if self.data_type == "h5":
            r = self.h5f[f"grid/xl_xb"][:]

            phi_MA = self.h5f[f"fields/Phi_MA"][:]
            self.set_profile(r, phi_MA["real"] + 1j * phi_MA["imag"], "Phi_MA")

            return

        path = os.path.join(self.data_path, self.mode_string, "fields", "Phi_MA.dat")
        try:
            phi_MA = np.loadtxt(path)
            self.set_profile(phi_MA[:, 0], phi_MA[:, 1] + 1j * phi_MA[:, 2], "Phi_MA")
        except FileNotFoundError:
            print("Phi_MA file not found in: " + path)

    def get_B0(self) -> None:

        if self.data_type == "h5":
            r = self.h5f[f"grid/xl_xb"][:]

            B0 = self.h5f[f"backs/B0"][:]
            self.set_profile(r, B0, "B0")

            return

        path = os.path.join(self.data_path, self.mode_string, "backs", "B0.dat")
        try:
            B0 = np.loadtxt(path)
            self.set_profile(B0[:, 0], B0[:, 1], "B0")
        except FileNotFoundError:
            print("B0 file not found in: " + path)

    def normalize_by_B0(self, dict_entry: str) -> None:
        if self.data["B0"] is not None:
            self.data[dict_entry + "/B0"] = self.data[dict_entry] / self.data["B0"]
        else:
            print("Envoke get_B0 first to normalize profile by B0")

    def set_r(self, r: NDArray[np.float64]) -> None:
        if self.data["r"] is None:
            self.data["r"] = np.asarray(r)
            self.r_is_set = True

    def set_profile(
        self, r: NDArray[np.float64], quant: NDArray[np.float64], quant_string: str
    ) -> None:
        if not self.r_is_set:
            self.set_r(r)
        if (np.asarray(r) == self.data["r"]).all():
            self.data[quant_string] = quant
        else:
            self.data[quant_string] = np.interp(self.data["r"], r, quant)
