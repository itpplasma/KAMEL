import os

import h5py
import numpy as np


class KiLCA_data:

    device = ""

    kinetic_profiles = {
        "r_eff": {"value": np.array([]), "unit": "cm"},
        "n_prof": {"value": np.array([]), "unit": "cm^-3"},
        "Te_prof": {"value": np.array([]), "unit": "eV"},
        "Ti_prof": {"value": np.array([]), "unit": "eV"},
        "Vz_prof": {"value": np.array([]), "unit": "cm/s"},
        "Er_prof": {"value": np.array([]), "unit": "statV / cm"},
    }

    equilibrium = {
        "psi_N": {"value": np.array([]), "unit": "1"},
        "Btor": {"value": 0.0, "unit": "G"},
        "R0": {"value": 0.0, "unit": "cm"},
        "q_prof": {"value": np.array([]), "unit": "1"},
    }

    resonance = {
        "r_res": {"value": 0.0, "unit": "cm", "type": "float"},
        "m_modes": {"value": 0, "unit": "1", "type": "int"},
        "n_mode": {"value": 0, "unit": "1", "type": "int"},
    }

    def __init__(self):
        pass

    def add_cylinder_equilibrium(self, dats):
        for key in dats.keys():
            equil_data_keys.append(key)
            setattr(self, key, dats[key])

    def write_to_h5(self, h5path):
        if os.path.exists(h5path):
            self.h5f = h5py.File(h5path, "a")
        else:
            self.h5f = h5py.File(h5path, "w")

        group_names = ["KiLCA/kinetic_profiles", "KiLCA/equilibrium", "KiLCA/resonance"]
        data_structs = [self.kinetic_profiles, self.equilibrium, self.resonance]

        for i, group_name in enumerate(group_names):
            grp = self.h5f.create_group(group_name)
            for key in self.data_structs[i].keys():
                ds = grp.create_dataset(key, data=self.data_structs[i][key]["value"])
                ds.attrs["unit"] = self.data_structs[i][key]["unit"]
        h5f.close()

    def read_from_h5(self, h5path):
        if os.path.exists(h5path):
            self.h5f = h5py.File(h5path, "r")
        else:
            raise ValueError(f"File {h5path} does not exist.")
