# Python framework for analytical bifurcation criterion

import os

import h5py
import numpy as np


class pyABC:

    shot = 0
    time_slice = 0

    m_mode = np.array([])  # poloidal mode numbers
    n_mode = np.array([])  # toroidal mode numbers

    path_equil = ""
    path_tmhd = ""
    path_profiles = ""

    path_run = ""
    path_h5inp = ""

    def __init__(self, shot, time_slice, path_run):
        """Constructor of the python class for the approximate bifurcation criterion."""

        self.path_run = path_run
        os.system("mkdir -p " + self.path_run)

        self.shot = shot
        self.time_slice = time_slice

        self.path_h5inp = (
            self.path_run + str(self.shot) + "_" + str(self.time_slice) + "_a2bc_inp.hdf5"
        )

        self.h5inp = h5py.File(self.path_h5inp)
        shot_ds = self.h5inp.create_dataset("shot", (1,), dtype="i")
        shot_ds[0] = self.shot
        time_slice_ds = self.h5inp.create_dataset("time_slice", (1,), dtype="i")
        time_slice_ds[0] = self.time_slice
        self.h5inp.close()

    def set_modes(self, m, n):
        if ~isinstance(m, np.ndarray):
            raise ValueError("Poloidal mode numbers need to be numpy arrays.")
        if isinstance(n, np.ndarray):
            if len(n) != len(m):
                raise ValueError("Mode number arrays must be of same size.")

        self.m_mode = m
        self.n_mode = n

    def set_coil(self, cfile):
        self.coil_file = cfile

    def set_profiles(self, path_to_profiles, type="PED_MMARKL_rho_pol.dat"):
        self.path_profiles = path_to_profiles
        # TODO: write profiles to h5 file
        # expt profiles
        expt_prof = ["ne", "Te", "Ti", "Vz"]

        with h5py.File(self.path_h5inp, "a") as f:
            for prof in expt_prof:
                # load profile data from text file
                dat = np.loadtxt(
                    self.path_profiles
                    + str(int(self.shot))
                    + "."
                    + str(int(self.time_slice))
                    + prof
                    + "_"
                    + type
                )
                # write to hdf5
                dset = f.create_dataset("/profiles/" + prof, data=dat)

    def set_equi(self, gfile, flux_data_path):
        self.file_equi = gfile
        self.flux_data_path = flux_data_path

    def set_Da_estimation(self):
        pass

    def calc_ABC(self, h5inp):
        """Calculate analytical bifurcation criterion.
        Input:
            h5inp ... h5 input file containing required data."""
        pass
