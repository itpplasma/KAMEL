import os
import shutil
import sys

import h5py
import numpy as np
from KiLCA_interface import KiLCA_interface
from KiLCA_interface.KiLCA_postprocessor import KiLCA_postprocessor
from utility import utility

from .balance_conf import *
from .balance_input_h5 import *


class QL_Balance_interface:

    machine = "AUG"  # default machine is AUG
    executable_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "build", "install", "bin", "ql-balance.x"
    )

    run_types = ["SingleStep", "TimeEvolution", "ParameterScan"]

    delta_r_antenna = 3.0  # distance of RMP antenna from plasma boundary

    def __init__(self, run_path, shot, time, name, input_file="", debug=True):
        """Constructor of the QL-Balance interface class.
        Args:
            run_path (str): Path to the run directory.
            shot (int): Shot number.
            time (float): Time of the shot.
            name (str): Name of the balance run, should indicate purpose.
            input_file (str): Name of the input file (hdf5).
        """

        self.run_path = run_path
        self.shot = shot
        self.time = time
        self.name = name
        self.debug = debug

        # if not os.path.exists(self.input_h5_file):
        # raise FileNotFoundError(f'Input file {self.input_h5_file} not found.')
        # else:
        # try:
        # h5_input = h5py.File(self.input_h5_file, 'r')
        # except:
        # raise ValueError(f'Error reading input file {self.input_file}.')

        os.makedirs(self.run_path, exist_ok=True)
        self.output_path = os.path.join(self.run_path, "out")
        os.makedirs(self.output_path, exist_ok=True)
        self.output_h5_file = os.path.join(
            self.output_path, f"{self.shot}_{self.time}_{self.name}.hdf5"
        )

        if input_file == "":
            self.input_h5_file = os.path.join(
                self.run_path, f"INPUT_{self.shot}_{self.time}_{self.name}.hdf5"
            )
        else:
            self.input_h5_file = input_file
        self.util = utility()

    def set_type_of_run(self, run_type="SingleStep"):
        """Set the type of the run, e.g. 'SingleStep', 'TimeEvolution' or 'ParameterScan"""
        assert run_type in self.run_types, f"Run type {run_type} not supported."
        self.run_type = run_type
        try:
            self.conf.conf["balancenml"]["type_of_run"] = run_type
        except:
            raise ValueError(f"Namelist config not read to write run type {run_type}.")

    def set_modes(self, m_mode, n_mode):
        self.m_mode = m_mode
        self.n_mode = n_mode

    def copy_profiles(self, profile_path):
        """Copy the profiles to the run directory."""
        os.makedirs(os.path.join(self.run_path, "profiles"), exist_ok=True)
        files = [
            f for f in os.listdir(profile_path) if os.path.isfile(os.path.join(profile_path, f))
        ]
        for f in files:
            shutil.copy2(os.path.join(profile_path, f), os.path.join(self.run_path, "profiles"))

    def link_profiles(self, profile_path):
        if not os.path.exists(profile_path):
            raise FileNotFoundError(f"Profile path {profile_path} not found.")
        # Resolve symlinks to prevent circular references
        resolved_profile_path = os.path.normpath(os.path.realpath(profile_path))
        resolved_run_path = os.path.normpath(os.path.realpath(self.run_path))
        # Check if profile_path is inside run_path (would create circular symlink)
        if (
            resolved_profile_path.startswith(resolved_run_path + os.sep)
            or resolved_profile_path == resolved_run_path
        ):
            raise ValueError(
                f"Profile path {profile_path} is inside run_path {self.run_path}. "
                f"This would create a circular symlink. Use an external profile directory."
            )
        link_path = os.path.join(self.run_path, "profiles")
        # Remove existing symlink, file, or directory at link_path
        if os.path.islink(link_path):
            os.unlink(link_path)
        elif os.path.isdir(link_path):
            # Check if directory is empty or only contains a 'profiles' symlink (from previous bug)
            try:
                contents = os.listdir(link_path)
                if contents == ["profiles"] and os.path.islink(os.path.join(link_path, "profiles")):
                    # Remove the buggy nested symlink and the directory
                    os.unlink(os.path.join(link_path, "profiles"))
                    os.rmdir(link_path)
                elif len(contents) == 0:
                    os.rmdir(link_path)
                else:
                    raise ValueError(
                        f"Cannot create profiles symlink: {link_path} is a non-empty directory. "
                        f"Please remove or rename it manually."
                    )
            except FileNotFoundError:
                # Directory disappeared between checks; nothing left to clean up
                pass
            except OSError as e:
                raise ValueError(
                    f"Cannot safely prepare profiles symlink at {link_path}: {e}"
                ) from e
        elif os.path.exists(link_path):
            os.unlink(link_path)
        os.symlink(resolved_profile_path, link_path)

    def set_factors(self, fac_n, fac_Te, fac_Ti, fac_vz):
        """Set the factors for the balance run."""
        self.facs = {"fac_n": fac_n, "fac_Te": fac_Te, "fac_Ti": fac_Ti, "fac_vz": fac_vz}

    def prepare_balance(self, Btor, a_minor):
        self.prepare_KiLCA(Btor, a_minor)
        # self.prepare_balance_input(self.input_h5_file)
        if self.debug:
            print("D: Prepare balance output")
        self.prepare_balance_output(self.output_h5_file)
        self.link_executable()

    def set_existing_balance_input(self, input_file):
        self.input_h5_file = input_file

    def prepare_balance_input(self, input_file):
        """ "Prepare the input files for the balance code."""
        self.input_h5_file = input_file

        if os.path.exists(self.input_h5_file):
            print(f"Input file {self.input_h5_file} already exists. Overwriting...")
            os.remove(self.input_h5_file)

        try:
            h5_input = h5py.File(self.input_h5_file, "w")
        except:
            raise ValueError(f"Error creating input file {self.input_h5_file}.")

        h5_input.create_dataset("shot", data=self.shot)
        h5_input.create_dataset("time", data=self.time)

        if self.debug:
            print("D: Git version: ", self.util.get_git_version())
        h5_input.create_dataset("git_version", data=self.util.get_git_version())
        h5_input.close()

    def prepare_balance_output(self, output_file):
        self.output_h5_file = output_file
        self.check_if_factors_set()
        self.prepare_input_h5()
        self.prepare_output_h5()

    def prepare_output_h5(self):
        h5f = h5py.File(self.output_h5_file, "w")
        self.input_h5.write_fac_with_bound_info(h5f, "/factors/fac_n", data=[self.facs["fac_n"]])
        self.input_h5.write_fac_with_bound_info(h5f, "/factors/fac_Te", data=[self.facs["fac_Te"]])
        self.input_h5.write_fac_with_bound_info(h5f, "/factors/fac_Ti", data=[self.facs["fac_Ti"]])
        self.input_h5.write_fac_with_bound_info(h5f, "/factors/fac_vz", data=[self.facs["fac_vz"]])
        h5f.close()

    def check_if_factors_set(self):
        if not hasattr(self, "facs"):
            print("No factors set, will use ones.")
            self.facs = {
                "fac_n": np.array([1.0]),
                "fac_Te": np.array([1.0]),
                "fac_Ti": np.array([1.0]),
                "fac_vz": np.array([1.0]),
            }

    def prepare_KiLCA(self, Btor, a_minor):
        if self.debug:
            print("D: Prepare KiLCA")
        self.kil_flre = KiLCA_interface(self.shot, self.time, self.run_path, "flre", self.machine)
        self.kil_flre.set_machine(delta_r_antenna=self.delta_r_antenna, a_minor=a_minor)
        self.kil_flre.background.data["Btor"] = Btor

        self.kil_flre.set_modes(self.m_mode, self.n_mode)
        self.kil_flre.antenna.data["flab"] = [1.0, 0.0]
        self.kil_flre.write()
        self.kil_flre.run()
        self.kil_flre_post = KiLCA_postprocessor(self.kil_flre)
        self.I_KiLCA = self.get_KiLCA_current()

        kil_vac = KiLCA_interface(self.shot, self.time, self.run_path, "vacuum", self.machine)
        kil_vac.background.data["Btor"] = Btor
        kil_vac.set_machine(delta_r_antenna=self.delta_r_antenna, a_minor=a_minor)
        kil_vac.set_modes(self.m_mode, self.n_mode)
        kil_vac.antenna.data["flab"] = [1.0, 0.0]
        kil_vac.write()
        kil_vac.run()
        if self.debug:
            print("D: Finished KiLCA preperation")

    def get_KiLCA_current(self):
        if not hasattr(self, "kil_flre"):
            raise ValueError("KiLCA not prepared.")
        linear_data_path = os.path.join(
            self.run_path, "flre", "linear-data", f"m_{self.m_mode}_n_{self.n_mode}_flab_[1,0]", ""
        )
        background_data_path = os.path.join(self.run_path, "flre", "background-data", "")
        self.r_kilca, self.jpar_kilca = self.kil_flre_post.calculate_parallel_current_density(
            self.m_mode, self.n_mode, linear_data_path, background_data_path
        )
        self.layer_width = self.kil_flre_post.calculate_layer_width(self.m_mode, self.n_mode)
        return self.kil_flre_post.integrate_par_current_dens()

    def link_executable(self):
        """Link the executable to the run directory."""
        self.executable = os.path.join(self.run_path, "ql-balance")
        if not os.path.exists(self.executable_path):
            raise FileNotFoundError(f"Executable {self.executable_path} not found.")
        else:
            try:
                os.unlink(self.executable)
            except:
                pass
            os.symlink(self.executable_path, self.executable)

    def set_default_config_nml(self):
        """Set the default balance configuration."""
        self.read_config_nml()
        self.write_config_nml(os.path.join(self.run_path, "balance_conf.nml"))

    def read_config_nml(self, path=""):
        """Read the balance configuration file."""
        self.conf = balance_conf(path=path)

    def set_config_nml(self):
        self.conf.conf["balancenml"]["flre_path"] = os.path.join(
            os.path.abspath(self.run_path), "flre/"
        )
        self.conf.conf["balancenml"]["vac_path"] = os.path.join(
            os.path.abspath(self.run_path), "vacuum/"
        )
        self.conf.conf["balancenml"]["path2out"] = os.path.abspath(self.output_h5_file)
        self.conf.conf["balancenml"]["path2inp"] = os.path.abspath(self.input_h5_file)

    def write_config_nml(self, path):
        """Write the balance configuration file."""
        self.conf.write_conf(path)

    def prepare_input_h5(self):
        self.input_h5 = Balance_Input_h5(
            self.input_h5_file, os.path.join(self.run_path, "profiles/")
        )
        self.input_h5.get_required_data()
        self.input_h5.write_data_to_h5(self.input_h5_file, self.facs)
        self.write_KiLCA_data_to_input_h5()

    def write_KiLCA_data_to_input_h5(self):
        h5f = h5py.File(self.input_h5_file, "a")
        grp = h5f.create_group("KiLCA")
        grp.create_dataset("I_KiLCA", data=[self.I_KiLCA])
        grp.create_dataset("r", data=self.r_kilca)
        grp.create_dataset("jpar", data=self.jpar_kilca)
        grp.create_dataset("layer_width", data=[self.layer_width])
        h5f.close()

    def run_balance(self, suppress_console_output=True):
        """Run the balance code."""
        print("")
        print(f"== Start balance run {self.name} ==")
        if suppress_console_output:
            options = ">/dev/null 2>&1"
        else:
            options = ""
        cwd = os.getcwd()
        os.chdir(self.run_path)
        out = os.system(f"./ql-balance | tee out/balance.log {options}")
        os.chdir(cwd)
        print(f"== Balance run {self.name} finished. ==")
        print("")
