import numpy as np
import matplotlib.pyplot as plt
import shutil
import h5py
import os
import sys
import subprocess

sys.path.append(os.path.dirname(__file__) + '/../fieldpy/')
from fieldpy import fieldpy

sys.path.append(os.path.dirname(__file__) + '/../neo2_for_Er/')
from neo2_for_Er import neo2_for_Er

sys.path.append(os.path.dirname(__file__) + '/../profile_processor/')
from profile_processor import profile_processor

sys.path.append(os.path.dirname(__file__) + '/../AnLoBiCr/')
from analytical_local_bif_criterion import *

sys.path.append(os.path.dirname(__file__) + '/../tMHD_current/')
from tMHD_current import *


class KQ_processor:
    """ Kilca QL-balance processor class to calculate equilibrium and profiles for 
        KiLCA and QL-Balance runs.
    """

    def __init__(self, shot, time, runpath):
        
        self.shot = shot
        self.time = time
        self.runpath = runpath
        self.save_profs = os.path.join(self.runpath, 'profiles/')
        self.scriptpath = os.getcwd()

    def process_equilibrium(self, gfile, pfile='', convex_wall='', flux_data='', skip=False):
        """Process the equilibrium EQDSK file with field_divB0 (without perturbations).
        gfile ... path to the EQDSK file
        flux_data ... path where the flux data will be moved to"""

        self.gfile = gfile 
        self.pfile = pfile
        self.convex_wall = convex_wall 
        self.flux_data = flux_data 

        if skip == True:
            print('Equilibrium already processed, skipping...')
            return
        print('Processing equilibrium...')

        self.fp = fieldpy(self.gfile, self.pfile, self.convex_wall, self.flux_data)
        self.fp.write_field_divB0_inp(self.fp.path_to_fourier_modes_exe + 'template_field_divB0.inp', self.fp.path_to_fourier_modes_exe + 'field_divB0.inp')
        self.fp.run_fourier_modes() 

    def process_profiles(self, prof_path, skip=True):
        """Process the kinetic profiles. Map the 4 input profiles (density, electron temperature, 
            ion temperature and toroidal rotation) onto the effective radius determined in the 
            equilibrium processing. Use NEO-2 to determine Er profile.
        """
        
        self.pp = profile_processor(self.runpath)

        self.pp.run_fieldpy(self.gfile, self.convex_wall, self.flux_data, skip=True)

        self.prof_path = prof_path 

        self.pp.map_profs_to_reff(self.prof_path, self.save_profs, self.flux_data, plot=False)

        self.pp.calc_Er_prof(recalc=not skip)

        self.pp.determine_anomalous_diff_coeff()

    def rescale_dens_prof(self, rescale_factor):
        """Rescale the density profile by a constant factor."""

        self.pp.ne = self.pp.ne * rescale_factor

    def rescale_el_temp_prof(self, rescale_factor):
        """Rescale the electron temperature profile by a constant factor."""

        self.pp.Te = self.pp.Te * rescale_factor

    def get_tMHD_current(self, curr_file, m_mode=np.array([10]), delta_phi=0.0, coil_curr_scale_l=1.0, coil_curr_scale_u=1.0):
        """Get the tMHD current for a given m_mode."""

        self.tmhd = tMHD_current()

        self.tmhd.set_equil(self.flux_data + 'equil_r_q_psi.dat', self.flux_data + 'btor_rbig.dat')
        #self.tmhd.load_curr_harmonics_MARSF(curr_file)
        self.tmhd.load_current_MARSF(curr_file)
        self.tmhd.mix_coil_rows(delta_phi=delta_phi, coil_curr_scale_l=coil_curr_scale_l, coil_curr_scale_u=coil_curr_scale_u)
        self.prepare_efit2boozer_inp()
        self.prepare_field_divB0_inp()
        self.tmhd.get_Jpar_over_B0_boozer_harmonics()
        self.tmhd.integrate_curr_dens(m_mode=m_mode)

        return self.tmhd

    def prepare_efit2boozer_inp(self):
        """Prepare the input files for efit2boozer."""

        print(os.path.abspath(os.path.join(os.path.dirname(__file__)+ '/../../../efit_to_boozer/RUN/efit_to_boozer.inp')))
        print(self.scriptpath + '/efit_to_boozer.inp')
        shutil.copy(os.path.join(os.path.dirname(__file__)+ '/../../../efit_to_boozer/RUN/efit_to_boozer.inp'), self.scriptpath + '/efit_to_boozer.inp')
        # get psi max
        with open(self.flux_data + 'equil_r_q_psi.dat', 'r') as f:
            lines = f.readlines()
            psi_max = float(lines[1].split()[4])
        # write psi max in efit file
        with open(self.scriptpath + '/efit_to_boozer.inp', 'a') as f:
            f.write(str(psi_max) + '   psimax - psi of plasma boundary\n')

    def prepare_field_divB0_inp(self):
        """Prepare the input files for field_divB0."""

        fp = fieldpy(self.gfile, self.pfile, self.convex_wall, self.flux_data)
        fp.write_field_divB0_inp(os.path.join(os.path.dirname(__file__)+ '/../../../efit_to_boozer/RUN/field_divB0.inp'), self.scriptpath + '/field_divB0.inp')

    def apply_analytical_criterion(self, m_mode=np.array([6]), n_mode=np.array([2]), Impar=1.0):

        self.crit = np.zeros((len(m_mode), len(n_mode)))

        for i, n in enumerate(n_mode):
            for j, m in enumerate(m_mode):
                #print(analytical_local_criterion(m,n, self.pp.r_eff, self.pp.Da, self.pp.R0, self.pp.Btor, self.pp.r_eff, self.pp.q, self.pp.Te, self.pp.Ti, self.pp.ne, self.pp.Er, Impar))
                self.crit[j,i] = analytical_local_criterion(m,n, self.pp.r_eff, self.pp.Da, self.pp.R0, self.pp.Btor, self.pp.r_eff, self.pp.q, self.pp.Te, self.pp.Ti, self.pp.ne, self.pp.Er, Impar)

        return self.crit

    def do_parameter_scan_analytical_crit(self, parameter, type, values):
        """Do parameter scans of the input profiles.
            input:
                parameter ... string, either ne, Te, Ti or Vz
                type ... string, either 'constant rescale' or , tbi
                values ... numpy array of values, e.g. rescaling values
        """

        if parameter == 'ne':
            pass
    
    def read_equil(self):
        pass

    def read_profiles(self):
        pass

    def prepare_balance_input(self, input_file):
        """"Prepare the input files for the balance code."""
        self.input_file = input_file

        if os.path.exists(self.input_file):
            raise Warning(f'Input file {self.input_file} already exists. Overwriting...')
            os.remove(self.input_file)
        else:
            try:
                h5_input = h5py.File(self.input_file, 'w')
            except:
                raise ValueError(f'Error creating input file {self.input_file}.')

        h5_input.create_dataset('shot', data=self.shot)
        h5_input.create_dataset('time', data=self.time)
        print('Git version: ', self.get_git_version())
        h5_input.create_dataset('git_version', data=self.get_git_version())
        
    def get_git_version():
        try:
            git_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip().decode('utf-8')
            return git_hash
        except subprocess.CalledProcessError:
            return None