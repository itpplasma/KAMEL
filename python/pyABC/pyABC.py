# Python framework for analytical bifurcation criterion

# Need:
#   - field div B0 run to calculate q profile
#   - NEO-2 run for k used in the Er calculation
#   - GPEC run for shielding current
#   - Da estimation from input power and radiation power

import numpy as np
import h5py
import os
import sys
# append path to kilca_interface python class, should be one level above
sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'kilca_interface')))
from InpOut import InpOut
#
from GPEC_interface import *
from field_divB0 import *
from fouriermodes import *

class pyABC:

    shot = 0
    time_slice = 0

    m_mode = np.array([]) # poloidal mode numbers
    n_mode = np.array([]) # toroidal mode numbers

    path_equil      = ''
    path_tmhd       = ''
    path_profiles   = ''
    coil_file       = ''

    path_run = ''
    path_h5inp = ''

    lib_pyABC = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
    fourier_path = os.path.abspath(os.path.join(lib_pyABC, '../matlab/fourier/')) + '/'
    convex_path = '/proj/plasma/RMP/DATA/convexwall.dat'
    blueprints_path = os.path.abspath(os.path.join(lib_pyABC, '../matlab/blueprints/')) + '/'


    def __init__(self, shot, time_slice, path_run, name):
        """Constructor of the python class for the approximate bifurcation criterion.
        Input:
        - shot       ... shot number
        - time_slice ... time of the time slice in ms
        - path_run   ... path of run of the balance code
        - name       ... name of the run"""

        self.path_run = path_run
        os.system('mkdir -p ' + self.path_run)
        
        self.shot = shot
        self.time_slice = time_slice
        self.name = name

        self.path_h5inp = self.path_run + str(self.shot) + '_' + str(self.time_slice) + '_a2bc_inp.hdf5'

        self.h5inp = h5py.File(self.path_h5inp, 'w')
        shot_ds = self.h5inp.create_dataset('shot', (1,), dtype='i')
        shot_ds[0] = self.shot
        time_slice_ds = self.h5inp.create_dataset('time_slice', (1,), dtype='i')
        time_slice_ds[0] = self.time_slice
        self.h5inp.close()

    def set_modes(self, m, n):
        '''Set the mode numbers. Also checks if mode numbers have same length.
        Input:
        - m ... poloidal mode number(s), must be numpy array
        - n ... toroidal mode number(s), must be numpy array'''

        if ~isinstance(m, np.ndarray):
            raise ValueError('Poloidal mode numbers need to be numpy arrays.')
        if isinstance(n, np.ndarray):
            if len(n)!=len(m):
                raise ValueError('Mode number arrays must be of same size.')

        self.m_mode = m
        self.n_mode = n

    def set_coil(self, cfile):
        '''Specify the coil file.
        Input:
        - cfile ... path to the coil file'''

        self.coil_file = cfile

    def set_profiles(self, path_to_profiles, type='PED_MMARKL_rho_pol.dat'):
        '''Set profile paths.
        Input:
        - path_to_profiles ... location of the profile text files
        - type             ... type of profiles, i.e. the extension in the file name'''

        self.path_profiles = path_to_profiles
        # TODO: write profiles to h5 file
        # expt profiles
        expt_prof = ['ne', 'Te', 'Ti', 'Vz']

        with h5py.File(self.path_h5inp, 'a') as f:
            for prof in expt_prof:
                dat = np.loadtxt(self.path_profiles + str(int(self.shot)) + '.' + str(int(self.time_slice)) + prof + '_' + type)
                dset = f.create_dataset('/profiles/' + prof, data=dat)

    def run_NEO2(self):
        """Run NEO-2 to calculate k coefficient needed for poloidal rotation velocity
        in the calculation of the E_0r profile."""
        pass


    def set_equi(self, gfile, flux_data_path):
        '''Set the equilibrium file and run field_divB0.
        Input:
        - gfile          ... path to the equilibrium file (gfile format)
        - flux_data_path ... path of the flux data'''

        self.file_equi = gfile
        self.flux_data_path = flux_data_path # path where the equi data should be put

        if (not os.path.exists(self.file_equi)): raise ValueError('Equi file not found in ' + self.file_equi)
        if (not os.path.exists(self.coil_file)): raise ValueError('Coil file not found in ' + self.coil_file)
        if (not os.path.exists(self.convex_path)): raise ValueError('Convex file not found in ' + self.convex_path)

        f0 = field_divB0(self.file_equi, self.coil_file, self.convex_path, self.flux_data_path) # TODO: implementation of this class
        f0.write(self.blueprints_path + f0.BLUEPRINT, self.fourier_path) # TODO: implementation of this method

        # run
        fouriermodes(self.fourier_path, self.flux_data_path)

        raw = np.loadtxt(self.flux_data_path + 'btor_rbig.dat')
        self.b_tor = raw[0]
        self.r_big = raw[1]

        raw = np.loadtxt(self.flux_data_path + 'equil_r_q_psi.dat', skiprow=3)
        self.r_sep_real = raw(-1,0)
        # TODO: get q from equil_r_q_psi


    def run_TMHD_code(self, path_tmhd_run, type='GPEC', flag_force_tmhd = False):
        '''Run the toroidal MHD code, e.g. GPEC (could be MEPHIT in the future).
        Input:
        - path_tmhd_run   ... path where the TMHD code will be run
        - type            ... TMHD code, only GPEC so far
        - flag_force_tmhd ... force the run of the TMHD code'''

        self.path_tmhd = path_tmhd_run
        # relevant output file:
        self.file_tmhd = self.path_tmhd + 'gpec_profile_output_n2.nc'

        # run if forced or if file does not exist
        if flag_force_tmhd or not os.path.exists(self.file_tmhd):
            run_tmhd = True

 
        if type == 'GPEC':
            
            if run_tmhd is True:
               pass
               self.tmhd = GPEC_interface(self.path_tmhd, self.shot, self.time_slice, self.name)
               self.tmhd.setEqui(self.file_equi)
               self.tmhd.setCoil(self.coil_file)
               self.tmhd.write()
               self.tmhd.run()
        else:
            raise ValueError('TMHD code type not implemented')

    def set_Da_estimation(self):
        pass


    def calc_ABC(self, h5inp):
        """Calculate analytical bifurcation criterion.
            Input:
                h5inp ... h5 input file containing required data."""
        pass