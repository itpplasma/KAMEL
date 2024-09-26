import numpy as np
import os
import h5py
import sys
import shutil

from postproc_class import utility_class

from KiLCA_interface import *

from .balance_conf import *
from .balance_input_h5 import *

class QL_Balance_interface():

    machine = 'AUG' # default machine is AUG
    executable_path = os.path.join(os.path.dirname(__file__) + '/../../ql-balance/build/ql-balance')

    run_types = ['SingleStep', 'TimeEvolution', 'ParameterScan']

    def __init__(self, run_path, shot, time, name, input_file=""):
        ''' Constructor of the QL-Balance interface class.
        Args:
            run_path (str): Path to the run directory.
            shot (int): Shot number.
            time (float): Time of the shot.
            name (str): Name of the balance run, should indicate purpose.
            input_file (str): Name of the input file (hdf5).
        '''

        self.run_path = run_path
        self.shot = shot
        self.time = time
        self.name = name 

        #if not os.path.exists(self.input_h5_file):
            #raise FileNotFoundError(f'Input file {self.input_h5_file} not found.')
        #else:
            #try:
                #h5_input = h5py.File(self.input_h5_file, 'r')
            #except:
                #raise ValueError(f'Error reading input file {self.input_file}.')

        os.makedirs(self.run_path, exist_ok=True)
        self.output_path = os.path.join(self.run_path, 'out')
        os.makedirs(self.output_path, exist_ok=True)
        self.output_h5_file = os.path.join(self.output_path, f'{self.shot}_{self.time}_{self.name}.hdf5')

        if input_file == "":
            self.input_h5_file = os.path.join(self.run_path, f'INPUT_{self.shot}_{self.time}_{self.name}.hdf5')
        else:
            self.input_h5_file = input_file
        self.util = utility_class.utility()

    def set_type_of_run(self, run_type='SingleStep'):
        '''Set the type of the run, e.g. 'SingleStep', 'TimeEvolution' or 'ParameterScan'''
        assert run_type in self.run_types, f'Run type {run_type} not supported.'
        self.run_type = run_type
        
    def set_modes(self, m_mode, n_mode):
        self.m_mode = m_mode
        self.n_mode = n_mode
    
    def copy_profiles(self, profile_path):
        '''Copy the profiles to the run directory.'''
        os.makedirs(os.path.join(self.run_path, 'profiles'), exist_ok=True)
        files = [f for f in os.listdir(profile_path) if os.path.isfile(os.path.join(profile_path, f))]
        for f in files:
            shutil.copy2(profile_path + f, os.path.join(self.run_path, 'profiles'))

    def link_profiles(self, profile_path):
        if not os.path.exists(profile_path):
            raise FileNotFoundError(f'Profile path {profile_path} not found.')
        else:
            try:
                os.unlink(self.run_path + 'profiles')
            except:
                pass
            os.symlink(profile_path, self.run_path + 'profiles')
           
    def set_factors(self, fac_n, fac_Te, fac_Ti, fac_vz):
        """Set the factors for the balance run."""
        self.fac_n = fac_n
        self.fac_Te = fac_Te
        self.fac_Ti = fac_Ti
        self.fac_vz = fac_vz
        
    def prepare_balance(self, Btor, a_minor):
        self.prepare_KiLCA(Btor, a_minor)
        #self.prepare_balance_input(self.input_h5_file)
        self.prepare_balance_output(self.output_h5_file)
        self.link_executable()

    def set_existing_balance_input(self, input_file):
        self.input_h5_file = input_file
        
    def prepare_balance_input(self, input_file):
        """"Prepare the input files for the balance code."""
        self.input_h5_file = input_file

        if os.path.exists(self.input_h5_file):
            print(f'Input file {self.input_h5_file} already exists. Overwriting...')
            os.remove(self.input_h5_file)
        
        try:
            h5_input = h5py.File(self.input_h5_file, 'w')
        except:
            raise ValueError(f'Error creating input file {self.input_h5_file}.')

        h5_input.create_dataset('shot', data=self.shot)
        h5_input.create_dataset('time', data=self.time)
        
        print('Git version: ', self.util.get_git_version())
        h5_input.create_dataset('git_version', data=self.util.get_git_version())
        h5_input.close()
    
    def prepare_balance_output(self, output_file):
        self.output_h5_file = output_file
        self.check_if_factors_set()
        self.write_factors_to_h5(self.output_h5_file)
        self.prepare_input_h5()
    
    def check_if_factors_set(self):
        if not hasattr(self, 'fac_n'):
            print('No factors set, will use ones.')
            self.fac_n = np.array([1.0])
            self.fac_Te = np.array([1.0])
            self.fac_Ti = np.array([1.0])
            self.fac_vz = np.array([1.0])
        
    def write_factors_to_h5(self, h5_file):
        f = h5py.File(h5_file, 'w')
        ds = f.create_dataset('/factors/fac_n', data=self.fac_n)
        ds.attrs['lbounds'] = 0
        ds.attrs['ubounds'] = len(self.fac_n)

        ds = f.create_dataset('/factors/fac_Te', data=self.fac_Te)
        ds.attrs['lbounds'] = 0
        ds.attrs['ubounds'] = len(self.fac_Te)

        ds = f.create_dataset('/factors/fac_Ti', data=self.fac_Ti)
        ds.attrs['lbounds'] = 0
        ds.attrs['ubounds'] = len(self.fac_Ti)

        ds = f.create_dataset('/factors/fac_vz', data=self.fac_vz)
        ds.attrs['lbounds'] = 0
        ds.attrs['ubounds'] = len(self.fac_vz)

        f.close()
    
    def prepare_KiLCA(self, Btor, a_minor):
        kil = KiLCA_interface(self.shot, self.time, self.run_path, 'flre', self.machine)
        kil.background.data['Btor'] = Btor
        kil.a_minor = a_minor
        kil.set_machine()
        kil.set_modes(self.m_mode,self.n_mode)
        kil.antenna.data['flab'] = [1.0, 0.0]
        kil.write()
        self.I_KiLCA = self.get_KiLCA_current()
        kil = KiLCA_interface(self.shot, self.time, self.run_path, 'vacuum', self.machine)
        kil.background.data['Btor'] = Btor
        kil.a_minor = a_minor
        kil.set_machine()
        kil.set_modes(self.m_mode,self.n_mode)
        kil.antenna.data['flab'] = [1.0, 0.0]
        kil.write()

    def get_KiLCA_current(self):
        if not hasattr(self, kil):
            raise ValueError('KiLCA not prepared.')
        kil.calculate_parallel_current_density(self.m_mode, self.n_mode, self.run_path + f'flre/linear-data/m_{self.m_mode}_n_{self.n_mode}_flab_[1,0]', self.run_path + f'flre/background-data/')
        kil.calculate_layer_width()
        return kil.integrate_par_current_dens()
    
    def link_executable(self):
        """Link the executable to the run directory."""
        self.executable = os.path.join(self.run_path, 'ql-balance')
        if not os.path.exists(self.executable_path):
            raise FileNotFoundError(f'Executable {self.executable_path} not found.')
        else:
            try:
                os.unlink(self.executable)
            except:
                pass
            os.symlink(self.executable_path, self.executable)
    
    def set_default_config_nml(self):
        """Set the default balance configuration."""
        self.read_config_nml()
        self.write_config_nml(os.path.join(self.run_path, 'balance_conf.nml'))

    def read_config_nml(self, path=""):
        """Read the balance configuration file."""
        self.conf = balance_conf(path=path)
    
    def set_config_nml(self):
        self.conf.conf['balancenml']['flre_path'] = os.path.join(os.path.abspath(self.run_path), 'flre/')
        self.conf.conf['balancenml']['vac_path'] = os.path.join(os.path.abspath(self.run_path), 'vacuum/')
        self.conf.conf['balancenml']['path2out'] = os.path.abspath(self.output_h5_file)
        self.conf.conf['balancenml']['path2inp'] = os.path.abspath(self.input_h5_file)
    
    def write_config_nml(self, path):
        """Write the balance configuration file."""
        self.conf.write_conf(path)
    
    def prepare_input_h5(self):
        self.input_h5 = Balance_Input_h5(self.input_h5_file, os.path.join(self.run_path, 'profiles/'))
        self.input_h5.get_required_data()
        self.input_h5.write_data_to_h5(self.input_h5_file)
    
    def run_balance(self, suppress_console_output=True):
        """Run the balance code."""
        if suppress_console_output:
            options = '>/dev/null 2>&1'
        else:
            options = ''
        os.chdir(self.run_path)
        out = os.system(f'./ql-balance | tee out/balance.log {options}')
        os.chdir(os.path.dirname(__file__))
        print(f'Balance run {self.name} finished.')
    