import numpy as np
import os
import h5py
import sys

sys.path.append(os.path.dirname(__file__) + '../../postproc_py_class/')
from utility_class import *

class QL_Balance_interface():

    def __init__(self, run_path, shot, time, name, input_file):
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
        self.input_file = input_file

        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f'Input file {self.input_file} not found.')
        else:
            try:
                h5_input = h5py.File(self.input_file, 'r')
            except:
                raise ValueError(f'Error reading input file {self.input_file}.')

        os.path.makedirs(self.run_path, exist_ok=True)
        self.output_path = os.path.join(self.run_path, 'out')
        os.path.makedirs(self.output_path, exist_ok=True)


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
        h5_input.close()