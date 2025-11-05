import numpy as np
from numpy.typing import NDArray
from typing import Dict, Optional
import os

class KIMData:

    def __init__(self, data_path: str, m=-6, n=2, collision_model='FokkerPlanck') -> None:
        self.r_is_set = False
        self.data_path = data_path
        self.mode_string = f'm{m}_n{n}'
        self.collision_model = collision_model

        if self.collision_model == 'FokkerPlanck':
            self.collision_model_acronym = 'fp'

        # main data dictionary
        self.data : Dict[str, Optional[NDArray[np.float64]]] = \
            {
                'r': None,
                'Phi_m': None,
                'B0': None
            }

    def __del__(self):
        self.r_is_set = False
        self.data = {}

    def get_Phi_m(self) -> None:
        '''
           Read the electrostatic potential perturbation $\Phi_{\bm}$, 
           which is the solution to the Poisson equation in KIM.
        '''
        path = os.path.join(self.data_path, self.mode_string, 'fields', 'Phi_m.dat')
        try:
            phi = np.loadtxt(path)
            self.set_profile(phi[:, 0], phi[:, 1] + 1j * phi[:, 2], 'Phi_m')
        except FileNotFoundError:
            print('Phi_m file not found in: ' + path)

    def get_Phi_aligned(self) -> None:
        '''
           Read the aligned electrostatic potential $\Phi_{A}$, 
           which is used as the boundary condition in the Poisson equation in KIM.
        '''
        path = os.path.join(self.data_path, self.mode_string, 'fields', 'Phi_aligned.dat')
        try:
            phi_A = np.loadtxt(path)
            self.set_profile(phi_A[:, 0], phi_A[:, 1] + 1j * phi_A[:, 2], 'Phi_aligned')
        except FileNotFoundError:
            print('Phi_aligned file not found in: ' + path)

    def get_B0(self) -> None:
        path = os.path.join(self.data_path, self.mode_string, 'backs', 'B0.dat')
        try:
            B0 = np.loadtxt(path)
            self.set_profile(B0[:, 0], B0[:, 1], 'B0')
        except FileNotFoundError:
            print('B0 file not found in: ' + path)

    def normalize_by_B0(self, dict_entry: str) -> None:
        if self.data['B0'] is not None:
            self.data[dict_entry + '/B0'] = self.data[dict_entry] / self.data['B0']
        else:
            print('Envoke get_B0 first to normalize profile by B0')
        
    def set_r(self, r: NDArray[np.float64]) -> None:
        if self.data['r'] is None:
            self.data['r'] = np.asarray(r)
            self.r_is_set = True

    def set_profile(self, r: NDArray[np.float64], quant: NDArray[np.float64], quant_string: str) -> None:
        if not self.r_is_set:
            self.set_r(r)
        if (np.asarray(r) == self.data['r']).all():
            self.data[quant_string] = quant
        else:
            self.data[quant_string] = np.interp(self.data['r'], r, quant)


