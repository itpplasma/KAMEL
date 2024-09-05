from KIMDispersionEquation import KIMDispersionEquation
from plasmapy.dispersion import plasma_dispersion_func as plasma_disp
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
from constants import *

class KIMDispersion_Horton(KIMDispersionEquation):
    
    def initialize(self, options: dict, species: dict, spec_dat: dict, general_dat: dict, equil_dat: dict):
        """
            Input are dictionary objects from the KIM_WKB class.
        """
        self.options = options
        self.species = species
        self.spec_dat = spec_dat
        self.general_dat = general_dat
        self.equil_dat = equil_dat
    
    def dispersion_equation(self, kr: complex, r_indx: int) -> complex:
        self.general_dat['kperp'] = np.sqrt(self.general_dat['ks']**2 + kr**2)
        ky = self.general_dat['ks'][r_indx]
        om_prime = self.options['omega'] - self.general_dat['om_E'][r_indx]
        
        dispersion_equation = 0 + 0j
        denom = 0 + 0j
        nom = 0 + 0j
        for spec in self.species:
            self.spec_dat[spec]['om_n'] = ky * sol / (self.spec_dat[spec]['charge'] * self.equil_dat['B0'][r_indx] \
                * self.spec_dat[spec]['n'][r_indx]) * self.spec_dat[spec]['T'][r_indx] * ev * self.spec_dat[spec]['dndr'][r_indx]
            self.spec_dat[spec]['om_T'] = ky * sol / (self.spec_dat[spec]['charge'] * self.equil_dat['B0'][r_indx]) * self.spec_dat[spec]['dTdr'][r_indx] * ev

            self.spec_dat[spec]['z'] = om_prime / (np.sqrt(2) * self.general_dat['kp'][r_indx] * self.spec_dat[spec]['vT'][r_indx])
            self.spec_dat[spec]['W'] = (self.spec_dat[spec]['om_n'] - om_prime + self.spec_dat[spec]['om_T'] * (self.spec_dat[spec]['z']**2 + 0.5)) \
                * self.spec_dat[spec]['z'] / om_prime * plasma_disp(self.spec_dat[spec]['z']) + self.spec_dat[spec]['om_T'] * self.spec_dat[spec]['z']**2 / om_prime

            self.spec_dat[spec]['nom'] = self.spec_dat[spec]['lambda_D'][r_indx]**-2 * (self.spec_dat[spec]['W'] - 1.0 \
                -self.spec_dat[spec]['om_T'] * self.spec_dat[spec]['z'] / om_prime * plasma_disp(self.spec_dat[spec]['z']))
            self.spec_dat[spec]['denom'] = self.spec_dat[spec]['W'] * self.spec_dat[spec]['vT'][r_indx]**2 / (self.spec_dat[spec]['omega_c'][r_indx]**2 \
                * self.spec_dat[spec]['lambda_D'][r_indx]**2)

            nom += self.spec_dat[spec]['nom']
            denom += self.spec_dat[spec]['denom']

        nom -= self.general_dat['kp'][r_indx]**2
        dispersion_equation = nom / (1 + denom) - self.general_dat['kperp'][r_indx]**2
        
        return dispersion_equation