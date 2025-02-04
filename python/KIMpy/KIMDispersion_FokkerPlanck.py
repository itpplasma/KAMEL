from .KIMDispersionEquation import KIMDispersionEquation
from .Bessel_calculation import calc_needed_bessel_of_mphi
import numpy as np

class KIMDispersion_FokkerPlanck(KIMDispersionEquation):
    
    def initialize(self, options:dict, species:dict, spec_dat: dict, general_dat: dict, equil_dat: dict):
        self.options = options
        self.species = species
        self.spec_dat = spec_dat
        self.general_dat = general_dat
        self.equil_dat = equil_dat
    
    def dispersion_equation(self, kr, r_indx):
        self.general_dat['kperp'] = np.sqrt(self.general_dat['ks']**2 + kr**2)
        dispersion_equation = self.general_dat['kperp'][r_indx]**2 + self.general_dat['kp'][r_indx]**2 + 0j
        #print(dispersion_equation)
        for spec in self.species:
            eval_b = self.general_dat['kperp'][r_indx]**2 * self.spec_dat[spec]['rho_TL'][r_indx]**2
            for m_phi in range(0, self.options['max_cyclotron_harmonic']+1):
                BesselProd0, BesselProd1 = calc_needed_bessel_of_mphi(mphi=m_phi, eval_b=eval_b)
                if np.isnan(BesselProd0) or np.isnan(BesselProd1):
                    print(f'FokkerPlanck: BesselProd0 or BesselProd1 contains NaNs, eval_b = {eval_b}, m_phi = {m_phi}')
                    print(f'BesselProd0: {BesselProd0}, BesselProd1: {BesselProd1}')

                dispersion_equation_temp = \
                   (1j * self.spec_dat[spec]['vT'][r_indx]**2.0 * self.general_dat['ks'][r_indx]) \
                       / (self.spec_dat[spec]['lambda_D'][r_indx]**2 * self.spec_dat[spec]['omega_c'][r_indx] * self.spec_dat[spec]['nu'][r_indx]) \
                    * (
                        self.spec_dat[spec]['A1'][r_indx] * self.spec_dat[spec]['I00'][m_phi][r_indx] * BesselProd0 \
                        + self.spec_dat[spec]['A2'][r_indx] \
                        * (
                           self.spec_dat[spec]['I00'][m_phi][r_indx] * (
                                 BesselProd0 * (1 + eval_b + m_phi) + eval_b * BesselProd1
                            ) 
                           + 0.5 * self.spec_dat[spec]['I20'][m_phi][r_indx] * BesselProd0
                        )
                    )
                dispersion_equation -= dispersion_equation_temp
        return dispersion_equation