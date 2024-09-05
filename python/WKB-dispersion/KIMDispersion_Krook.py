from KIMDispersionEquation import KIMDispersionEquation

class KIMDispersion_Krook(KIMDispersionEquation):
    """
        Class to calculate the dispersion equation for the collisionless model and Krook collisions. 
        The former is determined by simply setting the collision frequency to zero.
    """
    
    def init_dispersion_model(self, options: dict, species: dict, spec_dat: dict, general_dat: dict, equil_dat: dict):
        self.options = options
        self.species = species
        self.spec_dat = spec_dat
        self.general_dat = general_dat
        self.equil_dat = equil_dat
    
    def dispersion_equation(self, kr: complex, r_indx: int):
        dispersion_equation = self.general_dat['kperp'][r_indx]**2 + self.general_dat['kp'][r_indx]**2
        for spec in self.species:
            eval_b = self.general_dat['kperp'][r_indx]**2 * self.spec_dat[spec]['rho_TL'][r_indx]**2
            BesselProd0, BesselProd1 = calc_needed_bessel_of_mphi(0, eval_b)
            if np.isnan(BesselProd0) or np.isnan(BesselProd1):
                print(f'BesselProd0 or BesselProd1 contains NaNs')

            dispersion_equation += self.spec_dat[spec]['lambda_D'][r_indx]**-2 * \
                (
                    1.0 - self.general_dat['ks'][r_indx] * self.spec_dat[spec]['rho_TL'][r_indx] \
                / (self.general_dat['kp'][r_indx] * np.sqrt(2)) * (self.spec_dat[spec]['A1'][r_indx] * BesselProd0 * plasma_disp(self.spec_dat[spec]['z0'][r_indx]) + \
                self.spec_dat[spec]['A2'][r_indx] * (plasma_disp(self.spec_dat[spec]['z0'][r_indx]) * (1 + eval_b + self.spec_dat[spec]['z0'][r_indx]**2)\
                * np.exp(-eval_b) + BesselProd1 * eval_b + self.spec_dat[spec]['z0'][r_indx] * BesselProd0))
                )
        #print(dispersion_equation)
        return dispersion_equation