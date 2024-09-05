from DispersionEquationFactory import DispersionEquationFactory

class KIMDispersion_Krook(DispersionEquationFactory):
    """
        Class to calculate the dispersion equation for the collisionless model and Krook collisions. 
        The former is determined by simply setting the collision frequency to zero.
    """
    
    def __init__(self, spec_dat: dict, general_dat: dict):
        self.spec_dat = spec_dat
        self.general_dat = general_dat
    
    def get_dispersion_equation(self):
        pass