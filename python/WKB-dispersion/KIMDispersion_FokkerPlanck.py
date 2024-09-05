from DispersionEquationFactory import DispersionEquationFactory

class KIMDispersion_FokkerPlanck(DispersionEquationFactory, KIM_Krook):
    
    def __init__(self, spec_dat: dict, general_dat: dict):
        self.spec_dat = spec_dat
        self.general_dat = general_dat
    
    def get_dispersion_equation(self):
        pass