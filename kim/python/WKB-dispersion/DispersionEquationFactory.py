from abc import ABC, abstractmethod
from KIMDispersion_Horton import KIMDispersion_Horton
from KIMDispersion_Krook import KIMDispersion_Krook
from KIMDispersion_FokkerPlanck import KIMDispersion_FokkerPlanck
from KIMDispersionEquation import KIMDispersionEquation
    

    
class DispersionEquationFactory(ABC):

    @abstractmethod
    def get_dispersion_model(self, options, species, spec_dat, general_dat) -> complex:
        pass

class KIMDispersion_Horton_Factory(DispersionEquationFactory):
    
    def get_dispersion_model() -> KIMDispersionEquation:
        return KIMDispersion_Horton()

class KIMDispersion_Krook_Factory(DispersionEquationFactory):
    
    def get_dispersion_model() -> KIMDispersionEquation:
        return KIMDispersion_Krook()

class KIMDispersion_FokkerPlanck_Factory(DispersionEquationFactory):
    
    def get_dispersion_model() -> KIMDispersionEquation:
        return KIMDispersion_FokkerPlanck()