from abc import ABC, abstractmethod

# Abstract product:
class KIMDispersionEquation(ABC):

    @abstractmethod
    def initialize(self, options, species, spec_dat, general_dat):
        pass
    
    @abstractmethod
    def dispersion_equation(self, kr: complex, r_indx: int) -> complex:
        pass