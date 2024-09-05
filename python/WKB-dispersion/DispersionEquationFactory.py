from abc import ABC, abstractmethod

class DispersionEquationFactory(ABC):

    @abstractmethod
    def get_dispersion_equation(self):
        pass
    
    