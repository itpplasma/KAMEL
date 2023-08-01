import numpy as np


class KiLCA_coil():
    """Class that handles the coil current value. Note that KiLCA/QL-Balance uses the
    actual current value in statA, while GPEC uses the value per turn."""

    available_machines = ['AUG', 'MASTU']
    windings_machines = [5, 4] # only for information, not used in KiLCA/QL-Balance
    available_types = ['default', 'file']
    I0 = None

    def __init__(self, type='default', machine='AUG'):
        """Constructor of the KiLCA coil class.
            Input:
                type ... either "default" for the default value of the machine, 
                or "file" for reading it from a file, e.g. usual ASDEX coil file"""
        self.machine = machine
        if not self.machine in self.available_machines:
            raise ValueError(f'Coil information for machine {self.machine} not available')
        
        self.type = type
        if not self.type in self.available_types:
            raise ValueError(f'Type {self.type} not available.')

        if self.type == self.available_types[0]:
            self.set_default_value()

    def read_from_coil_file(self, path_to_file):
        """Read from a coil file"""
        if self.machine == 'AUG':
            pass

    def set_default_value(self):
        if self.machine == 'AUG':
            self.I0 = 5.1e12
        if self.machine == 'MASTU':
            self.I0 = 6.3e12

