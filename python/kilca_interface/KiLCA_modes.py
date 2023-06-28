from save_file import save_file
import numpy as np

class KiLCA_modes:
    """"
    Description:
        Class containing the information of the antenna modes in KiLCA.
    Variables:
        m,n (default = 3,2)
        ind, BLUEPRINT
    Methods:
        KiLCA_modes(m,n)
        get_mode_str()
        write(path_to)
    """

    ind = []
    BLUEPRINT = 'modes.in'

    m = 3
    n = 2

    def __init__(self, m=3, n=2):
        """
        Constructor of the modes class. Optional one can set m,n initially
        """
        self.m = m
        self.n = n

    def get_mode_str(self):
        """
        Description:
            Generates from mode variables m and n a list of strings in the form
                (m0, n0)
                (m1, n1)
                ...
            that is used to write in KiLCA input file "modes.in"
        """
        s = []
        if type(self.m) == list or type(self.m) == np.ndarray:
            for k in range(0, len(self.m)):
                s.append('('+ str(int(self.m[k])) + ', ' + str(int(self.n[k])) + ')')
        else:
            s = ['('+ str(int(self.m)) + ', ' + str(int(self.n)) + ')']
        return s



    def write(self, path_to):
        """
        Description:
            Writes modes into KiLCA input file. Doesn't need a blueprint file.
        """
        self.s = self.get_mode_str()
        save_file(self.s, path_to + self.BLUEPRINT)