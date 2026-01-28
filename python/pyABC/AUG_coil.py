import f90nml
import numpy as np


class AUG_coil:
    """class AUG_coil
    Description:
    Class to manage coil files from AUG experiments. In these files the vector [Iu, Il] is written in a row (upper and lower coil currents).
    """

    def __init__(self, file):
        """Description:
        Constructor of the AUG_coil class
        Input:
        - file ... path to coil file (experiment)"""

        self.path = file

    def read(self, flip=False):
        """Description:
        Reads the coil file. If needed, the sequence of Iu and Il can be flipped: instead of [Iu, Il] -> [Il,Iu] is read in.
        """

        raw = np.loadtxt(self.path)
        len_raw = len(raw)
        print(len_raw)

        if flip == False:
            self.Iu = raw[0 : int(len_raw / 2)]
            self.Il = raw[int(len_raw / 2) : len_raw]
        else:
            self.Il = raw[0 : int(len_raw / 2)]
            self.Iu = raw[int(len_raw / 2) : len_raw]
