# TODO:
# - Add database part, i.e. write a database entry for everytime KiLCA is run. If the kilca interface
# is used by the balance interface, this should be slimed. Database entry should also contain git hash.
# - add remote connection, such that when cloned on different machine and a itp account is available
# kilca can be run remotely.
#


import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import subprocess as sp

from KiLCA_antenna import KiLCA_antenna
from KiLCA_background import KiLCA_background
from KiLCA_eigmode import KiLCA_eigmode
from KiLCA_modes import KiLCA_modes
from KiLCA_output import KiLCA_output

class KiLCA_interface:

    EXEC_PATH = ''
    BLUE_PATH = 'blueprints/'
    PROF_PATH = 'profiles/'

    antenna = KiLCA_antenna
    background = KiLCA_background
    eigmode = KiLCA_eigmode
    modes = KiLCA_modes
    output = KiLCA_output

    zones = {}

    path = ''
    path_of_interface_file = ''
    path_of_profiles = ''
    path_of_run = ''

    run_type = ''


    def __init__(self, shot, time, path, rtype, machine='AUG'):
        """Constructor of KiLCA interface.
        input:
                shot ... shot number
                time ... time of the time slice in experiment
                path ... path to run directory, must end with /
                rtype... run type, i.e. vacuum or flre
                machine ... e.g. AUG or LHD
                """

        self.path = path
        self.path_of_profiles = path + self.PROF_PATH
        self.path_of_interface_file = os.getcwd()
        self.shot = shot # shot number of the experiment
        self.time = time # time of the time slice
        self.machine = machine # machine, e.g. AUG or LHD

        if rtype == 'vacuum' or rtype == 'flre' or rtype == 'imhd':
            self.path_of_run = self.path + rtype + '/'
            self.run_type = rtype
        else:
            raise ValueError("Runtype not supported")


    def set_antenna(self, ra=67, nmod=1):
        """Initializes the antenna with the minimum amount of needed 
        parameters to run.
        input:
                ra  ... position of antenna in cm (small radius)
                nmod... number of modes to calculate"""
        self.antenna = KiLCA_antenna(ra,nmod)

    def set_background(self, Rtor=170, rpl=63):
        """Initializes the background with the minimum amount of needed
        parameters to run.
        input:
                Rtor ... big torus radius in cm
                rpl  ... plasma radius in cm (small radius)"""

        self.background = KiLCA_background(Rtor, rpl)

    def set_zones(self, r,b,m):
        """Initializes the zones of the run.
        input:
                r ... position of zone boundaries
                b ... type of zone boundaries (equally long to r)
                    (center, infinity, interface, antenna)
                m ... media between zone boundaries (1 element less than r)
                    (vacuum, medium, imhd, rmhd, flre)"""

        # check if input is array
        if not (type(r)==np.ndarray):
            raise ValueError("r must be a vector")
        # check if size is suitable
        if not (len(r) == len(b)):
            raise ValueError("Size of r and b does not match")
        if not (len(r) == len(m)+1):
            raise ValueError("Size of r and m does not match")

        self.zones = []
        for k in range(0,len(m)):
            self.zones[k] = KiLCA_zone(k, r[k], b[k], m[k], r[k+1], b[k+1])

    def write(self):
        pass

    def run(self):

        cdir = os.getcwd()
        os.chdir(self.path_of_run)
        self.output = sp.Popen("./run_local", stdout=sp.PIPE)
        os.chdir(cdir)

        print(self.output.stdout.read())
        

    def run_remote(self):
        pass

    def run_condor(self):
        pass


