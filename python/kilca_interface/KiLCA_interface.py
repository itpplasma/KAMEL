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
import inspect
import subprocess as sp

from KiLCA_antenna import KiLCA_antenna
from KiLCA_background import KiLCA_background
from KiLCA_eigmode import KiLCA_eigmode
from KiLCA_modes import KiLCA_modes
from KiLCA_output import KiLCA_output
from KiLCA_zone import KiLCA_zone

class KiLCA_interface:

    EXEC_PATH = '/proj/plasma/soft/KiLCA-2.4.2/exe/KiLCA_Normal_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit'
    BLUE_PATH = 'blueprints/'
    PROF_PATH = 'profiles/'

    antenna = KiLCA_antenna()
    background = KiLCA_background()
    eigmode = KiLCA_eigmode()
    modes = KiLCA_modes
    output = KiLCA_output()

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
        self.machine = machine # machine, e.g. AUG or MASTU

        self.BLUE_PATH = inspect.getfile(KiLCA_interface)[0:-18] + self.BLUE_PATH

        if rtype == 'vacuum' or rtype == 'flre' or rtype == 'imhd':
            self.path_of_run = self.path + rtype + '/'
            self.run_type = rtype
        else:
            raise ValueError("Runtype not supported")

    def set_modes(self, m=3, n=2):
        """
        Description:    
            Set mode numbers m and n. Must either be scalar or lists/numpy arrays.
        """
        self.m = m
        self.n = n
        try:
            self.n_modes = len(m)
        except:
            print('single RMP mode')
        self.modes = KiLCA_modes(m, n)

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
        if not (type(r)==list):
            raise ValueError("r must be a vector")
        # check if size is suitable
        if not (len(r) == len(b)):
            raise ValueError("Size of r and b does not match")
        if not (len(r) == len(m)+1):
            raise ValueError("Size of r and m does not match")

        self.zones = []
        for k in range(0,len(m)):
            self.zones.append(KiLCA_zone(k, r[k], b[k], m[k], r[k+1], b[k+1]))

    def set_ASDEX(self, nmodes=1):
        """
        Description:
            Initializes the class for a standard run on ASDEX parameters.
        Input:
            nmodes ... number of modes to calculate, default=0
        """

        self.set_antenna(70, nmodes)
        self.set_background(170.05, 67.0)
        self.background.Btor = -17563.3704

        # set zones
        r = [3.0, 67.0, 70.0, 80.0]
        b = ['center', 'interface', 'antenna', 'idealwall']
        m = [self.run_type, 'vacuum', 'vacuum']
        self.set_zones(r,b,m)
    
    def set_MASTU(self, nmodes=0):
        """
        Description:
            Initializes the class for a standard run on MASTU parameters
        Input:
            nmodes ... number of modes to calculate, default = 0
        """

        raise ValueError('Not implemented yet')

        pass


    def write(self):
        """
        Description:
            Creates directory structure, copies all needed files. Needs to be called before .run()
        """

        # create paths if they are not there
        os.system('mkdir -p ' + self.path)
        os.system('mkdir -p ' + self.path_of_profiles)
        os.system('mkdir -p ' + self.path_of_run)

        # make local copy of profiles
        os.system('cp ' + self.PROF_PATH + '* ' + self.path_of_profiles + ' 2>/dev/null') # suppress warnings
        # delete all existing input files
        os.system('rm -f ' + self.path + 'background.in')
        os.system('rm -f ' + self.path + 'eigmode.in')
        os.system('rm -f ' + self.path + 'modes.in')
        os.system('rm -f ' + self.path + 'output.in')
        os.system('rm -f ' + self.path_of_run + 'antenna.in')
        os.system('rm -f ' + self.path_of_run + 'zone*.in')

        # create sym link to exe
        os.system('ln -sf ' + self.EXEC_PATH + ' ' + self.path_of_run + 'run_local')
        # create symbolic link to profiles directory
        os.system('ln -sfT ' + self.path_of_profiles + ' ' + self.path_of_run + 'profiles')

        self.antenna.write(self.BLUE_PATH + self.antenna.BLUEPRINT, self.path_of_run)

        self.background.write(self.BLUE_PATH + self.background.BLUEPRINT, self.path)
        os.system('ln -sf ' + self.path + self.background.BLUEPRINT + ' ' + self.path_of_run + self.background.BLUEPRINT)

        self.eigmode.write(self.BLUE_PATH + self.eigmode.BLUEPRINT, self.path)
        os.system('ln -sf ' + self.path + self.eigmode.BLUEPRINT + ' ' + self.path_of_run + self.eigmode.BLUEPRINT)

        self.output.write(self.BLUE_PATH + self.output.BLUEPRINT, self.path)
        os.system('ln -sf ' + self.path + self.output.BLUEPRINT + ' ' + self.path_of_run + self.output.BLUEPRINT)

        self.modes.write(self.path)
        os.system('ln -sf ' + self.path + self.modes.BLUEPRINT + ' ' + self.path_of_run + self.modes.BLUEPRINT)

        for k in range(0, len(self.zones)):
            self.zones[k].write(self.BLUE_PATH, self.path_of_run)

        self.ready_to_run = True

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


