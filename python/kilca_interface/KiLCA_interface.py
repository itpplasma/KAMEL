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
import warnings
import time

from KiLCA_antenna import KiLCA_antenna
from KiLCA_background import KiLCA_background
from KiLCA_eigmode import KiLCA_eigmode
from KiLCA_modes import KiLCA_modes
from KiLCA_output import KiLCA_output
from KiLCA_zone import KiLCA_zone

class KiLCA_interface:
    """
    Description:
        KiLCA interface class to manage input/output files of KiLCA and run the code.
    Variables:
        EXEC_PATH              ... path of the KiLCA executable
        BLUE_PATH              ... path containing the blueprints for the input files
        PROF_PATH              ... path to the input profiles
        antenna                ... object of the KiLCA_antenna class
        background             ... object of the KiLCA_background class
        eigmode                ... object of the KiLCA_eigmode class
        modes                  ... object of the KiLCA_modes class
        output                 ... object of the KiLCA_output class
        zones                  ... list of KiLCA_zone objects
        path                   ... path of the directory containing the run path
        path_of_interface_file ... path of the KiLCA_interface class file
        path_of_profiles       ... path where the profiles will be copied to
        path_of_run            ... run path. This is path + '/' + run type
        run_type               ... type of the run (flre, vacuum)
        machine                ... which machine, currently only AUG (ASDEX Upgrade)
    Methods:
        KiLCA_interface(shot: int, time: int, path: str, rtype: str, machine='AUG')
            Constructor of the class.
            Arguments:
                shot ... int of the shot number
                time ... int of the time slice
                path ... path to the directory containing the run path
                rtype ... run type (flre or vacuum)
                machine ... for machine settings, currently only AUG (ASDEX Upgrade)
        set_modes(m, n)
            Set the RMP mode numbers.
            Arguments:
                m ... poloidal mode number, either int or list/numpy array of ints
                n ... toroidal mode number, either int or list/numpy array of ints,
                      must be same length as m
        set_antenna(ra: float, nmod: int)
            Initializes KiLCA_antenna object
            Arguments:
                ra   ... radius of the RMP antenna
                nmod ... number of modes
        set_background(Rtor: float, rpl: float)
            Initializes KiLCA_background object
            Arguments:
                Rtor ... Major radius of the device (measured at the magnetic axis)
                rpl  ... Minor plasma radius of the device
        set_zones(r: list, b: list, m: list)
            Initializes list of zone objects.
            Arguments:
                r ... list of positions of the zone boundaries
                b ... list of type of zone boundaries (center, infinity, interface, antenna)
                m ... list of media between zone boundaries
        set_ASDEX(nmodes: int)
            Sets AUG machine properties (background and antenna)
            Arguments:
                nmodes ... number of RMP modes
        set_MASTU() !!! not implemented yet !!!
        write()
            Creates directory structure, copies all needed files

        run()
            Run KiLCA in path of run

        run_remote() !!! not implemented yet !!!
        run_condor() !!! not implemented yet !!!
    ################################################################################
    Minimum requirement for a run:
        obj = KiLCA_interface(shot, time, path, rtype)
        obj.PROF_PATH = 'path_to_existing_profiles'
        obj.set_modes(m,n)
        obj.write()
        obj.run()
    ################################################################################
    """

    EXEC_PATH = '/proj/plasma/soft/KiLCA-2.4.2/exe/KiLCA_Normal_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit'
    BLUE_PATH = 'blueprints/'
    PROF_PATH = 'profiles/'

    antenna = KiLCA_antenna()
    background = KiLCA_background()
    eigmode = KiLCA_eigmode()
    modes = KiLCA_modes
    output = KiLCA_output()

    zones = []

    path = ''
    path_of_interface_file = ''
    path_of_profiles = ''
    path_of_run = ''

    run_type = ''
    machine = ''


    def __init__(self, shot: int, time: int, path: str, rtype: str, machine: str='AUG'):
        """Constructor of KiLCA interface.
        input:
                shot ... shot number
                time ... time of the time slice in experiment
                path ... path to run directory, must end with /. Should contain profile directory '/profiles'. Otherwise, change PROF_PATH
                rtype... run type, i.e. vacuum or flre
                machine ... e.g. AUG or LHD
                """

        self.path = path
        try:
            ls = os.listdir(self.path)
        except:
            raise ValueError('Path ' + self.path + ' does not exist.')

        if not self.PROF_PATH[0:-1] in ls:
            warnings.warn('No profile directory found in ' + self.path + '\nMake sure to change PROF_PATH of class to path where the profiles are.')
            
        self.path_of_profiles = path + self.PROF_PATH
        self.path_of_interface_file = os.getcwd()

        if not type(shot) == int:
            raise ValueError('shot argument is not of type integer')
        self.shot = shot # shot number of the experiment

        if not type(time) == int:
            raise ValueError('time argument is not of type integer')
        self.time = time # time of the time slice


        if rtype == 'vacuum' or rtype == 'flre' or rtype == 'imhd':
            self.path_of_run = self.path + rtype + '/'
            self.run_type = rtype
        else:
            raise ValueError('Runtype not supported')

        self.machine = machine # machine, e.g. AUG (not implemented yet: MASTU)
        if self.machine == 'AUG':
            self.set_ASDEX()
        elif self.machine == 'MASTU':
            self.set_MASTU()
        else:
            raise ValueError('Machine not supported')

        self.BLUE_PATH = inspect.getfile(KiLCA_interface)[0:-18] + self.BLUE_PATH

        
        
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

    def set_antenna(self, ra: float = 67.0, nmod: int = 1):
        """Initializes the antenna with the minimum amount of needed 
        parameters to run.
        input:
                ra  ... position of antenna in cm (small radius)
                nmod... number of modes to calculate"""
        
        if ra < 0:
            raise ValueError('Radius of antenna must be above 0')
        if nmod < 0:
            raise ValueError('Number of modes invalid')

        self.antenna = KiLCA_antenna(ra,nmod)

    def set_background(self, Rtor: float = 170.0, rpl: float = 63.0):
        """Initializes the background with the minimum amount of needed
        parameters to run.
        input:
                Rtor ... big torus radius in cm
                rpl  ... plasma radius in cm (small radius)"""
        if Rtor < 0:
            raise ValueError('Rtor must be above 0')
        if rpl < 0:
            raise ValueError('rpl must be above 0')

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
        """
        Description:
            Run KiLCA in path of run. Needs write first.
        """
        
        if self.ready_to_run == False:
            raise ValueError('Not ready to run!')

        cdir = os.getcwd()
        os.chdir(self.path_of_run)
        start = time.time()
        self.output = sp.Popen("./run_local", stdout=sp.PIPE)
        self.output.wait()
        end = time.time()
        os.chdir(cdir)
        print(f'The KiLCA run took {end-start}s')

        print(self.output.stdout.read())
        

    def run_remote(self):
        pass

    def run_condor(self):
        pass


