import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import inspect
import subprocess as sp
import warnings
import time
import copy
import inspect
import sys

from .KiLCA_antenna import *#KiLCA_antenna
from .KiLCA_background import KiLCA_background
from .KiLCA_eigmode import KiLCA_eigmode
from .KiLCA_modes import KiLCA_modes
from .KiLCA_output import KiLCA_output
from .KiLCA_zone import KiLCA_zone

from device_config import MASTU_config, AUG_config

#sys.path.append(os.path.abspath(inspect.getfile(KiLCA_antenna)[0:-16] + '../../postproc_py_class/'))
from postproc_class import utility_class


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
        set_MASTU(nmodes: int) 
            Sets MASTU machine properties (major radius, plasma radius, toroidal B field)
            Arguments:
                nmodes ... number of RMP modes
        write()
            Creates directory structure, copies all needed files

        run()
            Run KiLCA in path of run

        run_remote() !!! not implemented yet !!!
        run_condor() !!! not implemented yet !!!

        create_parabolic_profiles_from_res_surf(path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, m_mode, n_mode, rmin, rmax, num, a, const='')
            Creates parabolic profiles starting from the value of the pressure at the resonant surface.

        check_profile_consistency(path_profiles)
            Checks the consistency of the profiles, i.e. firstly, if they even exist and secondly,
            if density and temperatures are greater than zero.
            
        plt_profiles(path_profiles)
            Plot profiles (ne, Te, Ti, Vz, Vth, Er, q) for fast check.

    ################################################################################
    Minimum requirement for a run:
        obj = KiLCA_interface(shot, time, path, rtype)
        obj.PROF_PATH = 'path_to_existing_profiles'
        obj.set_modes(m,n)
        obj.write()
        obj.run()
    ################################################################################
    """

    EXEC_PATH = os.path.join(os.path.dirname(__file__) + '/../../KiLCA/build/exe/KiLCA_Normal_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit')
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

    uc = utility_class.utility() # utility class for colors and adding grid lines to plots

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
            warnings.warn('Path ' + self.path + ' does not exist. I will make dir(s).')
            os.makedirs(self.path)

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
        #self.set_machine()

        self.BLUE_PATH = inspect.getfile(KiLCA_interface)[0:-18] + self.BLUE_PATH

    def set_machine(self, delta_r_antenna: float = 3.0):
        if self.machine == 'AUG':
            self.set_ASDEX()
        elif self.machine == 'MASTU':
            self.set_MASTU(delta_r_antenna=delta_r_antenna)
        elif self.machine == None:
            # don't set the machine
            print('Will not set a machine at constructor')
        else:
            raise ValueError('Machine not supported')
        
        
    def set_modes(self, m: int = 3, n: int = 2):
        """
        Description:    
            Set mode numbers m and n. Must either be scalar or lists/numpy arrays.
        """
        self.m = m
        self.n = n
        try:
            self.n_modes = len(m)
        except:
            #print(f'Single RMP mode: m = {int(m)}, n = {int(n)}')
            pass
        self.modes = KiLCA_modes(m, n)

    def set_antenna(self, ra: float = 67.0, nmod: int = 1, I0: float = 4.5e12):
        """Initializes the antenna with the minimum amount of needed 
        parameters to run.
        input:
                ra  ... position of antenna in cm (small radius)
                nmod... number of modes to calculate"""
        
        if ra < 0:
            raise ValueError('Radius of antenna must be above 0')
        if nmod < 0:
            raise ValueError('Number of modes invalid')

        self.antenna = KiLCA_antenna(ra,nmod, I0)

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

        self.zones = [copy.deepcopy(KiLCA_zone(k, r[k], b[k], m[k], r[k+1], b[k+1])) for k in range(0, len(m))]

    def set_ASDEX(self, nmodes: int = 1, ra: float = 70, I0: float = 5.1e12):
        """
        Description:
            Initializes the class for a standard run on ASDEX parameters.
        Input:
            nmodes ... number of modes to calculate, default=0
            ra ... radius of the antenna
            I0 ... RMP coil current value in statA, AUG: 1.7kA = 5.1e12 statA
        """
        self.machine = 'AUG'
        #print(f'Machine setting: AUG, type: {self.run_type}')

        if ra < 67:
            raise ValueError('Radius of antenna inside plasma')
        if ra > 80:
            raise ValueError('Radius of antenna outside wall')

        self.a_minor = 67.0
        self.R0 = 165.0

        self.set_antenna(ra, nmodes, I0 = I0)
        self.set_background(Rtor=self.R0, rpl=self.a_minor)
        self.background.data['Btor'] = -17563.3704

        # set zones
        r = [3.0, 67.0, 70.0, 80.0]
        b = ['center', 'interface', 'antenna', 'idealwall']
        m = [self.run_type, 'vacuum', 'vacuum']
        self.set_zones(r,b,m)
    
    def set_MASTU(self, nmodes: int = 1, I0: float = 6.0e12, delta_r_antenna: float = 3.0):
        """
        Description:
            Initializes the class for a standard run on MASTU parameters
        Input:
            nmodes ... number of modes to calculate, default = 0
            I0 ... RMP coil current value in statA, MASTU: 2.0kA=6.0e12statA (MAST-U has 4 turns 
            per RMP coil with a maximum current of 8kAt)
            delta_r_antenna ... effective radius distance of antenna to plasma
        """
        self.machine = 'MASTU'
        #print(f'Machine setting: MASTU, type: {self.run_type}')
        self.machine_config = MASTU_config()
        self.R0 = self.machine_config.R0
        self.a_minor = self.machine_config.r_eff_plasma

        self.r_antenna = self.machine_config.r_eff_plasma + delta_r_antenna
        print('r antenna: ', self.r_antenna)

        self.set_antenna(self.r_antenna, nmodes, I0 = I0)
        self.set_background(self.machine_config.R0, self.machine_config.r_eff_plasma)
        #self.background.data['Btor'] = -6400.0

        # set zones
        r = [3.0, self.machine_config.r_eff_plasma, self.r_antenna, self.machine_config.r_eff_wall]
        b = ['center', 'interface', 'antenna', 'idealwall']
        m = [self.run_type, 'vacuum', 'vacuum']
        self.set_zones(r,b,m)



    def write(self):
        """
        Description:
            Creates directory structure, copies all needed files. Needs to be called before .run()
        """

        self.check_profile_consistency()

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
        os.system('ln -sf ' + os.path.abspath(self.path_of_profiles) + ' ' + self.path_of_run + 'profiles')

        self.antenna.write(self.BLUE_PATH + self.antenna.BLUEPRINT, self.path_of_run)

        self.background.write(self.BLUE_PATH + self.background.BLUEPRINT, self.path)
        os.system('ln -sf ' + os.path.abspath(self.path) + '/' + self.background.BLUEPRINT + ' ' + self.path_of_run + self.background.BLUEPRINT)

        self.eigmode.write(self.BLUE_PATH + self.eigmode.BLUEPRINT, self.path)
        os.system('ln -sf ' + os.path.abspath(self.path) + '/' + self.eigmode.BLUEPRINT + ' ' + self.path_of_run + self.eigmode.BLUEPRINT)

        self.output.write(self.BLUE_PATH + self.output.BLUEPRINT, self.path)
        os.system('ln -sf ' + os.path.abspath(self.path) + '/' + self.output.BLUEPRINT + ' ' + self.path_of_run + self.output.BLUEPRINT)

        self.modes.write(self.path)
        os.system('ln -sf ' + os.path.abspath(self.path) + '/' + self.modes.BLUEPRINT + ' ' + self.path_of_run + self.modes.BLUEPRINT)

        for k in range(0, len(self.zones)):
            self.zones[k].write(self.BLUE_PATH, self.path_of_run)

        self.ready_to_run = True
    
    def get_r_res(self, m_mode, n_mode):
        q_prof = np.loadtxt(self.path_of_profiles + 'q.dat')
        self.r_res = np.interp(m_mode/n_mode, np.abs(q_prof[:,1]), q_prof[:,0])
        print(self.r_res)
        return self.r_res
        

    def run(self):
        """
        Description:
            Run KiLCA in path of run. Needs write first.
        """
        
        if self.ready_to_run == False:
            raise ValueError('Not ready to run!')

        print(f'Machine setting: {self.machine}, type: {self.run_type}')
        try:
            self.n_modes = len(self.m)
        except:
            print(f'Single RMP mode: m = {int(self.m)}, n = {int(self.n)}')

        cdir = os.getcwd()
        os.chdir(self.path_of_run)
        start = time.time()
        self.run_out = sp.Popen("./run_local", stdout=sp.PIPE, stderr=sp.PIPE)
        self.run_out.wait()
        end = time.time()
        os.chdir(cdir)
        print(f'The KiLCA run took {end-start}s')
        
        if self.run_out.stderr.read() == b'':
            if not (self.run_out.stdout.read() == b''):
                print('Output:')
                for line in self.run_out.stdout:
                    print(line.decode('utf8'))
            print('')
        else:
            print('KiLCA ran into problems!')
            print('Errors:')
            for line in self.run_out.stderr:
                print(line.decode('utf8'))
            print('')
        #print(self.output.stdout.read())
        

    def run_remote(self):
        pass

    def run_condor(self):
        pass


    def create_parabolic_profiles_from_res_surf(self, path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, m_mode, n_mode, rmin, rmax, num, a, const=''):
        """ Create parabolic profiles for fixed density and electron
        temperature values at the rational surface. """
        r = np.linspace(rmin, rmax, num)
        q = -(1.05 + q0 * (r/a)**2)
        rres = np.interp(m_mode/n_mode, np.abs(q), r)
    
        fac_par = 1 - (r/a)**2
        n   = n0   * fac_par / (1-(rres/a)**2)
        Te  = Te0  * fac_par / (1-(rres/a)**2)
        Ti  = Ti0  * fac_par
        Vz  = Vz0  * fac_par
        Er  = Er0  * fac_par
        Vth = Vth0 * fac_par

        np.savetxt(path + 'q.dat', np.array((r,q)).transpose())
        np.savetxt(path + 'Te.dat', np.array((r, Te)).transpose())
        np.savetxt(path + 'Ti.dat', np.array((r, Ti)).transpose())
        np.savetxt(path + 'n.dat', np.array((r, n)).transpose())
        np.savetxt(path + 'Vz.dat', np.array((r, Vz)).transpose())
        np.savetxt(path + 'Er.dat', np.array((r, Er)).transpose())
        np.savetxt(path + 'Vth.dat', np.array((r, Vth)).transpose())


    def check_profile_consistency(self, path_profiles=''):
        """Check consistency of profiles, i.e. if temperatures and density
        are above zero and if all profile files exist."""

        if path_profiles=='':
            path_profiles = self.path_of_profiles

        req_profs = ['n.dat', 'Te.dat', 'Ti.dat', 'Vz.dat', 'q.dat', 'Er.dat', 'Vth.dat']
        ls = os.listdir(path_profiles)

        count = 0
        for prof_type in req_profs:
            if not prof_type in ls:
                raise ValueError(f'Profile file {prof_type} is missing in {path_profiles}')
            if count < 3:
                dat = np.loadtxt(path_profiles + prof_type)
                if np.any(dat[:,1] <0):
                    raise ValueError(f'Profile file {prof_type} contains negative values, but it should not!')
            count = count +1


    def plt_profiles(self, path_profiles=''):
        """Plot all profiles used in KiLCA."""
        
        if path_profiles=='':
            path_profiles = self.path_of_profiles

        self.check_profile_consistency()
        req_profs = ['n.dat', 'Te.dat', 'Ti.dat', 'Vz.dat', 'q.dat', 'Er.dat', 'Vth.dat']

        prof_data = {}
        for prof in req_profs:
            prof_data[prof] = np.loadtxt(path_profiles + prof)
        
        fig,ax = plt.subplots(4,2, figsize=(8,6), sharex=True)
        ax[0,0].plot(prof_data['n.dat'][:,0], prof_data['n.dat'][:,1])
        ax[0,0].set_ylabel(r'n [cm$^{-3}$]')

        ax[1,0].plot(prof_data['Te.dat'][:,0], prof_data['Te.dat'][:,1])
        ax[1,0].set_ylabel(r'T$_e$ [eV]')           

        ax[2,0].plot(prof_data['Ti.dat'][:,0], prof_data['Ti.dat'][:,1])
        ax[2,0].set_ylabel(r'T$_i$ [eV]')

        ax[3,0].plot(prof_data['Vz.dat'][:,0], prof_data['Vz.dat'][:,1])
        ax[3,0].set_ylabel(r'V$_z$ [cm s$^{-1}$]')
        ax[3,0].set_xlabel('r [cm]')

        ax[1,1].plot(prof_data['Vth.dat'][:,0], prof_data['Vth.dat'][:,1])
        ax[1,1].set_ylabel(r'V$_\theta$ [cm s$^{-1}$]')

        ax[2,1].plot(prof_data['Er.dat'][:,0], prof_data['Er.dat'][:,1])
        ax[2,1].set_ylabel(r'E$_r$ [statV cm$^{-1}$]')

        ax[3,1].plot(prof_data['q.dat'][:,0], prof_data['q.dat'][:,1])
        ax[3,1].set_ylabel(r'q [1]')
        ax[3,1].set_xlabel('r [cm]')

        list(map(lambda x: self.uc.add_grid_to_axis(x), ax[:,0]))
        list(map(lambda x: self.uc.add_grid_to_axis(x), ax[:,1]))

        fig.delaxes(ax[0,1])
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.15)
        plt.show()

