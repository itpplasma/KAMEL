import numpy as np
import matplotlib.pyplot as plt
from .KiLCA_interface import KiLCA_interface
import os
import re
import sys
import inspect
#sys.path.append(inspect.getfile(KiLCA_interface)[0:-18] + '../../postproc_py_class/')
from postproc_class import utility_class
from scipy.optimize import curve_fit

class KiLCA_postprocessor:
    """
    Description:
        KiLCA postprocessing class to (mainly) plot the output data.
    Methods:
        KiLCA_postprocessor(args): constructor takes either KiLCA_interface object, or the same input as the KiLCA_interface class, i.e. shot, time, runpath and run type
        read_Eb(path_to_file): reads in EB.dat output file
    """

    uc = utility_class.utility() # for colors and plotting stuff
    debug = False

    def __init__(self, *args):
        """
        Constructor of the postprocessor class.
        """
        self.m = []
        self.n =[]
        self.EBdat = {}

        if isinstance(args[0], KiLCA_interface):
            self.kil_in = args[0]
            self.shot = self.kil_in.shot
            self.time_slice = self.kil_in.time
            self.path_of_run = self.kil_in.path_of_run
            self.run_type = self.kil_in.run_type
            self.machine = self.kil_in.machine
            if self.debug:
                print('KiLCA interface was parsed')
                print('    shot      : ' + str(self.shot))
                print('    time slice: ' + str(self.time_slice))
                print('    run path  : ' + self.path_of_run)
                print('    run type  : ' + self.run_type)
                print('    machine   : ' + self.machine)
        elif isinstance(args[0], int):
            self.shot = args[0]
            self.time_slice = args[1]
            self.path_of_run = args[2] + '/' + args[3] + '/'
            self.run_type = args[3]
            if not len(args) > 4:
                self.machine = 'AUG'
            else:
                self.machine = args[4]
            if self.debug:
                print('Path was parsed')
                print('    shot      : ' + str(self.shot))
                print('    time slice: ' + str(self.time_slice))
                print('    run path  : ' + self.path_of_run)
                print('    run type  : ' + self.run_type)
                print('    machine   : ' + self.machine)
        else:
            raise ValueError(str(args[0]) + ' is not a valid input.')

        fs = os.listdir(self.path_of_run + '/linear-data/')
        self.fs = fs
        if self.debug:
            print('available modes:')
        for el in fs:
            if re.search('m_*', el):
                nums = re.findall('\d+', el)
                if el[2] == '-':
                    self.m = np.append(self.m ,-int(nums[0]))
                    if self.debug:
                        print(f'    m = -{nums[0]}    n = {nums[1]}')
                else:
                    self.m = np.append(self.m ,int(nums[0]))
                    if self.debug:
                        print(f'    m = {nums[0]}    n = {nums[1]}')
                self.n = np.append(self.n, int(nums[1]))
                
                #self.read_EB(self.path_of_run + 'linear-data/'+ el + '/')

    def read_EB(self, path_to_file='', m=0, n=0):
        """
        Description:
            Read in EB.dat file which is the linear data output of KiLCA. Returns the data as numpy array if no mode number is given and as a dict otherwise. It also safes it in the EBdat variable of the class.
        """
        self.EBdat = {}
        found_file = False
        if m==0 and n==0 and not path_to_file=='':
            self.EBdat = np.loadtxt(path_to_file + 'EB.dat')
        elif path_to_file=='':
            for f in self.fs[:]:
                nums = re.findall('\d+', f)
                if f[2] == '-':
                    st = '-' + nums[0]
                else:
                    st = nums[0]
                if str(m) == st:
                    self.EBdat[f'({m}, {n})'] = np.loadtxt(self.path_of_run + 'linear-data/'+ f + '/EB.dat')
                    found_file = True
            if not found_file:
                raise ValueError('mode number not available')

        return self.EBdat
    
    def EB_to_fields(self, EBdat, m=0, n=0):
        """
        Description:
            Takes EB dat numpy array as input an puts it into Br, Bth, Bz,... variables.
        """

        if type(EBdat)==np.ndarray:
            self.r          = EBdat[:,0]
            self.Er_real    = EBdat[:,1]
            self.Er_imag    = EBdat[:,2]
            self.Eth_real   = EBdat[:,3]
            self.Eth_imag   = EBdat[:,4]
            self.Ez_real    = EBdat[:,5]
            self.Ez_imag    = EBdat[:,6]
            self.Br_real    = EBdat[:,7]
            self.Br_imag    = EBdat[:,8]
            self.Bth_real   = EBdat[:,9]
            self.Bth_imag   = EBdat[:,10]
            self.Bz_real    = EBdat[:,11]
            self.Bz_imag    = EBdat[:,12]
        elif type(EBdat)==dict:
            self.r          = EBdat[f'({m}, {n})'][:,0]
            self.Er_real    = EBdat[f'({m}, {n})'][:,1]
            self.Er_imag    = EBdat[f'({m}, {n})'][:,2]
            self.Eth_real   = EBdat[f'({m}, {n})'][:,3]
            self.Eth_imag   = EBdat[f'({m}, {n})'][:,4]
            self.Ez_real    = EBdat[f'({m}, {n})'][:,5]
            self.Ez_imag    = EBdat[f'({m}, {n})'][:,6]
            self.Br_real    = EBdat[f'({m}, {n})'][:,7]
            self.Br_imag    = EBdat[f'({m}, {n})'][:,8]
            self.Bth_real   = EBdat[f'({m}, {n})'][:,9]
            self.Bth_imag   = EBdat[f'({m}, {n})'][:,10]
            self.Bz_real    = EBdat[f'({m}, {n})'][:,11]
            self.Bz_imag    = EBdat[f'({m}, {n})'][:,12]
        else:
            raise ValueError('Wrong type of EBdat in EB_to_fields().')

    def plot_E_field(self, m=0, n=0):
        """
        Description:
            Plots the real and imag E field components in 3 subplots for a certain mode (if multiple are given.).
        """

        found_file = False
        for f in self.fs[:]:
            if str(m) == f[2]:
                self.read_EB(self.path_of_run + 'linear-data/'+ f + '/')
                found_file = True
        if not found_file:
            raise ValueError('mode number not available')
        self.EB_to_fields(self.EBdat, m,n)

        fig, ax = plt.subplots(3, figsize=(4,6), sharex=True)

        ax[0].plot(self.r, self.Er_real, c=self.uc.col_blue, label='real')
        ax[0].plot(self.r, self.Er_imag, c=self.uc.col_yellow, label='imag')
        ax[0].set_ylabel(r'E$_r$ [statV cm$^{-1}$]')
        ax[0].legend()

        ax[1].plot(self.r, self.Eth_real, c=self.uc.col_blue, label='real')
        ax[1].plot(self.r, self.Eth_imag, c=self.uc.col_yellow, label='imag')
        ax[1].set_ylabel(r'E$_\theta$ [statV cm$^{-1}$]')
        ax[1].legend()

        ax[2].plot(self.r, self.Ez_real, c=self.uc.col_blue, label='real')
        ax[2].plot(self.r, self.Ez_imag, c=self.uc.col_yellow, label='imag')
        ax[2].set_xlabel('r [cm]')
        ax[2].set_ylabel(r'E$_z$ [statV cm$^{-1}$]')
        ax[2].legend()

        list(map(lambda x: self.uc.add_grid_to_axis(x), ax))

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        plt.show()
        return fig

    def plot_B_field(self, m=0, n=0):
        """
        Description:
            Plots the real and imag E field components in 3 subplots for a certain mode (if multiple are given.).
        """
        found_file = False
        for f in self.fs[:]:
            if str(m) == f[2]:
                self.read_EB(self.path_of_run + 'linear-data/'+ f + '/')
                found_file = True
        if not found_file:
            raise ValueError('mode number not available')
            
        self.EB_to_fields(self.EBdat, m,n)

        fig, ax = plt.subplots(3, figsize=(4,6), sharex=True)

        ax[0].plot(self.r, self.Br_real, c=self.uc.col_blue, label='real')
        ax[0].plot(self.r, self.Br_imag, c=self.uc.col_yellow, label='imag')
        ax[0].set_ylabel(r'B$_r$ [G]')
        ax[0].legend()

        ax[1].plot(self.r, self.Bth_real, c=self.uc.col_blue, label='real')
        ax[1].plot(self.r, self.Bth_imag, c=self.uc.col_yellow, label='imag')
        ax[1].set_ylabel(r'B$_\theta$ [G]')
        ax[1].legend()

        ax[2].plot(self.r, self.Bz_real, c=self.uc.col_blue, label='real')
        ax[2].plot(self.r, self.Bz_imag, c=self.uc.col_yellow, label='imag')
        ax[2].set_xlabel('r [cm]')
        ax[2].set_ylabel(r'B$_z$ [G]')
        ax[2].legend()

        list(map(lambda x: self.uc.add_grid_to_axis(x), ax))

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        plt.show()
        return fig
    
    def calculate_parallel_current_density(self, m_mode, n_mode, path_to_linear, path_to_background):
        self.read_EB(self, path_to_file=path_to_linear, m=0, n=0)
        self.read_background_field(path_to_background)

        self.kth = m_mode / self.r_eff
        self.kz = self.n_mode / self.R0

        self.Br = self.Br_real + 1j * self.Br_imag
        self.Bth = self.Bth_real + 1j * self.Bth_imag
        self.Bz = self.Bz_real + 1j * self.Bz_imag
        self.dBr = np.gradient(self.Br, self.r)
        self.dBth = np.gradient(self.Bth, self.r)
        self.dBz = np.gradient(self.Bz, self.r)

        #self.Jr = 1j / (4*np.pi) * (self.kth * self.Bz - self.kz * self.Bth)
        self.Jth = 1 /(4*np.pi) * (1j * self.kz * self.Br - self.dBz)
        self.Jz = 1 / (4 * np.pi) * (self.Bth / self.r + self.dBth - 1j * self.kth * self.Br)

        self.Jpar = (self.Jth * self.B0th + self.Jz * self.B0z) / self.B0

    def calculate_layer_width(self):
        model = lambda r, b, c: c * (1/np.sqrt(2*np.pi*b**2)) * np.exp(-(r-self.r_res)**2/(2*b**2))
        popt, popcov = curve_fit(model, self.r, self.Jpar)
        self.d = np.real(5 * np.sqrt(popt[0]))
        
    def integrate_par_current_dens(self):
        ind = np.where(self.r >= self.r_res - self.d/2 and self.r <= self.r_res + self.d/2)
        self.Ipar = np.abs(2 * np.pi * np.trapz(self.Jpar[ind] * self.r[ind], self.r[ind]))
        return self.Ipar
        
    def read_background_field(self, path):
        self.B0 = np.loadtxt(path + 'b0.dat')[:,1]
        self.B0th = np.loadtxt(path + 'b0th.dat')[:,1]
        self.B0z = np.loadtxt(path + 'b0z.dat')[:,1]