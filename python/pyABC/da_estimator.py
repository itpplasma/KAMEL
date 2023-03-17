import numpy as np
import os
from smooth2level import *
import h5py
import scipy.signal as sig
import scipy.interpolate as inter
import pandas as pd

class da_estimator:
    '''This class is used to make an estimation for the anomalous diffusion coefficient.
    It uses equilibrium data and (already preprocessed) profiles. It loads data depending
    on the input of the load method: if astra files are found in the path, astra output is
    loaded for Da. This can also be specified to be used as P in the estimation for Da
    which is called "combined" here. If no astra dat is found, it lookds for power data in
    form of 2 files for the radiated power and the input power. If nothing is found the
    "universal constant" is used.'''

    DA_UNIV = 1e4 # universal constant
    e = 1.602e-19 # electron charge in SI (to calc W to eV/s)

    profpath = ''
    fluxdatapath = ''

    type_Da = '' # type of Da: const/est/combined/astra

    r = None # small radius
    q = None # safety factor
    psi_pol_norm = None # normalized poloidal flux
    V = None # flux surface volume
    S = None # flux surface area

    ne = None # density profile
    Te = None # electron temperature profile
    dTe = None # electron temperature gradient

    Pinp = None # input power
    Prad = None # radiated power in plasma
    Peff = None # effective power on electrons

    chi_e = None # astra chi of electrons
    Q_e = None # astra Q of electrons

    Da = None # anomalous diffusion coefficient


    def __init__(self, profpath, fluxdatapath):
        
        if (not os.path.exists(profpath)): raise ValueError('Profilepath not found in ' + profpath)
        if (not os.path.exists(fluxdatapath)): raise ValueError('Fluxdatapath not found in ' + fluxdatapath)

        self.profpath = profpath
        self.fluxdatapath = fluxdatapath

        # load profiles and equi
        equilrqpsi = np.loadtxt(self.fluxdatapath + 'equil_r_q_psi.dat', skiprows=3)
        self.r = equilrqpsi[:,0] # radius
        self.q = equilrqpsi[:,1] # safety factor
        self.psi_pol_norm = equilrqpsi[:,2] / equilrqpsi[-1,2] # normalized pol flux
        self.V = equilrqpsi[:,6] # flux surface volume

        #smo = smooth2level(self.V, self.r, 2, 5)
        #self.V = smo[:,0]
        #self.S = smo[:,1]

        self.V = sig.savgol_filter(self.V, 5, 3, 0)
        self.S = np.gradient(self.V, self.r)
        

        
        # load raw profiles
        print(profpath)
        ne_raw = np.loadtxt(profpath + 'n.dat')
        Te_raw = np.loadtxt(profpath + 'Te.dat')

        self.ne = np.interp(self.psi_pol_norm, ne_raw[:,0], ne_raw[:,1])
        self.Te = np.interp(self.psi_pol_norm, Te_raw[:,0], Te_raw[:,1])
        self.dTe = np.gradient(self.Te, self.r)
    
    def load_estimation(self, shot, time, path, combine):
        # description:
        # Loads input for the estimation and calculates Da. 
        # needs shot number and time to search for files in path:
        #   1) if astra file is found (ASTRA_PB_AUG_#<shot>_t<time in s>) it will be prioritized
        # and the flag combined is used. If not, files SHOT_BPD_Prad.dat and <shot>_TOT_P_TOT.dat are searched for.
        # If found, they are loaded and the time slices are interpolated in time.
        # If no useful file is found in path, the universal constant is used
        #
        # input:
        # shot      ... shot number
        # time      ... shot time in ms
        # path      ... path of data
        # combine   ... flag if astra data is found to use pure astra or combined Da

        astraname = path + 'ASTRA_PB_AUG_#' + str(int(shot)) + '_t' + "{:.1f}".format(time/1000)
        radpowername = path + str(int(shot)) + '_BPD_Prad.dat'
        inppowername = path + str(int(shot)) + '_TOT_P_TOT.dat'

        # check if astra name exists otherwise go for powernames
        print(astraname)
        print(os.path.exists(astraname))
        if (os.path.exists(astraname)):
            print('Da from ASTRA.')
            #dat = np.genfromtxt(astraname, skip_header=1, dtype=None)
            self.dat = pd.read_csv(astraname, skiprows=1, header=None, delimiter=' ')
            self.dat = self.dat.values.transpose()

            rho_pol = self.dat[:,2]
            self.chi_e = inter.pchip_interpolate(rho_pol**2, self.dat[:,3], self.psi_pol_norm ) * 1e4
            self.Q_e = inter.pchip_interpolate(rho_pol**2, self.dat[:,5], self.psi_pol_norm) * 1e6 / self.e

            if not combine:
                self.Da = self.chi_e / 1.5
                self.type_Da = 'astra'
            else:
                self.Da = -self.Q_e / self.ne / self.dTe / self.S
                self.type_Da = 'combined'
        elif (os.path.exists(radpowername) and os.path.exists(inppowername)):
            print('Calculate Da with power estimation.')
            # load radiated power time trace
            P_rad_raw = np.loadtxt(radpowername)
            # interpolate on time and convet to eV/s
            self.Prad = np.interp(time/1000, P_rad_raw[:,0], P_rad_raw[:,1]) / self.e

            # load input power time trace
            P_inp_raw = np.loadtxt(inppowername)
            # interpolate on time and convert to eV/s
            self.Pinp = np.interp(time/1000, P_inp_raw[:,0], P_inp_raw[:,1]) / self.e

            # distribute evently to ions and electrons
            self.Peff = 0.5 * (self.Pinp - self.Prad)

            self.Da = - self.Peff / self.ne / self.dTe / self.S
            self.type_Da = 'estimation'

        else:
            print('Da from universal constant')
            self.Da = self.DA_UNIV * np.ones(len(self.r))
            self.type_Da = 'const'


    def write(self, fpath):
        # write processed data to file. Requires load_estimation before
        # input:
        # fpath ... path + name of file

        try:
            np.savetxt(fpath, np.array([self.r, self.Da]).transpose())
        except:
            raise ValueError('Cannot write Da estimation to text file')

    def export2hdf5(self, fname):
        self.writehdf5(fname, '/Da_est/r', self.r, 'small radius', 'cm')
        self.writehdf5(fname, '/Da_est/V', self.V, 'flux surface volume', 'cm^3')
        self.writehdf5(fname, '/Da_est/S', self.S, 'flux surface surface', 'cm^2')

        self.writehdf5(fname, '/Da_est/ne', self.ne, 'density profile', 'cm^{-3}')
        self.writehdf5(fname, '/Da_est/Te', self.Te, 'electron temperature profile', 'eV')
        self.writehdf5(fname, '/Da_est/dTe', self.dTe, 'radial derivative of the electron temperature profile', 'eV cm^{-1}')
        self.writehdf5(fname, '/Da_est/Da', self.Da, 'anomalous diffusion coefficient', 'cm^2 s^{-1}')

        if not (self.Pinp is None):
            self.writehdf5(fname, '/Da_est/Pinp', self.Pinp, 'input power', 'eV s^{-1}')
        if not (self.Prad is None):
            self.writehdf5(fname, '/Da_est/Prad', self.Prad, 'radiated power in plasma', 'eV s^{-1}')
        if not (self.Peff is None):
            self.writehdf5(fname, '/Da_est/Peff', self.Pinp, 'effective power on electrons', 'eV s^{-1}')
        if not (self.chi_e is None):
            self.writehdf5(fname, '/Da_est/chi_e', self.chi_e, 'astra chi of electrons', 'cm^2 s^{-1}')
        if not (self.Q_e is None):
            self.writehdf5(fname, '/Da_est/Q_e', self.Q_e, 'astra Q of electrons', 'eV s^{-1}')







    def writehdf5(self, fname, dataset, data, comment, unit):
        '''Write dataset entry into h5 file. Add comment and unit.'''
        with h5py.File(fname, 'a') as f:
            dset = f.create_dataset(dataset, data=data)
            dset.attrs['comment'] = comment
            dset.attrs['unit'] = unit