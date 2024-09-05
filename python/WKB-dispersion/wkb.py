import cxroots as cx
import numpy as np
import matplotlib.pyplot as plt

from plasmapy.dispersion import plasma_dispersion_func as plasma_disp
from wkb_grid import WKB_Grid

import os
import sys
from scipy.integrate import solve_ivp
from scipy.special import iv as Bessel
from scipy.special import gamma as Gamma
from math import factorial
from jax import grad
import jax
import jax.numpy as jnp
import h5py
import logging


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../KiLCA-QB/python/susc_functions/'))
import susc_funcs


##Define constants
pi = np.pi
e_mass = 9.1094e-28
p_mass = 1.6726e-24
e_charge = 4.8032e-10
ev = 1.6022e-12
sol = 29979245800.0
kB = 1.380649e-16


class KIM_WKB():
    ##Setup
    btor = -17977.413  # AUG config
    R0 = 165.0  # major radius in cm
    r_plas = 64.7  # plasma radius

    m_mode = 6  # poloidal mode number
    n_mode = 2  # toroidal mode number

    Zi = 1  # Ion charge number
    Ai = 2  # Ion mass number

    prof_path = ''

    bessel_large_arg_limit = 5
    contour_limit = 0 # 50 works for H, 20 for D
    
    ##Options
    options = {'prof': 'parab',
               'n_points': 50,
               'omega': 0,
               'radius': 250,
               'r_per': 0.9,
               'r_range_start': 10,
               'omegaOfk': False,
               'kOfr': True,
               'mode': 'KIM',
               'int_method': 'quad',
               'der': False,
               'save': True,
               'log': False,
               'new_grid': False,
               'number_of_roots_to_find': 2,
               'Collisions': 'collisionless', # Possible: collisionless, Krook, FokkerPlanck
               'max_cyclotron_harmonic': 0}
    options['omega_range'] = np.linspace(0, 50, options['n_points'])

    possible_coll_models = ['collisionless', 'Krook', 'FokkerPlanck']

    possible_operation_modes = ['KIM', 'horton']

    species = ['e', 'D']
    ion_species = ['D']

    general_dat = {}
    general_dat_keys = ['r', 'Er', 'ks', 'kp', 'om_E', 'q', 'prof_length']
    general_dat_units = {'r': 'cm', 'Er': 'statV/cm', 'ks': '1/cm', 'kp': '1/cm', 'om_E': '1/s', 'q': '1', 'prof_length': '1'}

    spec_keys = ['n', 'T', 'vT', 'nu', 'mass', 'Ai', 'charge', 'Zi']
    spec_dat_units = {'n': 'cm', 'T': 'erg', 'vT': 'cm/s', 'nu': '1/s', 'mass': 'g', 'Ai': '1', 'charge': 'statC', 'Zi': '1'}
    spec_dat = {}

    equil_keys = ['u', 'B0z', 'B0th', 'B0', 'hz', 'hth']
    equil_dat_units = {'u': 'G^2', 'B0z': 'G', 'B0th': 'G', 'B0': 'G', 'hz': '1', 'hth': '1'}
    equil_dat = {}

    def __init__(self, species=['e', 'D'], spec_mass = [e_mass, 2*p_mass], spec_charge_num = [-1, 1]):
        self.species = species
        for spec in self.species:
            self.spec_dat[spec] = {}
            if spec != 'e':
                self.spec_dat[spec]['Ai'] = spec_mass[self.species.index(spec)] / p_mass
            self.spec_dat[spec]['mass'] = spec_mass[self.species.index(spec)]

            self.spec_dat[spec]['charge'] = spec_charge_num[self.species.index(spec)] * e_charge
            self.spec_dat[spec]['Zi'] = spec_charge_num[self.species.index(spec)]
        self.ion_species = [spec for spec in self.species if spec != 'e']

        self.set_plot_options()

    def set_plot_options(self):
        plt.rcParams.update({"font.size": 12})
        plt.rcParams.update({"lines.markersize": 10})

    def load_prof(self, name, mode="const"):
        loadR = False
        if name == "r":
            loadR = True
            name = "Er.dat"
        if mode == "const":
            with open(os.path.join("profiles_const", name), "r") as file:
                data = file.readlines()
        elif mode == "parab":
            with open(os.path.join("profiles_parab", name), "r") as file:
                data = file.readlines()
        else:
            with open(os.path.join("profiles", name), "r") as file:
                data = file.readlines()
        out = np.zeros(len(data))
        if loadR:
            for k in range(len(data)):
                if "\t" in data[k]:
                    line = data[k].split("\t")
                else:
                    line = data[k].split(" ")
                out[k] = float(line[0])
            return out
        for k in range(len(data)):
            if "\t" in data[k]:
                line = data[k].split("\t")
            else:
                line = data[k].split(" ")
            out[k] = float(line[1].split("\n")[0])
        return out


    def calcEquilibrium(self):
        # Pressure due to electrons
        press_prof = np.zeros(len(self.general_dat['r']))
        for spec in self.species:
            press_prof += ev * self.spec_dat[spec]['n'] * self.spec_dat[spec]['T']
        dpress_prof = np.gradient(press_prof, self.general_dat['r'])
        # solve ivp
        f = lambda r,u: -2*r*u/(self.R0**2 * self.general_dat['q'][self.val2ind(r, self.general_dat['r'])]**2 + r**2)-8*pi*dpress_prof[self.val2ind(r, self.general_dat['r'])]
        u0 = self.btor**2 * (1 + self.general_dat['r'][0] ** 2 / (self.R0**2 * self.general_dat['q'][0] ** 2))
        u = solve_ivp(f, t_span=[self.general_dat['r'][0], self.general_dat['r'][-1]], y0=[u0], t_eval=self.general_dat['r'])
        u = u.y[0]
        self.equil_dat['B0z'] = np.sign(self.btor) * np.sqrt(u / (1 + self.general_dat['r']**2 / (self.R0**2 * self.general_dat['q']**2)))
        self.equil_dat['B0th'] = self.equil_dat['B0z'] * self.general_dat['r'] / (self.general_dat['q'] * self.R0)
        self.equil_dat['B0'] = np.sqrt(self.equil_dat['B0th']**2 + self.equil_dat['B0z']**2)
        self.equil_dat['hz'] = self.equil_dat['B0z'] / self.equil_dat['B0']
        self.equil_dat['hth'] = self.equil_dat['B0th'] / self.equil_dat['B0']

    
    def calc_parameters(self):
        self.calcEquilibrium()
        if self.options['Collisions'] == 'collisionless':
            for spec in self.species:
                self.spec_dat[spec]['nu'] = np.zeros(self.general_dat['prof_length'])
        else:
            self.calc_collision_frequency()
        self.calc_all_derivs()

        self.general_dat['ks'] = (self.m_mode * self.equil_dat['hz'] - self.n_mode * self.equil_dat['hth'] / self.R0) / self.general_dat['r'] 
        self.general_dat['kp'] = self.m_mode / self.general_dat['r'] * self.equil_dat['hth'] + self.n_mode / self.R0 * self.equil_dat['hz'] 
        self.general_dat['om_E'] = -sol * self.general_dat['ks'] * self.general_dat['Er'] / self.equil_dat['B0'] 

        for spec in self.species:
            self.spec_dat[spec]['vT'] = np.sqrt(self.spec_dat[spec]['T'] * ev / self.spec_dat[spec]['mass'])
            self.spec_dat[spec]['omega_c'] = np.abs(self.spec_dat[spec]['charge']) * self.equil_dat['B0'] / (self.spec_dat[spec]['mass'] * sol)
            self.spec_dat[spec]['lambda_D'] = np.sqrt(self.spec_dat[spec]['T'] * ev / (4.0 * np.pi * self.spec_dat[spec]['n'] * self.spec_dat[spec]['charge']**2))
            self.spec_dat[spec]['A1'] = self.spec_dat[spec]['dndr'] / self.spec_dat[spec]['n'] \
                - self.spec_dat[spec]['charge'] / (self.spec_dat[spec]['T'] * ev) * self.general_dat['Er'] \
                - 3 / (2 * self.spec_dat[spec]['T']) * self.spec_dat[spec]['dTdr']
            self.spec_dat[spec]['A2'] = self.spec_dat[spec]['dTdr'] / self.spec_dat[spec]['T']


            self.spec_dat[spec]['z0'] = -(self.general_dat['om_E'] - self.options['omega'] - 1j * self.spec_dat[spec]['nu']) \
                / (np.abs(self.general_dat['kp']) * np.sqrt(2) * self.spec_dat[spec]['vT'])
            self.spec_dat[spec]['rho_TL'] = self.spec_dat[spec]['vT'] / np.abs(self.spec_dat[spec]['omega_c'])

        
    def calc_collision_frequency(self):

        Lee = 23.5 - np.log(np.sqrt(self.spec_dat['e']['n']) * self.spec_dat['e']['T'] ** (-5 / 4)) \
            - np.sqrt(1e-5 + ((np.log(self.spec_dat['e']['T']) - 2) ** 2) / 16)
        self.spec_dat['e']['nu'] = 5.8e-6 * self.spec_dat['e']['n'] * Lee * self.spec_dat['e']['T'] ** (-3 / 2)  # Collision frequency
        Lei = 24 - np.log(np.sqrt(self.spec_dat['e']['n']) / self.spec_dat['e']['T'])

        # test particle (first index) collides with field particles (second index)
        for test in self.ion_species:
            # nue ei: collision frequency electrons with ions
            self.spec_dat['e']['nu'] = self.spec_dat['e']['nu'] + 7.7e-6 * self.spec_dat['e']['T']**-1.5 * self.spec_dat[test]['n'] * Lei * self.spec_dat[test]['Zi']**2
            # nui ie: collision frequency ions with electrons
            self.spec_dat[test]['nu'] = 1.8e-7 * self.spec_dat[test]['Ai'] ** (-1 / 2) * self.spec_dat[test]['T'] ** (-3 / 2) * self.spec_dat['e']['n'] * self.spec_dat[test]['Zi'] ** 2 * Lei
            for field in self.ion_species:
                Liiprime = 23 - np.log(self.spec_dat[test]['Zi']*self.spec_dat[field]['Zi'] * (self.spec_dat[test]['Ai'] + self.spec_dat[field]['Ai']) \
                    / (self.spec_dat[test]['T'] * self.spec_dat[field]['Ai'] + self.spec_dat[field]['T'] * self.spec_dat[test]['Ai']) \
                    * np.sqrt(self.spec_dat[test]['n'] * self.spec_dat[test]['Zi'] ** 2 / self.spec_dat[test]['T'] \
                    + self.spec_dat[field]['n'] * self.spec_dat[field]['Zi'] ** 2 / self.spec_dat[field]['T'] ))

                self.spec_dat[test]['nu'] = self.spec_dat[test]['nu'] + 1.8e-7 * self.spec_dat[test]['Ai']**-0.5 * self.spec_dat[test]['T']**-1.5 * self.spec_dat[field]['n'] * self.spec_dat[test]['Zi']**2 * self.spec_dat[field]['Zi']**2 * Liiprime
                
        return [self.spec_dat[spec]['nu'] for spec in self.species]



    # For convenience
    def val2ind(self, value, array):
        idx = np.argmin(np.abs(array - value))
        return idx

    def save2txt(self, r_or_omega, k1, k2, name):
        if self.options['kOfr']:
            with open(name, "w") as file:
                file.write("radius [cm]; Re(k1); Im(k1); Re(k2); Re(k2)\n")
                for k in range(len(r_or_omega)):
                    file.write(
                        "%f; %f; %f; %f; %f\n"
                        % (r_or_omega[k], k1.real[k], k1.imag[k], k2.real[k], k2.imag[k])
                    )
        else:
            with open(name, "w") as file:
                file.write("omega [1/s]; Re(k1); Im(k1); Re(k2); Re(k2)\n")
                for k in range(len(r_or_omega)):
                    file.write(
                        "%f; %f; %f; %f; %f\n"
                        % (r_or_omega[k], k1.real[k], k1.imag[k], k2.real[k], k2.imag[k])
                    )
        return


    def load_txt(self, name):
        with open(name, "r") as file:
            header = file.readline()
            data = file.readlines()
        if "radius" in header:
            r_or_omega = list()
            k1 = list()
            k2 = list()
            for line in data:
                r_or_omega.append(float(line.split(";")[0]))
                k1.append(float(line.split(";")[1]) + 1j * float(line.split(";")[2]))
                k2.append(float(line.split(";")[3]) + 1j * float(line.split(";")[4]))
        else:
            r_or_omega = list()
            k1 = list()
            k2 = list()
            for line in data:
                r_or_omega.append(float(line.split(";")[0]))
                k1.append(float(line.split(";")[1]) + 1j * float(line.split(";")[2]))
                k2.append(float(line.split(";")[3]) + 1j * float(line.split(";")[4]))
        return np.array(r_or_omega), np.array(k1), np.array(k2)

    def load_all_profs(self):
        try:
            self.r_prof = np.loadtxt(self.prof_path + 'n.dat')[:,0]
        except:
            raise ValueError("No n.dat found in path " + self.prof_path)
        self.E_prof = np.loadtxt(self.prof_path + 'Er.dat')[:,1]
        self.n_prof = np.loadtxt(self.prof_path + 'n.dat')[:,1]
        self.q_prof = np.loadtxt(self.prof_path + 'q.dat')[:,1]
        self.Te_prof = np.loadtxt(self.prof_path + 'Te.dat')[:,1]
        self.Ti_prof = np.loadtxt(self.prof_path + 'Ti.dat')[:,1]
        self.prof_length = len(self.r_prof)

        self.spec_dat['e']['n'] = self.n_prof
        self.spec_dat['e']['T'] = self.Te_prof
        #self.spec_dat['e']['vT'] = np.sqrt(self.spec_dat['e']['T'] * ev / self.spec_dat['e']['mass'])
        self.general_dat['Er'] = self.E_prof
        self.general_dat['q'] = self.q_prof
        self.general_dat['r'] = self.r_prof
        self.general_dat['prof_length'] = len(self.r_prof)

        for spec in self.ion_species:
            self.spec_dat[spec]['T'] = self.Ti_prof

        self.set_ion_densities()

    def set_ion_densities(self):
        # TODO: needs correct implementation for more than 1 ion species
        for spec in self.species:
            if spec != 'e':
                self.ni_prof = self.n_prof
                self.spec_dat[spec]['n'] = self.ni_prof

    
    def calc_all_derivs(self):
        self.calc_species_derivs()
        self.calc_general_derivs()
    
    def calc_species_derivs(self):
        derivs_in_species_for = ['n', 'T']
        for spec in self.species:
            for deriv in derivs_in_species_for:
                self.spec_dat[spec][f'd{deriv}dr'] = np.gradient(self.spec_dat[spec][deriv], self.general_dat['r'])
    
    def calc_general_derivs(self):
        derivs_in_general_for = ['q']
        for deriv in derivs_in_general_for:
            self.general_dat[f'd{deriv}dr'] = np.gradient(self.general_dat[deriv], self.general_dat['r'])
        
    def find_res_surface(self):
        # Find resonant surface
        self.idx_res = self.val2ind(-self.m_mode / self.n_mode, self.general_dat['q'])
        self.res_surf_val = np.interp(self.m_mode / self.n_mode, np.abs(self.general_dat['q']), self.general_dat['r'])


    def calc_dispersion_equation_for_kr_values_at_r(self, kr, position, mode="KIM"):
        if type(position) == int:
            r_eval = self.general_dat['r'][position]
        if type(kr) == complex or type(kr) == np.complex128:
            #disp_prof = self.create_dispersion_equation_profile(kr, mode)
            dispersion_equation = self.create_dispersion_equation_single(kr, position, mode)
            #dispersion_equation = np.interp(r_eval, self.general_dat['r'], disp_prof)
        else:
            raise ValueError("Other types of kr not implemented yet")
        return dispersion_equation
        

    def create_dispersion_equation_single(self, kr, r_indx, mode="KIM"):
        assert mode in self.possible_operation_modes, f"Mode {mode} not supported"

        self.general_dat['kperp'] = np.sqrt(self.general_dat['ks']**2 + kr**2)
        if mode == 'KIM':
            dispersion_equation = self.calc_dispersion_equation_KIM(kr, r_indx)
        elif mode == "horton":
            dispersion_equation = self.calc_dispersion_equation_horton_single(kr, r_indx)
        return dispersion_equation


    def calc_dispersion_equation_KIM(self, kr, r_indx):
        if self.options['Collisions'] == 'collisionless' or self.options['Collisions'] == 'Krook':
            return self.calc_dispersion_equation_KIM_Krook(kr, r_indx)
        elif self.options['Collisions'] == 'FokkerPlanck':
            return self.calc_dispersion_equation_KIM_FokkerPlanck(kr, r_indx)
        else:
            raise ValueError(f'Collision model {self.options["Collisions"]} not supported')
        

    def calc_dispersion_equation_KIM_Krook(self, kr, r_indx):
        dispersion_equation = self.general_dat['kperp'][r_indx]**2 + self.general_dat['kp'][r_indx]**2
        for spec in self.species:
            eval_b = self.general_dat['kperp'][r_indx]**2 * self.spec_dat[spec]['rho_TL'][r_indx]**2
            BesselProd0, BesselProd1 = self.calc_needed_bessel_single(eval_b)
            if np.isnan(BesselProd0) or np.isnan(BesselProd1):
                print(f'BesselProd0 or BesselProd1 contains NaNs')

            dispersion_equation += self.spec_dat[spec]['lambda_D'][r_indx]**-2 * \
                (
                    1.0 - self.general_dat['ks'][r_indx] * self.spec_dat[spec]['rho_TL'][r_indx] \
                / (self.general_dat['kp'][r_indx] * np.sqrt(2)) * (self.spec_dat[spec]['A1'][r_indx] * BesselProd0 * plasma_disp(self.spec_dat[spec]['z0'][r_indx]) + \
                self.spec_dat[spec]['A2'][r_indx] * (plasma_disp(self.spec_dat[spec]['z0'][r_indx]) * (1 + eval_b + self.spec_dat[spec]['z0'][r_indx]**2)\
                * np.exp(-eval_b) + BesselProd1 * eval_b + self.spec_dat[spec]['z0'][r_indx] * BesselProd0))
                )
        #print(dispersion_equation)
        return dispersion_equation

    def calc_needed_bessel_single(self, eval_b):
        #if np.real(eval_b)>self.bessel_large_arg_limit:
            #BesselProd0 = 1.0 / (np.sqrt(2*np.pi*eval_b)) + 0j
            #BesselProd1 = 1.0 / (np.sqrt(2*np.pi*eval_b)*(1 + 1 / eval_b**2)**(1/4)) * np.exp(np.arcsinh(-1 / eval_b) + eval_b \
                #* np.sqrt(1 + 1 / eval_b**2) - eval_b) + 0j
        #else:
            #BesselProd0 = Bessel(0, eval_b)*np.exp(-eval_b) + 0j
            #BesselProd1 = Bessel(-1, eval_b)*np.exp(-eval_b) + 0j

        return self.calc_needed_bessel_of_mphi( 0, eval_b)
        #return BesselProd0, BesselProd1
    

    def calc_dispersion_equation_KIM_FokkerPlanck(self, kr, r_indx):
        dispersion_equation = self.general_dat['kperp'][r_indx]**2 + self.general_dat['kp'][r_indx]**2
        #print(dispersion_equation)
        for spec in self.species:
            eval_b = self.general_dat['kperp'][r_indx]**2 * self.spec_dat[spec]['rho_TL'][r_indx]**2
            for m_phi in range(0, self.options['max_cyclotron_harmonic']+1):
                BesselProd0, BesselProd1 = self.calc_needed_bessel_of_mphi(mphi=m_phi, eval_b=eval_b)
                if np.isnan(BesselProd0) or np.isnan(BesselProd1):
                    print(f'FokkerPlanck: BesselProd0 or BesselProd1 contains NaNs, eval_b = {eval_b}, m_phi = {m_phi}')
                    print(f'BesselProd0: {BesselProd0}, BesselProd1: {BesselProd1}')

                #dispersion_equation_temp = \
                #    (1j * self.spec_dat[spec]['vT'][r_indx]**2.0 * self.general_dat['ks'][r_indx]) / (self.spec_dat[spec]['lambda_D'][r_indx]**2 \
                #        * self.spec_dat[spec]['omega_c'][r_indx] * self.spec_dat[spec]['nu'][r_indx]) \
                        #* (self.spec_dat[spec]['I00'][m_phi][r_indx] * (BesselProd0 * (self.spec_dat[spec]['A1'][r_indx] + \
                            #self.spec_dat[spec]['A2'][r_indx] * (1.0 + eval_b + m_phi)) + self.spec_dat[spec]['A2'][r_indx] * eval_b * BesselProd1)\
                        #+ self.spec_dat[spec]['I20'][m_phi][r_indx] * 0.5 * self.spec_dat[spec]['A2'][r_indx]*BesselProd0)
                #print(dispersion_equation_temp)
                
                dispersion_equation_temp = \
                   (1j * self.spec_dat[spec]['vT'][r_indx]**2.0 * self.general_dat['ks'][r_indx]) \
                       / (self.spec_dat[spec]['lambda_D'][r_indx]**2 * self.spec_dat[spec]['omega_c'][r_indx] * self.spec_dat[spec]['nu'][r_indx]) \
                    * (
                        self.spec_dat[spec]['A1'][r_indx] * self.spec_dat[spec]['I00'][m_phi][r_indx] * BesselProd0 \
                        + self.spec_dat[spec]['A2'][r_indx] \
                        * (
                           self.spec_dat[spec]['I00'][m_phi][r_indx] * (
                                 BesselProd0 * (1 + eval_b + m_phi) + eval_b * BesselProd1
                            ) 
                           + 0.5 * self.spec_dat[spec]['I20'][m_phi][r_indx] * BesselProd0
                        )
                    )
                
                dispersion_equation -= dispersion_equation_temp
        #print(dispersion_equation)
        
        return dispersion_equation

    #def calc_needed_bessel_of_mphi(self, mphi, eval_b):
        #if np.abs(np.real(eval_b))>self.bessel_large_arg_limit:
            ## asymptotic form for Bessel function:

            #BesselProd0 = 1.0 / (np.sqrt(2*np.pi*eval_b + 0j) * (1+ mphi**2 / eval_b**2)**(1/4)) * np.exp(- mphi * np.arcsinh(mphi / eval_b) + eval_b \
                #* np.sqrt(1 + mphi**2 / eval_b**2) - eval_b) + 0j
            #BesselProd1 = 1.0 / (np.sqrt(2*np.pi*eval_b + 0j) * (1+ (mphi-1)**2 / eval_b**2)**(1/4)) * np.exp(- (mphi-1) * np.arcsinh((mphi-1) / eval_b) + eval_b \
                #* np.sqrt(1 + (mphi-1)**2 / eval_b**2) - eval_b) + 0j
        #else:
            #BesselProd0 = Bessel(mphi, eval_b)* np.exp(-eval_b)
            #BesselProd1 = Bessel(mphi-1, eval_b)* np.exp(-eval_b)
        #return BesselProd0, BesselProd1 
    
    def calc_needed_bessel_of_mphi(self, mphi, eval_b):
        
        #if np.abs(eval_b)>self.bessel_large_arg_limit:
        if True:
            # use definition via sum
            BesselProd0 = 0 + 0j
            BesselProd1 = 0 + 0j
            BesselProd0_store = 1.0 + 0j
            BesselProd1_store = 1.0 + 0j
            k = 0
            while (np.abs(BesselProd0 - BesselProd0_store) > 1e-10 or np.abs(BesselProd1 - BesselProd1_store) > 1e-10):
                BesselProd0_store = 1.0 * BesselProd0
                BesselProd1_store = 1.0 * BesselProd1
                BesselProd0 += (0.5 * eval_b)**(2*k) / (factorial(k) * Gamma(k+mphi+1) + 0j) * np.exp(-eval_b)
                BesselProd1 += (0.5 * eval_b)**(2*k) / (factorial(k) * Gamma(k+mphi) + 0j) * np.exp(-eval_b)
                k +=1
                if k > 50:
                    break
            BesselProd0 *= (0.5 * eval_b)**mphi
            BesselProd1 *= (0.5 * eval_b)**(mphi-1)
        #else:
        #    BesselProd0 = Bessel(mphi, eval_b)* np.exp(-eval_b)
        #    BesselProd1 = Bessel(mphi-1, eval_b)* np.exp(-eval_b)
        
        return BesselProd0, BesselProd1 
    
    def test_bessel_of_mphi(self):
        mphi_list = [0]
        eval_b = np.linspace(-7,7,100)
        plt.figure(figsize=(8,6))
        Bess0 = np.zeros(len(eval_b), dtype=complex)
        Bess1 = np.zeros(len(eval_b), dtype=complex)
        self.bessel_large_arg_limit = 1
        for i, mphi in enumerate(mphi_list):
            for j, b in enumerate(eval_b):
                Bess0[j], Bess1[j] = self.calc_needed_bessel_of_mphi(mphi, b)
            plt.plot(eval_b, np.real(Bess0), label=f'Re, Bessel0, mphi = {mphi}', marker='x')
            plt.plot(eval_b, np.real(Bess1), label=f'Re, Bessel1, mphi = {mphi}', marker='x')
            plt.plot(eval_b, np.imag(Bess0), label=f'Im, Bessel0, mphi = {mphi}', ls=':', marker='.')
            plt.plot(eval_b, np.imag(Bess1), label=f'Im, Bessel1, mphi = {mphi}', ls=':', marker='.')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.grid()
        #plt.yscale('log')
        plt.show()
        
        

    def plot_dispersion_equation_KIM_FokkerPlanck(self):

        self.load_all_profs()
        self.calc_parameters()
        self.find_res_surface()
        self.calc_needed_susc_funcs()

        kr = 1j*np.linspace(-100, 100, 100)
        r_indx = 0
        
        disp_eq = np.zeros(len(kr), dtype=complex)
        for i,k in enumerate(kr):
            self.general_dat['kperp'] = np.sqrt(self.general_dat['ks']**2 + k**2)
            disp_eq[i] = self.calc_dispersion_equation_KIM_FokkerPlanck(k, r_indx)
        print(disp_eq)
        plt.figure()
        plt.plot(np.imag(kr), np.real(disp_eq))
        plt.plot(np.imag(kr), np.imag(disp_eq))
        plt.grid()
        plt.show()

    def calc_dispersion_equation_horton_single(self, kr, r_indx):
        ky = self.general_dat['ks'][r_indx]
        om_prime = self.options['omega'] - self.general_dat['om_E'][r_indx]
        
        dispersion_equation = 0 + 0j
        denom = 0 + 0j
        nom = 0 + 0j
        for spec in self.species:
            self.spec_dat[spec]['om_n'] = ky * sol / (self.spec_dat[spec]['charge'] * self.equil_dat['B0'][r_indx] \
                * self.spec_dat[spec]['n'][r_indx]) * self.spec_dat[spec]['T'][r_indx] * ev * self.spec_dat[spec]['dndr'][r_indx]
            self.spec_dat[spec]['om_T'] = ky * sol / (self.spec_dat[spec]['charge'] * self.equil_dat['B0'][r_indx]) * self.spec_dat[spec]['dTdr'][r_indx] * ev

            self.spec_dat[spec]['z'] = om_prime / (np.sqrt(2) * self.general_dat['kp'][r_indx] * self.spec_dat[spec]['vT'][r_indx])
            self.spec_dat[spec]['W'] = (self.spec_dat[spec]['om_n'] - om_prime + self.spec_dat[spec]['om_T'] * (self.spec_dat[spec]['z']**2 + 0.5)) \
                * self.spec_dat[spec]['z'] / om_prime * plasma_disp(self.spec_dat[spec]['z']) + self.spec_dat[spec]['om_T'] * self.spec_dat[spec]['z']**2 / om_prime

            self.spec_dat[spec]['nom'] = self.spec_dat[spec]['lambda_D'][r_indx]**-2 * (self.spec_dat[spec]['W'] - 1.0 \
                -self.spec_dat[spec]['om_T'] * self.spec_dat[spec]['z'] / om_prime * plasma_disp(self.spec_dat[spec]['z']))
            self.spec_dat[spec]['denom'] = self.spec_dat[spec]['W'] * self.spec_dat[spec]['vT'][r_indx]**2 / (self.spec_dat[spec]['omega_c'][r_indx]**2 \
                * self.spec_dat[spec]['lambda_D'][r_indx]**2)

            nom += self.spec_dat[spec]['nom']
            denom += self.spec_dat[spec]['denom']

        nom -= self.general_dat['kp'][r_indx]**2
        dispersion_equation = nom / (1 + denom) - self.general_dat['kperp'][r_indx]**2
        
        return dispersion_equation

    def initialize_data(self):
        self.load_all_profs()
        #self.calc_all_derivs()
        self.find_res_surface()
        if self.options['new_grid']:
            self.generate_non_uniform_grid()
            self.interpolate_on_non_uniform_grid()
        self.calc_parameters()

        if self.options['Collisions'] == 'FokkerPlanck':
            #if self.options['max_cyclotron_harmonic'] > 0:
                #raise ValueError('Higher cylcotron harmonic than lowest order not yet implemented!')
            self.calc_needed_susc_funcs()
     
    def set_output_h5_file(self, file, append_or_write='w'):
        self.h5f = file
        self.h5_append_or_write = append_or_write
    
    def calc_dispersion_relation_k_of_r(self, mode):   

        self.set_model_mode(mode)
        self.initialize_data()

        idx = int(self.general_dat['prof_length'] * self.options['r_per'])
        idx_range = np.linspace(self.options['r_range_start'], self.general_dat['prof_length'] - 1, self.options['n_points'])

        res = []
        r_used = []
        res_sorted = []
         
        contour=cx.Rectangle(x_range=[-self.contour_limit,self.contour_limit], y_range=[-self.contour_limit,self.contour_limit])

        if self.options['log']:
            logging.basicConfig(level=logging.INFO)

        for r in idx_range:
            roots = self.find_roots(int(r), contour)
            if roots == 1:
                continue
            print(f"Found roots for index {int(r)} at {self.general_dat['r'][int(r)]}:")
            print(f"{roots}")
            #sorted_roots = self.sort_roots(roots)
            res.append(roots)
            r_used.append(self.general_dat['r'][int(r)])
            #res_sorted.append(sorted_roots)
        
        self.write_roots_to_h5(res, r_used)
        self.r_found, self.k_r1, self.k_r2 = self.store_roots_in_variables(res, r_used)

        #self.save_found_kr_of_r(mode)
        
    def find_roots(self, r_ind, contour):
        #equation_k = lambda k: self.createDispersionEquation(k, self.options['omega'], int(r), mode=self.options['mode'])
        equation_k = lambda k: self.calc_dispersion_equation_for_kr_values_at_r(k, r_ind, mode=self.options['mode'])
        iterations=0
        roots_number = 0
        search_scale_bigger = 2.0
        search_scale_smaller = 1.5
        iteration = 0
        try:
        #if True:
            #while roots_number != self.options['number_of_roots_to_find']:
            while roots_number > self.options['number_of_roots_to_find'] or roots_number < 2:
                roots_number=contour.count_roots(equation_k)
                print(f"roots found: {roots_number}")
                if roots_number == 0:
                    contour=cx.Rectangle(np.array(contour.x_range)*search_scale_bigger,np.array(contour.y_range)*search_scale_bigger)
                if roots_number > self.options['number_of_roots_to_find']:
                    contour=cx.Rectangle(np.array(contour.x_range)/search_scale_smaller,np.array(contour.y_range)/search_scale_smaller)
                iteration += 1
                if iteration > 100:
                    print(f"Abort. Too many iterations. Did not find exactly {self.options['number_of_roots_to_find']} roots")
                    return 1
        except:
            print(f"Error in count_roots for r_ind = {r_ind}")
            return 1
        if self.options['der']:
            roots = contour.roots(equation_k, df=lambda k: complex(grad(equation_k, holomorphic=True)(k)), guess_roots_symmetry=lambda z: [-z],verbose=True)
        else:
            try:
                roots = contour.roots(equation_k,guess_roots_symmetry=lambda z: [-z],int_method=self.options['int_method'],verbose=True,int_abs_tol=0.1,root_err_tol=1e-3)
                #roots = contour.roots(equation_k,int_method=self.options['int_method'],verbose=True,int_abs_tol=0.1,root_err_tol=1e-3)
            except:
                print(f"Error in roots for r = {r_ind}")
                return 1
        return roots

    def sort_roots(self, roots):
        roots = np.array(roots.roots)
        idx = np.argsort(np.real(roots))
        return roots[idx]


    def write_roots_to_h5(self, res, r_used):
        h5f = h5py.File(self.h5f, self.h5_append_or_write)
        for i, r in enumerate(r_used):
            try:
                h5f.create_dataset(f'roots/r_{r}', data=res[i].roots)
            except:
                print(f"Error in writing roots for r = {r}")
        h5f.close()
    
    def store_roots_in_variables(self, res, r_used):
        k_r1 = list()
        k_r2 = list()
        r_found = list()
        for k in range(len(res)):
            for root in res[k].roots:
                k_r1.append(root)
                k_r2.append(root)
                r_found.append(r_used[k])
            #if len(res[k])==2:
                #print(f"Two roots found for r = {r_used[k]}")
                #print(f'Roots: {res[k]}')

                    
                #if res[k].roots[0].real > 0:
                    #k_r1.append(res[k][0])
                    #k_r2.append(res[k][1])
                    #r_found.append(r_used[k])
                #else:
                    #k_r1.append(res[k][1])
                    #k_r2.append(res[k][0])
                    #r_found.append(r_used[k])
            #if len(res[k])>2:
                ##print(f"More than two roots found for r = {r_used[k]}")
                #if res[k].roots[0].real > 0:
                    #k_r1.append(res[k][0])
                    #k_r2.append(res[k][1])
                    #r_found.append(r_used[k])
                #else:
                    #k_r1.append(res[k][1])
                    #k_r2.append(res[k][0])
                    #r_found.append(r_used[k])
        return np.array(r_found), np.array(k_r1), np.array(k_r2)
    
    def calc_needed_susc_funcs(self):
        for spec in self.species:
            self.spec_dat[spec]['x1'] = self.general_dat['kp'] * self.spec_dat[spec]['vT'] / self.spec_dat[spec]['nu']
            self.spec_dat[spec]['I00'] = {}
            self.spec_dat[spec]['I20'] = {}
            for m_phi in range(0, self.options['max_cyclotron_harmonic']+1):
                self.spec_dat[spec]['x2'] = - (m_phi * self.spec_dat[spec]['omega_c'] + self.general_dat['om_E'] - self.options['omega']) \
                    / self.spec_dat[spec]['nu']
                self.spec_dat[spec]['I00'][m_phi] = np.zeros(self.general_dat['prof_length'], dtype=complex)
                self.spec_dat[spec]['I20'][m_phi] = np.zeros(self.general_dat['prof_length'], dtype=complex)
                for i in range(0, self.general_dat['prof_length']):
                    I = susc_funcs.calc_imn_array(self.spec_dat[spec]['x1'][i], self.spec_dat[spec]['x2'][i])
                    self.spec_dat[spec]['I00'][m_phi][i] = I[0,0]
                    self.spec_dat[spec]['I20'][m_phi][i] = I[2,0]

    def plot_x2(self):
        self.initialize_data()
        plt.figure()
        for spec in self.species:
            plt.plot(self.general_dat['r'], np.real(self.spec_dat[spec]['I00'][0]), label=spec, marker='.')
            #plt.plot(self.general_dat['r'], np.imag(self.spec_dat[spec]['I00'][0]), label=spec,ls='--')
            #plt.plot(self.general_dat['r'], np.real(self.spec_dat[spec]['I20'][0]), label=spec)
            #plt.plot(self.general_dat['r'], np.imag(self.spec_dat[spec]['I20'][0]), label=spec, ls='--')
        plt.legend()
        plt.show()
    
    def test_susc_funcs(self):
        self.load_all_profs()
        self.calc_parameters()
        self.calc_needed_susc_funcs()
        
        #plt.figure()
        #for spec in self.species:
            #plt.plot(self.general_dat['r'], self.spec_dat[spec]['x1'], label=spec)
        #plt.plot(self.general_dat['r'], self.spec_dat[spec]['x2'])
        plt.figure()
        for spec in self.species: 
            if spec == 'e':
                continue
            for m_phi in range(0, self.options['max_cyclotron_harmonic']+1):
                plt.plot(self.general_dat['r'], self.spec_dat[spec]['I00'][m_phi].real, label=spec)
        plt.legend()
        #plt.yscale('log')
        plt.show()
         

    def save_found_kr_of_r(self, mode):
        file_name = mode + "_" + self.options['prof'] + "_" + "k_of_r"
        if self.options['der']:
            self.save2txt(self.r_found, self.k_r1, self.k_r2, file_name + "_jax" + ".txt",)
        else:
            self.save2txt(self.r_found, self.k_r1, self.k_r2, file_name + ".txt")

    def plot_kr_of_r(self):
        plt.figure()
        ax = plt.gca()
        ax.tick_params(right=True, top=True, direction="in")
        plt.grid()
        plt.plot(self.r_found, self.k_r1.real, ".k")
        plt.plot(self.r_found, self.k_r2.real, ".k")
        plt.axvline(
            self.general_dat['r'][self.idx_res],
            color="dimgrey",
            linestyle="dashed",
            label="Resonant surface",
        )
        plt.title("Real")
        plt.xlabel(r"$r$ [cm]")
        plt.ylabel(r"$\Re(k_r)$ [cm$^{-1}$]")

        plt.figure()
        ax = plt.gca()
        ax.tick_params(right=True, top=True, direction="in")
        plt.grid()
        plt.plot(self.r_found, self.k_r1.imag, ".k")
        plt.plot(self.r_found, self.k_r2.imag, ".k")
        plt.axvline(
            self.general_dat['r'][self.idx_res],
            color="dimgrey",
            linestyle="dashed",
            label="Resonant surface",
        )
        plt.title("Imaginary")
        plt.xlabel(r"$r$ [cm]")
        plt.ylabel(r"$\Im(k_r)$ [cm$^{-1}$]")
        plt.show()
    
    def plot_from_file(self, file):
        self.load_all_profs()
        self.find_res_surface()
        dat = np.loadtxt(file, skiprows=1, delimiter=";")
        self.r_found = dat[:,0]
        self.k_r1 = dat[:,1] + 1j*dat[:,2]
        self.k_r2 = dat[:,3] + 1j*dat[:,4]
        self.plot_kr_of_r()

    

    def calculate_dispersion_relation_and_plot(self, mode):

        # %matplotlib auto
        plt.rcParams.update({"font.size": 22})
        plt.rcParams.update({"lines.markersize": 10})

        self.load_all_profs()
        self.calc_all_derivs()
        self.calc_equilibrium_fields()
        self.find_res_surface()
        
        idx = int(self.prof_length * self.r_per)
        idx_range = np.linspace(self.r_range_start, self.prof_length - 1, self.n_points)

        res = list()
        omega_used = list()
        if self.options['omegaOfk']:
            for w in self.options['omega_range']:
                equation_k = lambda k: self.createDispersionEquation(k, w, int(idx), mode=mode)
                contour = cx.Circle(0, self.options['radius'])
                if self.options['der']:
                    roots = contour.roots(
                        equation_k,
                        df=lambda k: complex(grad(equation_k, holomorphic=True)(k)),
                        guess_roots_symmetry=lambda z: [-z],
                    )
                else:
                    roots = contour.roots(
                        equation_k,
                        guess_roots_symmetry=lambda z: [-z],
                        int_method=self.options['int_method'],
                    )
                res.append(roots)
                omega_used.append(w)
            k_w1 = list()
            k_w2 = list()
            omega_found = list()
            for k in range(len(res)):
                try:
                    if res[k].roots[0].real > 0:
                        k_w1.append(res[k].roots[0])
                        k_w2.append(res[k].roots[1])
                        omega_found.append(omega_used[k])
                    else:
                        k_w1.append(res[k].roots[1])
                        k_w2.append(res[k].roots[0])
                        omega_found.append(omega_used[k])
                except:
                    pass
            k_w1 = np.array(k_w1)
            k_w2 = np.array(k_w2)

            if self.options['save']:
                if self.options['der']:
                    self.save2txt(
                        omega_found,
                        k_w1,
                        k_w2,
                        self.mode + "_" + self.options['prof']+ "_" + "w(k)" + "_jax" + ".txt",
                    )
                else:
                    self.save2txt(
                        omega_found, k_w1, k_w2, self.mode + "_" + self.options['prof']+ "_" + "w(k)" + ".txt"
                    )
            plt.figure()
            plt.grid()
            plt.plot(k_w1.real, omega_used, ".k")
            plt.plot(k_w2.real, omega_used, ".k")
            plt.title("Real")
            plt.xlabel("Re(k)")
            plt.ylabel("$\omega$")

            plt.figure()
            plt.grid()
            plt.plot(k_w1.imag, omega_used, ".k")
            plt.plot(k_w2.imag, omega_used, ".k")
            plt.title("Imaginary")
            plt.xlabel("Im(k)")
            plt.ylabel("$\omega$")

        res = list()
        r_used = list()

        contour=cx.Rectangle([-10,10],[-10,10])

        if self.options['kOfr']:
            for r in idx_range:
                equation_k = lambda k: self.createDispersionEquation(k, self.options['omega'], int(r), mode=self.options['mode'])
                #equation_k = lambda k: createDispersionEquationKIM(k, omega, r_prof[int(r)])
                iterations=0
                
                try:
                    roots_number=contour.count_roots(equation_k)
                except:
                    print("Error in count_roots")
                    continue
                while roots_number!=2:
                    if roots_number>4:
                        contour=cx.Rectangle(np.array(contour.x_range)/2,np.array(contour.y_range)/2)
                    elif roots_number>2:
                        if np.max(contour.x_range)<2:
                            break
                        contour=cx.Rectangle(np.array(contour.x_range)*3/4,np.array(contour.y_range)*3/4)
                    elif roots_number<2:
                        contour=cx.Rectangle(np.array(contour.x_range)*2,np.array(contour.y_range)*2)
                    roots_number=contour.count_roots(equation_k)
                    iterations+=1
                    if iterations>100:
                        break
                if self.options['der']:
                    roots = contour.roots(
                        equation_k,
                        df=lambda k: complex(grad(equation_k, holomorphic=True)(k)),
                        guess_roots_symmetry=lambda z: [-z],
                        verbose=True
                    )
                else:
                    roots = contour.roots(
                        equation_k,
                        guess_roots_symmetry=lambda z: [-z],
                        int_method=self.options['int_method'],
                        verbose=True,
                        int_abs_tol=0.1,
                        root_err_tol=1e-3
                    )
                res.append(roots)
                r_used.append(self.general_dat['r'][int(r)])
                
            k_r1 = list()
            k_r2 = list()
            r_found = list()
            for k in range(len(res)):
                if len(res[k].roots)==2:
                    if res[k].roots[0].real > 0:
                        k_r1.append(res[k].roots[0])
                        k_r2.append(res[k].roots[1])
                        r_found.append(r_used[k])
                    else:
                        k_r1.append(res[k].roots[1])
                        k_r2.append(res[k].roots[0])
                        r_found.append(r_used[k])
            k_r1 = np.array(k_r1)
            k_r2 = np.array(k_r2)

            if self.options['save']:
                if self.options['der']:
                    self.save2txt(
                        r_found,
                        k_r1,
                        k_r2,
                        self.options['mode'] + "_" + self.options['prof'] + "_" + "k(r)" + "_" + "_jax" + ".txt",
                    )
                else:
                    self.save2txt(r_found, k_r1, k_r2, self.options['mode'] + "_" + self.options['prof']+ "_" + "k(r)" + ".txt")

            plt.figure()
            plt.grid()
            plt.plot(r_found, k_r1.real, ".k")
            plt.plot(r_found, k_r2.real, ".k")
            plt.vlines(
                self.r_prof[self.idx_res],
                np.min(np.append(k_r1.real, k_r2.real)),
                np.max(np.append(k_r1.real, k_r2.real)),
                colors="r",
                linestyles="dashed",
                label="Resonant surface",
            )
            plt.title("Real")
            plt.xlabel("r")
            plt.ylabel("Re(k)")

            plt.figure()
            plt.grid()
            plt.plot(r_found, k_r1.imag, ".k")
            plt.plot(r_found, k_r2.imag, ".k")
            plt.vlines(
                self.r_prof[self.idx_res],
                np.min(np.append(k_r1.imag, k_r2.imag)),
                np.max(np.append(k_r1.imag, k_r2.imag)),
                colors="r",
                linestyles="dashed",
                label="Resonant surface",
            )
            plt.title("Imaginary")
            plt.xlabel("r")
            plt.ylabel("Im(k)")
            plt.show()

    def print_species_data(self, spec):
        print("")
        print(f"Species: {spec}")
        [print(f"{key}: {value}") for key, value in self.spec_dat[spec].items()]
    
    def plot_ks(self):
        self.load_all_profs()
        self.calc_parameters()
        plt.figure()
        plt.plot(self.general_dat['r'], self.general_dat['kp'])
        plt.grid()
        plt.show()

    def generate_non_uniform_grid(self):
        gobj = WKB_Grid(self.general_dat['r'], self.res_surf_val)
        print(self.res_surf_val)
        #self.r_new = gobj.generate_non_equidistant_grid(gobj.r_min, gobj.r_max, 30, gobj.r_res, 3)
        self.r_new = np.array(gobj.gen_grid(self.general_dat['r'], gobj.r_res, 5.0, 50))

    def interpolate_on_non_uniform_grid(self):
        self.general_dat['old_r'] = np.copy(self.general_dat['r'])
        self.general_dat['r'] = self.r_new

        for key in self.general_dat:
            try:
                if len(self.general_dat[key]) == self.general_dat['prof_length']:
                    if key == 'old_r':
                        continue
                    self.general_dat[key] = np.interp(self.r_new, self.general_dat['old_r'], self.general_dat[key])
            except:
                pass
        for spec in self.species:
            for key in self.spec_dat[spec]:
                if key =='n' or key =='T':
                    self.spec_dat[spec][key] = np.interp(self.general_dat['r'], self.general_dat['old_r'], self.spec_dat[spec][key])
        self.general_dat['old_prof_length'] = self.general_dat['prof_length']  
        self.general_dat['prof_length'] = len(self.r_new)
        
        print('New grid length is ', self.general_dat['prof_length'])
            
    def test_non_uniform_grid(self):
        self.load_all_profs()
        self.calc_parameters()
        self.find_res_surface()
        self.generate_non_uniform_grid()    
        self.interpolate_on_non_uniform_grid()
        plt.figure()
        plt.plot(self.r_new, np.zeros_like(self.r_new), marker='.')
        plt.axvline(self.res_surf_val, color='r')
        plt.grid()
        plt.show()
        
    def write_dispersion_relation_to_h5(self, file):
        self.disp_rel_to_h5(file, '', 'w')
    
    def append_dispersion_relation_to_h5(self, file, study='D'):
        self.disp_rel_to_h5(file, study, 'a')
        
    def disp_rel_to_h5(self, file, study, mode):
        h5f = h5py.File(file, mode)
        try:
            h5f.create_dataset(study + '/r_found', data=self.r_found)
        except:
            raise ValueError(f"Study {study} already exists")
        h5f.create_dataset(study + '/k_r1', data=self.k_r1)
        h5f.create_dataset(study + '/k_r2', data=self.k_r2)
        for option in self.options.keys():
            h5f.create_dataset(study + '/options/'+ option, data=self.options[option])
        h5f.create_dataset(study + '/r_res', data=self.res_surf_val)
        h5f.close()
    
    def species_data_to_h5(self, file, study, mode):
        h5f = h5py.File(file, mode)
        for spec in self.species:
            print(f"Species {spec}")
            for key in self.spec_dat[spec]:
                if key == 'I00' or key == 'I20':
                    for mphi in range(0, self.options['max_cyclotron_harmonic']+1):
                        #try:
                        h5f.create_dataset(study + '/'+ spec + '/'+ key + '/mphi_' + str(mphi), data=self.spec_dat[spec][key][mphi])
                        #except:
                            #raise ValueError(f"Species {spec} problem with key {key}")
                else:
                    try:
                        h5f.create_dataset(study + '/'+ spec + '/'+ key, data=self.spec_dat[spec][key])
                    except:
                        raise ValueError(f"Species {spec} with key {key} already exists")
        h5f.close()
    
    def general_data_to_h5(self, file, study, mode):
        h5f = h5py.File(file, mode)
        for key in self.general_dat:
            try:
                h5f.create_dataset(study + '/general/'+ key, data=self.general_dat[key])
            except:
                raise ValueError(f"Key {key} already exists")
        h5f.close()
        
    def write_all_data_to_h5(self, file, study='', mode='a'):
        if os.path.isfile(file) and mode == 'w':
            os.remove(file)
        elif not os.path.isfile(file) and mode == 'a':
            mode = 'w'
        mode2 = 'a'
        
        self.disp_rel_to_h5(file, study, mode)
        self.species_data_to_h5(file, study, mode2)
        self.general_data_to_h5(file, study, mode2)
    
    def set_collision_mode(self, collisions):
        assert collisions in self.possible_coll_models, "Collision model not supported"
        self.options['Collisions'] = collisions
    
    def set_model_mode(self, mode):
        assert mode in self.possible_operation_modes, "Mode not supported"
        self.options['mode'] = mode
    

def determine_dispersion_for_all_species():
    specs = {0: ['e', 'H'], 1: ['e', 'D']}
    specs = {0: ['e', 'D']}
    spec_mass = {0: [e_mass, p_mass], 1: [e_mass, 2*p_mass]}
    spec_mass = {0: [e_mass, 2*p_mass]}
    spec_charge_num = {0: [-1,1]}
    modes = {0: 'KIM', 1: 'horton'}

    for i, spec in specs.items():
        print(spec)
        for j, mode in modes.items():
            kwkb = KIM_WKB(species=spec, spec_mass=spec_mass[i], spec_charge_num=spec_charge_num[i])
            kwkb.options['n_points'] = 50
            kwkb.prof_path = '../../../kim-wkb/profiles_parab/'

            print('Mode is ', mode, ', Collisions: ', not kwkb.options['Collisions'])
            kwkb.calc_dispersion_relation_k_of_r(mode=mode)
            #kwkb.plot_kr_of_r()
            np.savetxt(f'./{mode}_{spec[1]}.dat',np.vstack((np.real(kwkb.r_found),np.real(kwkb.k_r1), np.imag(kwkb.k_r1), np.real(kwkb.k_r2), np.imag(kwkb.k_r2))).T)


def test_FokkerPlanck():
    
    specs = {0: ['e', 'D']}
    spec_mass = {0: [e_mass, 2*p_mass]}
    spec_charge_num = {0: [-1,1]}
    mode = 'KIM'

    kwkb = KIM_WKB(species=specs[0], spec_mass=spec_mass[0], spec_charge_num=spec_charge_num[0])
    kwkb.contour_limit = 10 # 50 works for H, 20 for D
    kwkb.prof_path = '../../../kim-wkb/profiles_parab/'
    #kwkb.plot_ks()
    #exit()
    mphi_max = 0
    kwkb.options['n_points'] = 50
    kwkb.options['number_of_roots_to_find'] = 8
    kwkb.set_collision_mode('FokkerPlanck')
    kwkb.options['max_cyclotron_harmonic'] = mphi_max
    kwkb.options['der'] = False
    kwkb.options['log'] = False
    kwkb.options['new_grid'] = True
    #kwkb.test_non_uniform_grid()
    #kwkb.plot_x2()
    #kwkb.test_bessel_of_mphi()
    #exit()
    #kwkb.plot_dispersion_equation_KIM_FokkerPlanck()
    print('Mode is ', mode, ', Collisions: ', kwkb.options['Collisions'])
    kwkb.set_output_h5_file(f'./{mode}_{kwkb.options["Collisions"]}.h5', append_or_write = 'w')
    kwkb.calc_dispersion_relation_k_of_r(mode=mode)
    kwkb.write_all_data_to_h5(f'./{mode}_{kwkb.options["Collisions"]}.h5', mode = 'a')

#    np.savetxt(,np.vstack((np.real(kwkb.r_found),np.real(kwkb.k_r1), np.imag(kwkb.k_r1), np.real(kwkb.k_r2), np.imag(kwkb.k_r2))).T)
    
if __name__ == "__main__":

    if len(sys.argv) > 1:
        assert sys.argv[1] in ['horton', 'KIM'], "Mode not supported"
        mode = sys.argv[1]
    else:
        mode = 'horton'


    test_FokkerPlanck()
