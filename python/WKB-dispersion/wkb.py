import cxroots as cx
import numpy as np
import matplotlib.pyplot as plt

from plasmapy.dispersion import plasma_dispersion_func as plasma_disp
from wkb_grid import WKB_Grid

import os
import sys
from scipy.integrate import solve_ivp
from scipy.special import iv as Bessel
from jax import grad
import jax
import jax.numpy as jnp
import h5py
import logging

from Bessel_calculation import calc_needed_bessel_of_mphi
from KIMDispersion_Horton import KIMDispersion_Horton
from KIMDispersion_Krook import KIMDispersion_Krook
from KIMDispersion_FokkerPlanck import KIMDispersion_FokkerPlanck
from DispersionEquationFactory import *

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../KiLCA-QB/python/susc_functions/'))
import susc_funcs

sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
from constants import *


class KIM_WKB():
    """
        Class to calculate dispersion relation for KIM integral kernels in the WKB 
        approximation. Uses the cxroots module to find the roots of the dispersion
        equation.
        
        Modes available:
        - Horton: for collisionless FLRE model from Horton 1999, Rev. of Mod. Phys.
        - KIM: kinetic integral model based on KiLCA. Comes in 3 different collision models:
            - Krook: simple Krook collision model
            - collisionless: Krook, but with collision frequencies set to zero
            - Fokker_Planck: FP-type collision model with an Ornstein-Uhlenebeck collision
                model including an energy-preserving integral term. This collision model
                results in the susceptibility functions (I00 and I20) imported from susc_funcs.
                The needed Bessel functions are calculated in the Bessel_calculation file.
        The modes are implemented with abstract base classes (ABC).
        
        Data is handled with dictionaries. Different dicts are: 
            - options: contains options for the run including possible collision models
            - species_dat: contains all the species specific profiles like density, 
                temperature, as well as derived ones like the thermodynamic forces.
            - general_dat: contains data that is not specific to a species, e.g. 
                radial electric field, radial grid,...
            - equil_dat: contains data for the magnetic equilibrium of the cylinder that
                is calculated in calc_equilibrium
                
       The output data is written to the hdf5 file format using the h5py module. 
    """
    
    btor = -17977.413  # AUG config
    R0 = 165.0  # major radius in cm
    r_plas = 64.7  # plasma radius

    m_mode = 6  # poloidal mode number
    n_mode = 2  # toroidal mode number

    prof_path = ''

    bessel_large_arg_limit = 5
    contour_limit = 10
    
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
               'new_grid': True,
               'number_of_roots_to_find': 4,
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

    spec_keys = ['n', 'T', 'vT', 'nu', 'mass', 'Ai', 'charge', 'Zi', 'I00', 'I20']
    spec_dat_units = {'n': 'cm', 'T': 'erg', 'vT': 'cm/s', 'nu': '1/s', 'mass': 'g', 'Ai': '1', 'charge': 'statC', 'Zi': '1', 'I00': '1', 'I20': '1'}
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
        
    def set_output_h5_file(self, file, append_or_write='w'):
        self.h5f = file
        self.h5_append_or_write = append_or_write
        
    def calc_dispersion_relation_k_of_r(self, mode, collisions):  

        self.set_model_mode(mode)
        self.set_collision_mode(collisions)
        self.initialize_data()
        
        self.init_dispersion_equation()

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

    def set_model_mode(self, mode):
        assert mode in self.possible_operation_modes, "Mode not supported"
        self.options['mode'] = mode
        
    def set_collision_mode(self, collisions):
        if self.options['mode'] != 'horton':
            assert collisions in self.possible_coll_models, "Collision model not supported"
        self.options['Collisions'] = collisions

    def initialize_data(self):
        self.load_all_profs()
        self.find_res_surface()
        if self.options['new_grid']:
            self.generate_non_uniform_grid()
            self.interpolate_on_non_uniform_grid()
        self.calc_parameters()

        if self.options['Collisions'] == 'FokkerPlanck':
            self.calc_needed_susc_funcs()

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

    def find_res_surface(self):
        # Find resonant surface
        self.idx_res = self.val2ind(-self.m_mode / self.n_mode, self.general_dat['q'])
        self.res_surf_val = np.interp(self.m_mode / self.n_mode, np.abs(self.general_dat['q']), self.general_dat['r'])
 
    def val2ind(self, value, array):
        idx = np.argmin(np.abs(array - value))
        return idx   

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

    def calc_parameters(self):
        self.calc_equilibrium()
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

    def calc_equilibrium(self):
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

    def init_dispersion_equation(self):
        if self.options['mode'] == 'horton':
            self.dispersion_model = KIMDispersion_Horton_Factory.get_dispersion_model()
        if self.options['mode'] == 'KIM':
            if self.options['Collisions'] == 'Krook':
                self.dispersion_model = KIMDispersion_Krook_Factory.get_dispersion_model()
            if self.options['Collisions'] == 'FokkerPlanck':
                self.dispersion_model = KIMDispersion_FokkerPlanck_Factory.get_dispersion_model()
        self.dispersion_model.initialize(self.options, self.species, self.spec_dat, self.general_dat, self.equil_dat)
        
    def find_roots(self, r_ind, contour):
        equation_k = lambda k: self.calc_dispersion_equation_for_kr_values_at_r(k, r_ind, mode=self.options['mode'])
        iterations=0
        roots_number = 0
        search_scale_bigger = 2.0
        search_scale_smaller = 1.5
        iteration = 0
        try:
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
            except:
                print(f"Error in roots for r = {r_ind}")
                return 1
        return roots
    
    def calc_dispersion_equation_for_kr_values_at_r(self, kr, position):
        if type(position) == int:
            r_eval = self.general_dat['r'][position]
        if type(kr) == complex or type(kr) == np.complex128:
            dispersion_equation = self.dispersion_model.dispersion_equation(kr, position)
        else:
            raise ValueError("Other types of kr not implemented yet")
        return dispersion_equation

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

    def plot_non_uniform_grid(self):
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

    def plot_bessel_of_mphi(self):
        mphi_list = [0]
        eval_b = np.linspace(-7,7,100)
        plt.figure(figsize=(8,6))
        Bess0 = np.zeros(len(eval_b), dtype=complex)
        Bess1 = np.zeros(len(eval_b), dtype=complex)
        self.bessel_large_arg_limit = 1
        for i, mphi in enumerate(mphi_list):
            for j, b in enumerate(eval_b):
                Bess0[j], Bess1[j] = calc_needed_bessel_of_mphi(mphi, b)
            plt.plot(eval_b, np.real(Bess0), label=f'Re, Bessel0, mphi = {mphi}', marker='x')
            plt.plot(eval_b, np.real(Bess1), label=f'Re, Bessel1, mphi = {mphi}', marker='x')
            plt.plot(eval_b, np.imag(Bess0), label=f'Im, Bessel0, mphi = {mphi}', ls=':', marker='.')
            plt.plot(eval_b, np.imag(Bess1), label=f'Im, Bessel1, mphi = {mphi}', ls=':', marker='.')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.grid()
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

    def plot_susc_funcs(self):
        self.load_all_profs()
        self.calc_parameters()
        self.calc_needed_susc_funcs()
        
        plt.figure()
        for spec in self.species: 
            if spec == 'e':
                continue
            for m_phi in range(0, self.options['max_cyclotron_harmonic']+1):
                plt.plot(self.general_dat['r'], self.spec_dat[spec]['I00'][m_phi].real, label=spec)
        plt.legend()
        plt.show()
        
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

    def plot_dispersion_equation_KIM_FokkerPlanck(self):
        self.initialize_data()

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


def test_FokkerPlanck():
    specs = {0: ['e', 'D']}
    spec_mass = {0: [e_mass, 2*p_mass]}
    spec_charge_num = {0: [-1,1]}
    mode = 'KIM'

    kwkb = KIM_WKB(species=specs[0], spec_mass=spec_mass[0], spec_charge_num=spec_charge_num[0])
    kwkb.contour_limit = 10 # 50 works for H, 20 for D
    kwkb.prof_path = '../../../kim-wkb/profiles_parab/'
    mphi_max = 0
    
    kwkb.options['n_points'] = 50
    kwkb.options['number_of_roots_to_find'] = 8
    kwkb.set_collision_mode('Krook')
    kwkb.options['max_cyclotron_harmonic'] = mphi_max
    kwkb.options['der'] = False
    kwkb.options['log'] = False
    kwkb.set_output_h5_file(f'./{mode}_{kwkb.options["Collisions"]}.h5', append_or_write = 'w')
    kwkb.calc_dispersion_relation_k_of_r(mode=mode)
    kwkb.write_all_data_to_h5(f'./{mode}_{kwkb.options["Collisions"]}.h5', mode = 'a')

def test_ABC():
    specs = {0: ['e', 'D']}
    spec_mass = {0: [e_mass, 2*p_mass]}
    spec_charge_num = {0: [-1,1]}
    mode = 'KIM'
    collisions = 'Krook'

    kwkb = KIM_WKB(species=specs[0], spec_mass=spec_mass[0], spec_charge_num=spec_charge_num[0])
    kwkb.contour_limit = 10 # 50 works for H, 20 for D
    kwkb.prof_path = '../../../kim-wkb/profiles_parab/'
    kwkb.options['n_points'] = 50
    kwkb.options['number_of_roots_to_find'] = 8
    kwkb.options['der'] = False # set True if jax is used for differentiation
    kwkb.options['log'] = False
    kwkb.options['new_grid'] = True

    study = 'ABC_test'
    kwkb.set_output_h5_file(f'./'+study+f'_{mode}_{kwkb.options["Collisions"]}.h5', append_or_write = 'w')
    kwkb.calc_dispersion_relation_k_of_r(mode=mode, collisions=collisions)
    kwkb.write_all_data_to_h5(f'./'+study+f'_{mode}_{kwkb.options["Collisions"]}.h5', mode = 'a')

if __name__ == "__main__":

    if len(sys.argv) > 1:
        assert sys.argv[1] in ['horton', 'KIM'], "Mode not supported"
        mode = sys.argv[1]
    else:
        mode = 'horton'

    test_ABC()