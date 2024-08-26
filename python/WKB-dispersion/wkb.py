# %%
import cxroots as cx
import numpy as np
import matplotlib.pyplot as plt

from plasmapy.dispersion import plasma_dispersion_func as plasma_disp

import os
import sys
from scipy.integrate import solve_ivp
from scipy.special import iv as Bessel
from jax import grad
import jax
import jax.numpy as jnp
import h5py


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
    contour_limit = 10
    
    ##Options
    options = {'prof': 'parab',
               'n_points': 50,
               'omega': 0,
               'radius': 250,
               'r_per': 0.9,
               'r_range_start': 10,
               'omegaOfk': False,
               'kOfr': True,
               'mode': 'KIM2',
               'int_method': 'quad',
               'der': False,
               'save': True,
               'noCollisions': True}
    options['omega_range'] = np.linspace(0, 50, options['n_points'])

    possible_operation_modes = ['KIM2', 'horton']

    # which profile should be loaded: const, parab or full
    prof = "parab"
    # number of points to evaluate
    n_points = 50
    # Define perturbation frequency
    omega = 0
    omega_range = np.linspace(0, 50, n_points)
    # Define radial position
    r_per = 0.9  # percentage of total
    r_range_start = 10  # start index
    # Define k-radius to search over
    radius = 250
    # Calculate w(k)
    omegaOfk = False
    # Calculate k(r)
    kOfr = True
    # Choose model: KIM2 or horton
    mode = "KIM2"
    # Choose method for integration: quad or romb
    int_method = "romb"
    # Should the derivative be approximated with jax.grad
    der = False
    # Should the results be saved
    save = True
    noCollisions=True

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
        #f = lambda r, u: dudr(r, u, r_prof, q_prof, dpress_prof)
        f = lambda r,u: -2*r*u/(self.R0**2 * self.general_dat['q'][self.val2ind(r, self.general_dat['r'])]**2 + r**2)-8*pi*dpress_prof[self.val2ind(r, self.general_dat['r'])]
        u0 = self.btor**2 * (1 + self.general_dat['r'][0] ** 2 / (self.R0**2 * self.general_dat['q'][0] ** 2))
        u = solve_ivp(f, t_span=[self.general_dat['r'][0], self.general_dat['r'][-1]], y0=[u0], t_eval=self.general_dat['r'])
        u = u.y[0]
        self.equil_dat['B0z'] = np.sign(self.btor) * np.sqrt(u / (1 + self.general_dat['r']**2 / (self.R0**2 * self.general_dat['q']**2)))
        self.equil_dat['B0th'] = self.equil_dat['B0z'] * self.general_dat['r'] / (self.general_dat['q'] * self.R0)
        self.equil_dat['B0'] = np.sqrt(self.equil_dat['B0th']**2 + self.equil_dat['B0z']**2)
        self.equil_dat['hz'] = self.equil_dat['B0z'] / self.equil_dat['B0']
        self.equil_dat['hth'] = self.equil_dat['B0th'] / self.equil_dat['B0']


    def dudr(self, r, u, r_prof, q_prof, dpress_prof):
        nlagr = 4
        nder = 0
        # Find index which satisfies: p(i-1)<xi<p(i)
        ir = (np.abs(r_prof - r)).argmin()
        if r_prof[ir] < r:
            ir += 1
        if (ir - nlagr / 2) > 0:
            ibeg = ir - nlagr / 2
        else:
            ibeg = 0
        iend = ibeg + nlagr
        if iend > prof_length:
            iend = prof_length
            ibeg = iend - nlagr
        ibeg = int(ibeg)
        iend = int(iend)
        coef = self.plag_coeff(nlagr, nder, r, r_prof[ibeg:iend])
        q = np.sum(coef[0, :] * q_prof[ibeg:iend])
        dpress = np.sum(coef[0, :] * dpress_prof[ibeg:iend])
        return -2 * r * u / (q**2 * self.R0**2 + r**2) - 8 * pi * dpress


    def plag_coeff(self, npoi, nder, x, xp):
        coef = np.zeros((1, npoi))
        for k in range(npoi):
            coef[0, k] = 1
            for m in range(npoi):
                if k == m:
                    continue
                coef[0, k] = coef[0, k] * (x - xp[m]) / (xp[k] - xp[m])
        # Maybe add more if nder !=0 later
        return coef

    
    def calc_parameters(self):
        self.calcEquilibrium()
        if self.options['noCollisions']:
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
            self.spec_dat[spec]['omega_c'] = self.spec_dat[spec]['charge'] * self.equil_dat['B0'] / (self.spec_dat[spec]['mass'] * sol)
            self.spec_dat[spec]['lambda_D'] = np.sqrt(self.spec_dat[spec]['T'] * ev / (4.0 * np.pi * self.spec_dat[spec]['n'] * self.spec_dat[spec]['charge']**2))
            self.spec_dat[spec]['A1'] = self.spec_dat[spec]['dndr'] / self.spec_dat[spec]['n'] \
                - self.spec_dat[spec]['charge'] / (self.spec_dat[spec]['T'] * ev) * self.general_dat['Er'] \
                - 3 / (2 * self.spec_dat[spec]['T']) * self.spec_dat[spec]['dTdr']
            self.spec_dat[spec]['A2'] = self.spec_dat[spec]['dTdr'] / self.spec_dat[spec]['T']


            self.spec_dat[spec]['z0'] = -(self.general_dat['om_E'] - self.options['omega'] - 1j * self.spec_dat[spec]['nu']) \
                / (self.general_dat['kp'] * np.sqrt(2) * self.spec_dat[spec]['vT'])
            self.spec_dat[spec]['rho_TL'] = self.spec_dat[spec]['vT'] / self.spec_dat[spec]['omega_c']

        
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


    def makePlot(self, y):
        plt.figure()
        plt.grid()
        plt.plot(self.r_prof, y)
        return


    def save4compare(self, r1, r2):
        out1 = list()
        out2 = list()
        for k in range(len(r1)):
            out1.append(r1[k])
            out2.append(r2[k])
        return np.array(out1), np.array(out2)


    def save2txt(self, r_or_omega, k1, k2, name):
        if self.kOfr:
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


    def plotHortonKIM(self, rh, rk, h1, h2, k1, k2, prof, omega=0):
        if omegaOfk:
            fig1, ax1 = plt.subplots()
            ax2 = ax1.twiny()
            plt.grid()
            ax1.plot(h1.real, omega, "k.", label="Horton")
            ax1.plot(h2.real, omega, "k.")
            ax2.plot(k1.real, omega, "rx", label="KIM")
            ax2.plot(k2.real, omega, "rx")
            ax1.set_ylabel("$\omega$ [$s^{-1}$]")
            ax1.set_xlabel("k Horton [$cm^{-1}$]")
            ax2.set_xlabel("k KIM [$cm^{-1}$]")
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines2 + lines, labels2 + labels, loc=2)
            plt.title("Real: " + prof)

            fig2, ax3 = plt.subplots()
            ax4 = ax3.twiny()
            plt.grid()
            ax3.plot(h1.imag, omega, "k.", label="Horton")
            ax3.plot(h2.imag, omega, "k.")
            ax4.plot(k1.imag, omega, "rx", label="KIM")
            ax4.plot(k2.imag, omega, "rx")
            ax3.set_ylabel("$\omega$ [$s^{-1}$]")
            ax3.set_xlabel("k horton [$cm^{-1}$]")
            ax4.set_xlabel("k KIM [$cm^{-1}$]")
            lines, labels = ax3.get_legend_handles_labels()
            lines2, labels2 = ax4.get_legend_handles_labels()
            ax3.legend(lines2 + lines, labels2 + labels, loc=2)
            plt.title("Imaginary: " + prof)
        elif kOfr:
            fig1, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            plt.grid()
            ax1.plot(rh, h1.real, "k.", label="Horton")
            ax1.plot(rh, h2.real, "k.")
            ax2.plot(rk, k1.real, "rx", label="KIM")
            ax2.plot(rk, k2.real, "rx")
            ax1.vlines(
                self.r_prof[self.idx_res],
                np.min(np.append(h1.real, h2.real)),
                np.max(np.append(h1.real, h2.real)),
                colors="b",
                linestyles="dashed",
                label="Resonant surface",
            )
            ax1.set_ylabel("k Horton [$cm^{-1}$]")
            ax2.set_ylabel("k KIM [$cm^{-1}$]")
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines2 + lines, labels2 + labels, loc=2)
            ax1.set_xlabel("r [$cm$]")
            plt.title("Real: " + prof)

            fig2, ax3 = plt.subplots()
            ax4 = ax3.twinx()
            plt.grid()
            ax3.plot(rh, h1.imag, "k.", label="Horton")
            ax3.plot(rh, h2.imag, "k.")
            ax4.plot(rk, k1.imag, "rx", label="KIM")
            ax4.plot(rk, k2.imag, "rx")
            ax3.vlines(
                self.r_prof[self.idx_res],
                np.min(np.append(h1.imag, h2.imag)),
                np.max(np.append(h1.imag, h2.imag)),
                colors="b",
                linestyles="dashed",
                label="Resonant surface",
            )
            ax3.set_ylabel("k Horton [$cm^{-1}$]")
            ax4.set_ylabel("k KIM [$cm^{-1}$]")
            lines, labels = ax3.get_legend_handles_labels()
            lines2, labels2 = ax4.get_legend_handles_labels()
            ax3.legend(lines2 + lines, labels2 + labels, loc=2)
            ax3.set_xlabel("r [$cm$]")
            plt.title("Imaginary: " + prof)
        return

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
        for spec in self.species:
            if spec != 'e':
                self.ni_prof = self.n_prof
                self.spec_dat[spec]['n'] = self.ni_prof

    
    def calc_all_derivs(self):
        # Calculate derivatives
        self.calc_species_derivs()
        self.calc_general_derivs()
        #dndr, dTedr, dqdr, dnidr, dTidr = calcDerivatives()
        self.dndr, self.dTedr, self.dqdr, self.dnidr, self.dTidr = list(map(lambda x: np.gradient(x,self.r_prof),[self.n_prof,self.Te_prof,self.q_prof,self.ni_prof,self.Ti_prof]))
    
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


    def calc_dispersion_equation_for_kr_values_at_r(self, kr, position, mode="KIM"):
        if type(position) == int:
            r_eval = self.general_dat['r'][position]
        if type(kr) == complex or type(kr) == np.complex128:
            #disp_prof = self.create_dispersion_equation_profile(kr, mode)
            dispersion_equation = self.create_dispersion_equation_single(kr, position, mode)
            #dispersion_equation = np.interp(r_eval, self.general_dat['r'], disp_prof)
        else:
            dispersion_equation = np.zeros(len(kr))
            for i, val in enumerate(kr):
                disp_prof = self.create_dispersion_equation_profile(val, mode)
                dispersion_equation[i] = np.interp(r_eval, self.r_prof, disp_prof)
        return dispersion_equation

    def create_dispersion_equation_single(self, kr, r_indx, mode="KIM"):
        assert mode in self.possible_operation_modes, f"Mode {mode} not supported"

        self.general_dat['kperp'] = np.sqrt(self.general_dat['ks']**2 + kr**2)
        if mode == 'KIM2':
            dispersion_equation = self.calc_dispersion_equation_KIM_of_r_index(kr, r_indx)
        elif mode == "horton":
            raise ValueError("Horton mode not supported for single calculation")
            #dispersion_equation = self.calc_dispersion_equation_horton_profile(kr)
        #if np.isnan(dispersion_equation).any():
        #    print(f'dispersion_equation contains NaNs')
        return dispersion_equation

    def create_dispersion_equation_profile(self, kr, mode="KIM"):
        
        assert mode in self.possible_operation_modes, f"Mode {mode} not supported"
        
        self.general_dat['kperp'] = np.sqrt(self.general_dat['ks']**2 + kr**2)
        if mode == 'KIM2':
            dispersion_equation = self.calc_dispersion_equation_KIM_profile(kr)
        elif mode == "horton":
            dispersion_equation = self.calc_dispersion_equation_horton_profile(kr)
        #if np.isnan(dispersion_equation).any():
        #    print(f'dispersion_equation contains NaNs')
        return dispersion_equation

    def calc_dispersion_equation_KIM_profile(self, kr):
        dispersion_equation = kr**2 + self.general_dat['ks']**2 + self.general_dat['kp']**2
        for spec in self.species:
            self.spec_dat[spec]['eval_b'] = self.general_dat['kperp']**2 * self.spec_dat[spec]['rho_TL']**2
            BesselProd0, BesselProd1 = self.calc_needed_bessel(spec, self.spec_dat[spec]['eval_b'])

            dispersion_equation += self.spec_dat[spec]['lambda_D']**-2 * \
                (
                    1.0 - self.general_dat['ks'] * self.spec_dat[spec]['rho_TL'] \
                / (self.general_dat['kp'] * np.sqrt(2)) * (self.spec_dat[spec]['A1'] * BesselProd0 * plasma_disp(self.spec_dat[spec]['z0']) + \
                self.spec_dat[spec]['A2'] * (plasma_disp(self.spec_dat[spec]['z0']) * (1 + self.spec_dat[spec]['eval_b'] + self.spec_dat[spec]['z0']**2)\
                * np.exp(-self.spec_dat[spec]['eval_b']) + BesselProd1 * self.spec_dat[spec]['eval_b'] + self.spec_dat[spec]['z0'] * BesselProd0))
                )
        #print(dispersion_equation)
        return dispersion_equation
    
    def calc_needed_bessel(self, spec, eval_b):
        BesselProd0 = np.zeros(len(self.general_dat['r']), dtype=complex)
        BesselProd1 = np.zeros(len(self.general_dat['r']), dtype=complex)
        for i, eval_b in enumerate(self.spec_dat[spec]['eval_b']):
            if np.abs(eval_b)>self.bessel_large_arg_limit:
                BesselProd0[i] = 1.0 / (np.sqrt(2*np.pi*eval_b))
                BesselProd1[i] = 1.0 / (np.sqrt(2*np.pi*eval_b)*(1 + 1 / eval_b**2)**(1/4)) * np.exp(np.arcsinh(-1 / eval_b) + eval_b \
                    * np.sqrt(1 + 1 / eval_b**2) - eval_b)
            else:
                BesselProd0[i] = Bessel(0, eval_b)*np.exp(-eval_b)
                BesselProd1[i] = Bessel(-1, eval_b)*np.exp(-eval_b)

        return BesselProd0, BesselProd1


    def calc_dispersion_equation_KIM_of_r_index(self, kr, r_indx):
        dispersion_equation = kr**2 + self.general_dat['ks'][r_indx]**2 + self.general_dat['kp'][r_indx]**2
        for spec in self.species:
            eval_b = self.general_dat['kperp'][r_indx]**2 * self.spec_dat[spec]['rho_TL'][r_indx]**2
            BesselProd0, BesselProd1 = self.calc_needed_bessel_single(spec, eval_b)

            dispersion_equation += self.spec_dat[spec]['lambda_D'][r_indx]**-2 * \
                (
                    1.0 - self.general_dat['ks'][r_indx] * self.spec_dat[spec]['rho_TL'][r_indx] \
                / (self.general_dat['kp'][r_indx] * np.sqrt(2)) * (self.spec_dat[spec]['A1'][r_indx] * BesselProd0 * plasma_disp(self.spec_dat[spec]['z0'][r_indx]) + \
                self.spec_dat[spec]['A2'][r_indx] * (plasma_disp(self.spec_dat[spec]['z0'][r_indx]) * (1 + eval_b + self.spec_dat[spec]['z0'][r_indx]**2)\
                * np.exp(-eval_b) + BesselProd1 * eval_b + self.spec_dat[spec]['z0'][r_indx] * BesselProd0))
                )
        #print(dispersion_equation)
        return dispersion_equation


    def calc_needed_bessel_single(self, spec, eval_b):
        if np.abs(eval_b)>self.bessel_large_arg_limit:
            BesselProd0 = 1.0 / (np.sqrt(2*np.pi*eval_b))
            BesselProd1 = 1.0 / (np.sqrt(2*np.pi*eval_b)*(1 + 1 / eval_b**2)**(1/4)) * np.exp(np.arcsinh(-1 / eval_b) + eval_b \
                * np.sqrt(1 + 1 / eval_b**2) - eval_b)
        else:
            BesselProd0 = Bessel(0, eval_b)*np.exp(-eval_b)
            BesselProd1 = Bessel(-1, eval_b)*np.exp(-eval_b)
        return BesselProd0, BesselProd1

    def calc_dispersion_equation_horton_profile(self, kr):
        ky = self.general_dat['ks']
        om_prime = self.options['omega'] - self.general_dat['om_E']
        
        dispersion_equation = 0 + 0j
        denom = 0 + 0j
        nom = 0 + 0j
        for spec in self.species:
            self.spec_dat[spec]['om_n'] = ky * sol / (self.spec_dat[spec]['charge'] * self.equil_dat['B0'] \
                * self.spec_dat[spec]['n']) * self.spec_dat[spec]['T'] * ev * self.spec_dat[spec]['dndr']
            self.spec_dat[spec]['om_T'] = ky * sol / (self.spec_dat[spec]['charge'] * self.equil_dat['B0']) * self.spec_dat[spec]['dTdr'] * ev

            self.spec_dat[spec]['z'] = om_prime / (np.sqrt(2) * self.general_dat['kp'] * self.spec_dat[spec]['vT'])
            self.spec_dat[spec]['W'] = (self.spec_dat[spec]['om_n'] - om_prime + self.spec_dat[spec]['om_T'] * (self.spec_dat[spec]['z']**2 + 0.5)) \
                * self.spec_dat[spec]['z'] / om_prime * plasma_disp(self.spec_dat[spec]['z']) + self.spec_dat[spec]['om_T'] * self.spec_dat[spec]['z']**2 / om_prime

            self.spec_dat[spec]['nom'] = self.spec_dat[spec]['lambda_D']**-2 * (self.spec_dat[spec]['W'] - 1.0 \
                -self.spec_dat[spec]['om_T'] * self.spec_dat[spec]['z'] / om_prime * plasma_disp(self.spec_dat[spec]['z']))
            self.spec_dat[spec]['denom'] = self.spec_dat[spec]['W'] * self.spec_dat[spec]['vT']**2 / (self.spec_dat[spec]['omega_c']**2 \
                * self.spec_dat[spec]['lambda_D']**2)

            nom += self.spec_dat[spec]['nom']
            denom += self.spec_dat[spec]['denom']

        nom -= self.general_dat['kp']**2
        dispersion_equation = nom / (1 + denom) - self.general_dat['kperp']**2
        
        return dispersion_equation

    # Create dispersion equation
    def createDispersionEquation(self, kr, omega, position, mode="KIM"):
        
        assert mode in self.possible_operation_modes, f"Mode {mode} not supported"
        
        # Do WKB
        (
            om_E_wkb,
            omce_wkb,
            omci_wkb,
            vTe_wkb,
            vTi_wkb,
            lambda_De_wkb,
            lambda_Di_wkb,
            nue_wkb,
            nui_wkb,
            kp_wkb,
            ks_wkb,
            A1_wkb,
            A2_wkb,
            A1i_wkb,
            A2i_wkb,
        ) = list(
            map(
                lambda x: np.interp(self.r_prof[position],self.r_prof,x),
                [self.om_E, self.omce, self.omci, self.vTe, self.vTi, self.lambda_De, self.lambda_Di, self.nue, self.nui, self.kp, self.ks, self.A1, self.A2, self.A1i, self.A2i],
            )
        )
        kperp_wkb = np.sqrt(ks_wkb**2 + kr**2)
        # Use either KIM or horton
        if mode == 'KIM2':
            z0e_wkb = -(om_E_wkb - omega - 1j * nue_wkb) / (kp_wkb * np.sqrt(2) * vTe_wkb)
            rho_TLe = vTe_wkb / omce_wkb
        
            # In WKB bp and bt are the same
            eval_b = kperp_wkb * rho_TLe**2
            if eval_b.real>100:
                BesselProd0=1/(np.sqrt(2*np.pi*eval_b))
                BesselProd1=1/(np.sqrt(2*np.pi*eval_b)*(1+1/eval_b**2)**(1/4))*np.exp(np.arcsinh(-1/eval_b)+eval_b*np.sqrt(1+1/eval_b**2)-eval_b)
            else:
                BesselProd0=Bessel(0,eval_b)*np.exp(-eval_b)
                BesselProd1=Bessel(-1,eval_b)*np.exp(-eval_b)
            dispersionEquation = (
                kr**2 + ks_wkb**2 + kp_wkb**2 - 1 / lambda_De_wkb**2 * (1 - ks_wkb * rho_TLe / (kp_wkb * np.sqrt(2)) * (A1_wkb * BesselProd0 * plasma_disp(z0e_wkb) + 
                    A2_wkb * (plasma_disp(z0e_wkb) * (1 + eval_b + z0e_wkb**2) * np.exp(-eval_b) + BesselProd1 * eval_b + z0e_wkb * BesselProd0))))

            rho_TLi = vTi_wkb / omci_wkb
            z0i_wkb = -(om_E_wkb - omega - 1j * nui_wkb) / (kp_wkb * np.sqrt(2) * vTi_wkb)
            eval_bi = kperp_wkb * rho_TLi**2
            if eval_bi.real>100:
                BesselProd0=1/(np.sqrt(2*np.pi*eval_bi))
                BesselProd1=1/(np.sqrt(2*np.pi*eval_bi)*(1+1/eval_bi**2)**(1/4))*np.exp(np.arcsinh(-1/eval_bi)+eval_bi*np.sqrt(1+1/eval_bi**2)-eval_bi)
            else:
                BesselProd0=Bessel(0,eval_bi)*np.exp(-eval_bi)
                BesselProd1=Bessel(-1,eval_bi)*np.exp(-eval_bi)
            dispersionEquation -= (
                1 / lambda_Di_wkb**2
                * (1 - ks_wkb * rho_TLi / (kp_wkb * np.sqrt(2))
                * (
                    A1i_wkb * BesselProd0 * plasma_disp(z0i_wkb) + 
                    A2i_wkb * (
                        plasma_disp(z0i_wkb) * (1 + eval_b + z0i_wkb**2) * np.exp(-eval_bi)
                        + BesselProd1 * eval_bi
                        + z0i_wkb * BesselProd0
                    )
                )
                )
            )
        elif mode == "horton":
            ky = ks_wkb
            om_prime = omega - om_E_wkb
            # Do WKB
            (
                B0_wkb,
                n_prof_wkb,
                ni_prof_wkb,
                Te_prof_wkb,
                Ti_prof_wkb,
                dndr_wkb,
                dnidr_wkb,
                dTedr_wkb,
                dTidr_wkb,
            ) = list(
                map(
                    lambda x: np.interp(r_prof[position],r_prof,x),
                    [B0, n_prof, ni_prof, Te_prof, Ti_prof, dndr, dnidr, dTedr, dTidr],
                )
            )
            om_ne = (
                -ky * sol / (e_charge * B0_wkb * n_prof_wkb) * Te_prof_wkb * ev * dndr_wkb
            )  #with -?
            om_ni = (
                ky
                * sol
                / (e_charge * Zi * B0_wkb * ni_prof_wkb)
                * Ti_prof_wkb
                * ev
                * dnidr_wkb
            )  
            om_Te = -ky * sol / (e_charge * B0_wkb) * dTedr_wkb * ev  #with -?
            om_Ti = ky * sol / (e_charge * Zi * B0_wkb) * dTidr_wkb * ev  
            ze = om_prime / (np.sqrt(2) * kp_wkb * vTe_wkb)
            zi = om_prime / (np.sqrt(2) * kp_wkb * vTi_wkb)
            W_e = (
                om_ne - om_prime + om_Te * (ze**2 + 1 / 2)
            ) * ze / om_prime * plasma_disp(ze) + om_Te * ze**2 / om_prime
            W_i = (
                om_ni - om_prime + om_Ti * (zi**2 + 1 / 2)
            ) * zi / om_prime * plasma_disp(zi) + om_Ti * zi**2 / om_prime
            nom_e = (
                1 / lambda_De_wkb**2 * (W_e - 1 - om_Te * ze / om_prime * plasma_disp(ze))
            )
            nom_i = (
                1 / lambda_Di_wkb**2 * (W_i - 1 - om_Ti * zi / om_prime * plasma_disp(zi))
            )
            denom_e = W_e * vTe_wkb**2 / (omce_wkb**2 * lambda_De_wkb**2)
            denom_i = W_i * vTi_wkb**2 / (omci_wkb**2 * lambda_Di_wkb**2)
            dispersionEquation = (
                (nom_e + nom_i - kp_wkb**2) / (1 + denom_e + denom_i) - ks_wkb**2 - kr**2 
            )
        return dispersionEquation

    
    def calc_dispersion_relation_k_of_r(self, mode=mode):

        self.load_all_profs()
        #self.calc_all_derivs()
        self.calc_parameters()
        self.find_res_surface()
        
        idx = int(self.general_dat['prof_length'] * self.options['r_per'])
        idx_range = np.linspace(self.options['r_range_start'], self.general_dat['prof_length'] - 1, self.options['n_points'])

        res = []
        r_used = []
        
        contour=cx.Rectangle([-self.contour_limit,self.contour_limit],[-self.contour_limit,self.contour_limit])

        for r in idx_range:
            #equation_k = lambda k: self.createDispersionEquation(k, self.options['omega'], int(r), mode=self.options['mode'])
            equation_k = lambda k: self.calc_dispersion_equation_for_kr_values_at_r(k, int(r), mode=mode)
            iterations=0
            roots_number=contour.count_roots(equation_k)
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
                    print("Abort. Too many iterations")
                    break
            if self.options['der']:
                roots = contour.roots(equation_k, df=lambda k: complex(grad(equation_k, holomorphic=True)(k)), guess_roots_symmetry=lambda z: [-z],verbose=True)
            else:
                roots = contour.roots(equation_k,guess_roots_symmetry=lambda z: [-z],int_method=self.options['int_method'],verbose=True,int_abs_tol=0.1,root_err_tol=1e-3)
            res.append(roots)
            r_used.append(self.r_prof[int(r)])
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

        self.r_found = r_found
        self.k_r1 = k_r1
        self.k_r2 = k_r2

        self.save_found_kr_of_r(mode)

    def save_found_kr_of_r(self, mode):
        file_name = mode + "_" + self.options['prof'] + "_" + "k_of_r"
        if self.der:
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


    def calculate_dispersion_relation_and_plot(self, mode=mode):

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
        if self.omegaOfk:
            for w in self.omega_range:
                equation_k = lambda k: self.createDispersionEquation(k, w, int(idx), mode=mode)
                contour = cx.Circle(0, self.options['radius'])
                if self.der:
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

            if self.save:
                if self.der:
                    self.save2txt(
                        omega_found,
                        k_w1,
                        k_w2,
                        self.mode + "_" + self.prof + "_" + "w(k)" + "_jax" + ".txt",
                    )
                else:
                    self.save2txt(
                        omega_found, k_w1, k_w2, self.mode + "_" + self.prof + "_" + "w(k)" + ".txt"
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

        if self.kOfr:
            for r in idx_range:
                equation_k = lambda k: self.createDispersionEquation(k, self.options['omega'], int(r), mode=self.options['mode'])
                #equation_k = lambda k: createDispersionEquationKIM(k, omega, r_prof[int(r)])
                iterations=0
                roots_number=contour.count_roots(equation_k)
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
                if self.der:
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
                r_used.append(self.r_prof[int(r)])
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

            if self.save:
                if self.der:
                    self.save2txt(
                        r_found,
                        k_r1,
                        k_r2,
                        self.options['mode'] + "_" + self.prof + "_" + "k(r)" + "_" + "_jax" + ".txt",
                    )
                else:
                    self.save2txt(r_found, k_r1, k_r2, self.options['mode'] + "_" + self.prof + "_" + "k(r)" + ".txt")

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

#######################################################################
    # dispose later:
    def calcParameters(self):
        ##Electrons
        Lee = (
            23.5
            - np.log(np.sqrt(self.n_prof) * self.Te_prof ** (-5 / 4))
            - np.sqrt(1e-5 + ((np.log(self.Te_prof) - 2) ** 2) / 16)
        )  # Coulomb logarithm
        vTe = np.sqrt(self.Te_prof * ev / e_mass)  # Thermal velocity
        omce = e_charge * self.B0 / (e_mass * sol)  # Cyclotron frequency btor or B0?
        nue = 5.8e-6 * self.n_prof * Lee * self.Te_prof ** (-3 / 2)  # Collision frequency
        lambda_De = np.sqrt(
            self.Te_prof * ev / (4 * np.pi * self.n_prof * e_charge**2)
        )  # Debye length
        A1 = (
            self.dndr / self.n_prof + e_charge / (self.Te_prof * ev) * self.E_prof - 3 / (2 * self.Te_prof) * self.dTedr
        )  # first thermodynamic force
        A2 = self.dTedr / self.Te_prof  # second thermodynamic force

        ks = (
            self.m_mode * self.hz - self.n_mode * self.hth / self.R0
        ) / self.r_prof  #'senkrecht' wavenumber ->look it up
        kp = self.m_mode / self.r_prof * self.hth + self.n_mode / self.R0 * self.hz  # parallel wavenumber ->look it up

        om_E = -sol * ks * self.E_prof / self.B0  # ExB rotation frequency
        ##Ions
        Lei = 24 - np.log(
            np.sqrt(self.n_prof) / self.Te_prof
        )  # Coulomb logarithm electrons-ions -> in KIM Ti_prof instead of Te_prof
        vTi = np.sqrt(ev * self.Ti_prof / (p_mass*self.Ai))  # Thermal velocity
        omci = e_charge * self.Zi * self.B0 / (p_mass * self.Ai * sol)  # Cyclotron frequency
        nue = nue + 7.7e-6 * self.ni_prof * Lei * self.Zi**2 * self.Te_prof ** (
            -3 / 2
        )  # Add ions to electrons collision frequency
        nui = (
            1.8e-7 * self.Ai ** (-1 / 2) * self.Ti_prof ** (-3 / 2) * self.n_prof * self.Zi**2 * Lei
        )  # Collision frequency ions with electrons
        Lii = 23 - np.log(
            self.Zi**2 * self.Ai * 2 / ((self.Ti_prof * self.Ai) * 2) * np.sqrt((self.ni_prof * self.Zi**2 / (self.Ti_prof)) * 2)
        )  # Coulomb logarithm ions-ions
        nui = nui + 1.8e-7 * self.ni_prof * self.Zi**4 * Lii * self.Ai ** (-1 / 2) * self.Ti_prof ** (
            -3 / 2
        )  # Add ion ion collisions
        lambda_Di = np.sqrt(
            self.Ti_prof * ev / (4 * np.pi * self.ni_prof * (e_charge * self.Zi) ** 2)
        )  # Debye length
        A1i = (
            self.dnidr / self.ni_prof
            - (e_charge * self.Zi) / (self.Ti_prof * ev) * self.E_prof
            - 3 / (2 * self.Ti_prof) * self.dTidr
        )  # first thermodynamic force
        A2i = self.dTidr / self.Ti_prof  # second thermodynamic force
        if self.noCollisions:
            nue=np.zeros(self.prof_length)
            nui=np.zeros(self.prof_length)
        return (
            vTe,
            omce,
            nue,
            lambda_De,
            A1,
            A2,
            ks,
            kp,
            om_E,
            vTi,
            omci,
            nui,
            lambda_Di,
            A1i,
            A2i,
        )

    def print_species_data(self, spec):
        print("")
        print(f"Species: {spec}")
        [print(f"{key}: {value}") for key, value in self.spec_dat[spec].items()]


if __name__ == "__main__":

    if len(sys.argv) > 1:
        assert sys.argv[1] in ['horton', 'KIM2'], "Mode not supported"
        mode = sys.argv[1]
    else:
        mode = 'horton'
    
    print('Mode is ', mode)

    kwkb = KIM_WKB()
    kwkb.prof_path = '../../../kim-wkb/profiles_parab/'

    #kwkb.load_all_profs()
    #kwkb.calc_all_derivs()
    #kwkb.calc_equilibrium_fields()
    #kwkb.calc_all_parameters()
    #kwkb.calc_parameters()
    #nue, nui = kwkb.calc_collision_frequency()

    kwkb.calc_dispersion_relation_k_of_r(mode=mode)
    kwkb.plot_kr_of_r()
    #kwkb.plot_from_file('horton_parab_k_of_r.txt')

    #kwkb.load_all_profs()
    #kwkb.calc_parameters()
    #kwkb.find_res_surface()
    #for spec in kwkb.species:
        #kwkb.print_species_data(spec)



    #kwkb.calculate_dispersion_relation_and_plot()
