# %%
import cxroots as cx
import numpy as np
import matplotlib.pyplot as plt

from plasmapy.dispersion import plasma_dispersion_func as plasma_disp

import os
from scipy.integrate import solve_ivp
from scipy.special import iv as Bessel
from jax import grad
import jax
import jax.numpy as jnp


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
    btor = -17977.413  # toroidal magnetic field in Gauss
    R0 = 165.0  # major radius in cm
    m_mode = 6  # poloidal mode number
    n_mode = 2  # toroidal mode number
    Zi = 1  # Ion charge number
    Ai = 2  # Ion mass number
    r_plas = 64.7  # plasma radius

    ##Options
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
    # Choose model: KIM or horton
    mode = "KIM2"
    # Choose method for integration: quad or romb
    int_method = "quad"
    # Should the derivative be approximated with jax.grad
    der = False
    # Should the results be saved
    save = True
    noCollisions=True

    species = ['e', 'i'] # could be extended to multiple species, e.g. ['e', 'D', 'He'], ...

    general_dat = {}
    general_dat_keys = ['Er', 'Vz']

    spec_keys = ['n', 'T', 'vT', 'nu']
    spec_dat_units = ['cm', 'erg', 'cm/s', '1/s']
    spec_dat = {}

    equil_keys = ['B0z', 'B0th', 'B0', 'hz', 'hth']
    equil_dat = {}

    def __init__(self):
        for spec in self.species:
            self.spec_dat[spec] = {}

    # Load profiles
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


    # Calc equilibrium field
    def calcEquilibrium(self):
        # Pressure due to electrons
        dpress_prof = ev * (self.dndr * self.Te_prof + self.n_prof * self.dTedr)
        # Pressure due to ions
        press_prof = ev * self.n_prof * (self.Te_prof + self.Ti_prof)
        dpress_prof += ev * (self.dnidr * self.Ti_prof + self.ni_prof * self.dTidr)
        # solve ivp
        #f = lambda r, u: dudr(r, u, r_prof, q_prof, dpress_prof)
        f = lambda r,u: -2*r*u/(self.R0**2 * self.q_prof[self.val2ind(r, self.r_prof)]**2 + r**2)-8*pi*dpress_prof[self.val2ind(r, self.r_prof)]
        u0 = self.btor**2 * (1 + self.r_prof[0] ** 2 / (self.R0**2 * self.q_prof[0] ** 2))
        u = solve_ivp(f, t_span=[self.r_prof[0], self.r_prof[-1]], y0=[u0], t_eval=self.r_prof)
        u = u.y[0]
        B0z = np.sign(self.btor) * np.sqrt(u / (1 + self.r_prof**2 / (self.R0**2 * self.q_prof**2)))
        B0th = B0z * self.r_prof / (self.q_prof * self.R0)
        B0 = np.sqrt(B0th**2 + B0z**2)
        hz = B0z / B0
        hth = B0th / B0
        return u, B0z, B0th, B0, hz, hth


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


    # Calc parameters and derivatives
    def calcDerivatives(self):
        dndr = np.zeros(prof_length)
        dTedr = np.zeros(prof_length)
        dqdr = np.zeros(prof_length)
        dnidr = np.zeros(prof_length)
        dTidr = np.zeros(prof_length)
        for k in range(prof_length - 1):
            dndr[k] = (n_prof[k + 1] - n_prof[k]) / (r_prof[k + 1] - r_prof[k])
            dTedr[k] = (Te_prof[k + 1] - Te_prof[k]) / (r_prof[k + 1] - r_prof[k])
            dqdr[k] = (q_prof[k + 1] - q_prof[k]) / (r_prof[k + 1] - r_prof[k])
            dnidr[k] = (ni_prof[k + 1] - ni_prof[k]) / (r_prof[k + 1] - r_prof[k])
            dTidr[k] = (Ti_prof[k + 1] - Ti_prof[k]) / (r_prof[k + 1] - r_prof[k])
        dndr[prof_length - 1] = dndr[prof_length - 2]
        dTedr[prof_length - 1] = dTedr[prof_length - 2]
        dqdr[prof_length - 1] = dqdr[prof_length - 2]
        dnidr[prof_length - 1] = dnidr[prof_length - 2]
        dTidr[prof_length - 1] = dTidr[prof_length - 2]
        return dndr, dTedr, dqdr, dnidr, dTidr


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
            self.Zi**2 * self.Ai * 2 / ((self.Ti_prof * self.Ai) * 2) * (self.ni_prof * self.Zi**2 / (self.Ti_prof)) * 2
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
        # Load profiles
        self.r_prof = self.load_prof("r", self.prof)
        self.E_prof = self.load_prof("Er.dat", self.prof)
        self.n_prof = self.load_prof("n.dat", self.prof)
        self.ni_prof = self.n_prof  # Quasineutrality for single ion species
        self.q_prof = self.load_prof("q.dat", self.prof)
        self.Te_prof = self.load_prof("Te.dat", self.prof)
        self.Ti_prof = self.load_prof("Ti.dat", self.prof)
        self.Vth_prof = self.load_prof("Vth.dat", self.prof)
        # Vz=load_prof('Vz.dat')
        self.prof_length = len(self.r_prof)
    
    def calc_all_derivs(self):
        # Calculate derivatives
        #dndr, dTedr, dqdr, dnidr, dTidr = calcDerivatives()
        self.dndr, self.dTedr, self.dqdr, self.dnidr, self.dTidr = list(map(lambda x: np.gradient(x,self.r_prof),[self.n_prof,self.Te_prof,self.q_prof,self.ni_prof,self.Ti_prof]))

    def calc_equilibrium_fields(self):
        # Calculate equilibrium fields
        self.u, self.B0z, self.B0th, self.B0, self.hz, self.hth = self.calcEquilibrium()
        # Calculate parameters
        self.vTe, self.omce, self.nue, self.lambda_De, self.A1, self.A2, \
            self.ks, self.kp, self.om_E, self.vTi, self.omci, self.nui, \
            self.lambda_Di, self.A1i, self.A2i = (
            self.calcParameters()
        )

    def find_res_surface(self):
        # Find resonant surface
        self.q = self.r_prof * self.B0z / (self.R0 * self.B0th)
        self.idx_res = self.val2ind(-self.m_mode / self.n_mode, self.q_prof)


    # Create dispersion equation
    def createDispersionEquation(self, kr, omega, position, mode="KIM"):
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
        m_phi = 0
        kperp_wkb = ks_wkb**2 + kr**2
        # Use either KIM or horton
        if mode == "KIM":
            z0e_wkb = -(om_E_wkb - omega - 1j * nue_wkb) / (kp_wkb * np.sqrt(2) * vTe_wkb)
            # In WKB bp and bt are the same
            eval_b = kperp_wkb**2 * vTe_wkb**2 / omce_wkb**2
            if eval_b.real>100:
                BesselProd0=1/(np.sqrt(2*np.pi*eval_b))
                BesselProd1=1/(np.sqrt(2*np.pi*eval_b)*(1+1/eval_b**2)**(1/4))*np.exp(np.arcsinh(-1/eval_b)+eval_b*np.sqrt(1+1/eval_b**2))*np.exp(-eval_b)
            else:
                BesselProd0=1/(np.sqrt(2*np.pi*eval_b))
                BesselProd1=1/(np.sqrt(2*np.pi*eval_b)*(1+1/eval_b**2)**(1/4))*np.exp(np.arcsinh(-1/eval_b)+eval_b*np.sqrt(1+1/eval_b**2)-eval_b)
            a0 = BesselProd0 * (
                -m_phi
                - om_E_wkb / omce_wkb
                + ks_wkb
                * vTe_wkb**2
                / omce_wkb**2
                * (A1_wkb + (1 + eval_b + m_phi) * A2_wkb)
            ) + ks_wkb * vTe_wkb**2 / omce_wkb**2 * A2_wkb * eval_b * BesselProd1
            a1 = -kp_wkb / omce_wkb *BesselProd0
            a2 = (
                ks_wkb
                / (2 * omce_wkb**2)
                * A2_wkb
                * BesselProd0
            )
            dispersionEquation = (
                kr**2
                + ks_wkb**2
                + kp_wkb**2
                - 1
                / lambda_De_wkb**2
                * (
                    BesselProd0
                    + omce_wkb
                    / (np.sqrt(2) * kp_wkb * vTe_wkb)
                    * (
                        plasma_disp(z0e_wkb)
                        * (
                            a0
                            + np.sqrt(2) * vTe_wkb * z0e_wkb * a1
                            + vTe_wkb**2 * a2 * 2 * z0e_wkb**2
                        )
                        + vTe_wkb * np.sqrt(2) * (a1 + vTe_wkb * np.sqrt(2) * z0e_wkb * a2)
                    )
                )
            )
            ##Calculate dispersion equation for ions
            z0i_wkb = -(om_E_wkb - omega - 1j * nui_wkb) / (kp_wkb * np.sqrt(2) * vTi_wkb)
            # In WKB bp and bt are the same
            eval_bi = kperp_wkb * vTi_wkb**2 / omci_wkb**2
            if eval_bi.real>100:
                BesselProd0=1/(np.sqrt(2*np.pi*eval_bi))
                BesselProd1=1/(np.sqrt(2*np.pi*eval_bi)*(1+1/eval_bi**2)**(1/4))*np.exp(np.arcsinh(-1/eval_bi)+eval_bi*np.sqrt(1+1/eval_bi**2)-eval_bi)
            else:
                BesselProd0=Bessel(0,eval_bi)*np.exp(-eval_bi)
                BesselProd1=Bessel(-1,eval_bi)*np.exp(-eval_bi)
            a0 = BesselProd0 * (
                -m_phi
                - om_E_wkb / omci_wkb
                + ks_wkb
                * vTi_wkb**2
                / omci_wkb**2
                * (A1i_wkb + (1 + eval_bi + m_phi) * A2i_wkb)
            ) + ks_wkb * vTi_wkb**2 / omci_wkb**2 * A2i_wkb * BesselProd1
            a1 = -kp_wkb / omci_wkb * BesselProd0
            a2 = (
                ks_wkb
                / (2 * omci_wkb**2)
                * A2i_wkb
                * BesselProd0
            )
            # Add ion part to dispersion equation
            dispersionEquation -= (
                1
                / lambda_Di_wkb**2
                * (
                    BesselProd0
                    + omci_wkb
                    / (np.sqrt(2) * kp_wkb * vTi_wkb)
                    * (
                        plasma_disp(z0i_wkb)
                        * (
                            a0
                            + np.sqrt(2) * vTi_wkb * z0i_wkb * a1
                            + vTi_wkb**2 * a2 * 2 * z0i_wkb**2
                        )
                        + vTi_wkb * np.sqrt(2) * (a1 + vTi_wkb * np.sqrt(2) * z0i_wkb * a2)
                    )
                )
            )
        elif mode == 'KIM2':
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
                kr**2
                + ks_wkb**2
                + kp_wkb**2
                - 1 / lambda_De_wkb**2
                * (1 - ks_wkb * rho_TLe / (kp_wkb * np.sqrt(2))
                * (
                    A1_wkb * BesselProd0 * plasma_disp(z0e_wkb) + 
                    A2_wkb * (
                        plasma_disp(z0e_wkb) * (1 + eval_b + z0e_wkb**2) * np.exp(-eval_b)
                        + BesselProd1 * eval_b
                        + z0e_wkb * BesselProd0
                    )
                )
                )
            )

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


    def createDispersionEquationKIM(self, kr, omega, r_pos):
        # Use rewritten kernel for dispersion equation derived by plugging in 
        # a_i coefficients.
        # Looks a bit simpler and does not use a_i abbreviations
        (
            om_E_wkb,
            vTe_wkb,
            vTi_wkb,
            lambda_De_wkb,
            lambda_Di_wkb,
            kp_wkb,
            ks_wkb,
            A1_wkb,
            A2_wkb,
            A1i_wkb,
            A2i_wkb,
        ) = list(
            map(
                lambda x: np.interp(r_pos, r_prof, x),
                [om_E, vTe, vTi, lambda_De, lambda_Di, kp, ks, A1, A2, A1i, A2i],
            )
        )

        omce = e_charge * B0 / (e_mass * sol)
        omci = e_charge * B0 / (p_mass * sol)
        omce_wkb = np.interp(r_pos, r_prof, omce)
        omci_wkb = np.interp(r_pos, r_prof, omci)
    
        m_phi = 0
        # Do WKB
        nue_wkb = 0#np.interp(r_pos, r_prof, nue)
        nui_wkb = 0#np.interp(r_pos, r_prof, nui)
        z0e_wkb = -(om_E_wkb - omega - 1j * nue_wkb) / (kp_wkb * np.sqrt(2) * vTe_wkb)
        rho_TLe = vTe_wkb / omce_wkb
        kperp_wkb = np.sqrt(ks_wkb**2 + kr**2)
    
        # In WKB bp and bt are the same
        eval_b = kperp_wkb**2 * rho_TLe**2
    
        dispersionEquation = (
            kr**2
            + ks_wkb**2
            + kp_wkb**2
            + 1 / lambda_De_wkb**2
            * (1 - np.exp(-eval_b)
            * ks_wkb * rho_TLe / (kp_wkb * np.sqrt(2))
            * (
                A1_wkb * Bessel(0, eval_b) * plasma_disp(z0e_wkb) + 
                A2_wkb * (
                    plasma_disp(z0e_wkb) * (1 + eval_b + z0e_wkb**2)
                    + Bessel(-1, eval_b) * eval_b
                    + z0e_wkb * Bessel(0, eval_b)
                )
            )
            )
        )

        rho_TLi = vTi_wkb / omci_wkb
        z0i_wkb = -(om_E_wkb - omega - 1j * nui_wkb) / (kp_wkb * np.sqrt(2) * vTi_wkb)
        eval_b = kperp_wkb**2 * rho_TLi**2

        dispersionEquation += (
            1 / lambda_De_wkb**2
            * (1 - np.exp(-eval_b)
            * ks_wkb * rho_TLi / (kp_wkb * np.sqrt(2))
            * (
                A1i_wkb * Bessel(0, eval_b) * plasma_disp(z0i_wkb) + 
                A2i_wkb * (
                    plasma_disp(z0i_wkb) * (1 + eval_b + z0i_wkb**2)
                    + Bessel(-1, eval_b) * eval_b
                    + z0i_wkb * Bessel(0, eval_b)
                )
            )
            )
        )
        return dispersionEquation
 
    

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
                contour = cx.Circle(0, radius)
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
                        int_method=self.int_method,
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
                equation_k = lambda k: self.createDispersionEquation(k, self.omega, int(r), mode=self.mode)
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
                        int_method=self.int_method,
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
                        self.mode + "_" + self.prof + "_" + "k(r)" + "_" + "_jax" + ".txt",
                    )
                else:
                    self.save2txt(r_found, k_r1, k_r2, self.mode + "_" + self.prof + "_" + "k(r)" + ".txt")

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


if __name__ == "__main__":
    kwkb = KIM_WKB()
    kwkb.calculate_dispersion_relation_and_plot()
