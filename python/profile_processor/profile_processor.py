import numpy as np
import matplotlib.pyplot as plt
import re
import os
import h5py
from scipy.interpolate import interp1d
import sys
sys.path.append(os.path.dirname(__file__) + '/../fieldpy/')
from fieldpy import fieldpy
sys.path.append(os.path.dirname(__file__) + '/../neo2_for_Er/')
from neo2_for_Er import neo2_for_Er

from scipy.integrate import solve_ivp


class profile_processor:

    flux_data = ''

    kB = 1.3807e-16
    eVK = 1.1604e4
    c = 29979245800.0
    echarge = 4.8031e-10
    eV_to_erg = 1.6022e-12

    m_mode = [11,12,13,14,15,16]
    n_mode = [2]

    def __init__(self, runpath):
        self.runpath = runpath
        self.plot_dir = self.runpath + '/plots/'

        os.makedirs(self.runpath + '/profiles/', exist_ok=True)

    def run_fieldpy(self, gfile, convex_wall, flux_data, pfile='', skip=False):
        '''Run field_divB0 to get equilibrium data containing r_eff, q, psi (pol. flux), phi (tor. flux) and more.
            input:
                gfile ... EQDSK-standard equilibrium file
                pfile ... field.dat, only necessary if vacuum field perturbation is required
                convex_wall ... convex wall for field calc. should be outside LCFS and inside RMP coils 
                flux_data ... output directory for the flux data
        '''

        self.gfile = gfile
        self.pfile = pfile
        self.convex_wall = convex_wall
        self.flux_data = flux_data

        if not skip:
            self.fp = fieldpy(self.gfile, self.pfile, self.convex_wall, self.flux_data)
            self.fp.write_field_divB0_inp(self.fp.path_to_fourier_modes_exe + 'template_field_divB0.inp', self.fp.path_to_fourier_modes_exe + 'field_divB0.inp')
            self.fp.run_fourier_modes()


    def map_profs_to_reff(self, prof_path, save_path='', flux_data='', plot=True):
        '''Read in the equil data with r_eff, q and psi (poloidal flux)
            input:
                prof_path ... path of the kinetic profiles
                flux_data ... path to the equil_r_q_psi.dat file 
                save_path ... path where the profiles should be saved to. If empty, save to prof_path
            '''

        def remove_non_monotonous_tail(data):
            non_mon_index = next((i for i, (x,y,) in enumerate(zip(data[:-1], data[1:])) if x > y), None)
            if non_mon_index is not None:
                return [data[:non_mon_index+1], non_mon_index]
                print('Motonous tail found and removed')
            else:
                return [data, len(data) - 1]

        if not flux_data == '' and self.flux_data == '':
            self.flux_data = flux_data
        
        if save_path == '':
            save_path = prof_path
        
        self.prof_path = prof_path
        self.save_path = save_path

        self.equil_data = np.loadtxt(flux_data + 'equil_r_q_psi.dat')

        self.r_eff = self.equil_data[:,0]
        [self.r_eff, ind] = remove_non_monotonous_tail(self.r_eff)

        #if ind == len(self.equil_data[:,0])-1:
        #    ind = ind-1

        self.q = self.equil_data[:ind+1,1]
        self.psi = self.equil_data[:ind+1,2]

        self.psi_n = self.psi / self.psi[-1]
        self.s = np.sqrt(self.psi_n)

        marsf_pattern = re.compile(r'^PROF*')
        kilca_pattern = r'\.dat$'
        

        for filename in os.listdir(prof_path):
            if filename == 'kprof':
                continue
            if marsf_pattern.match(filename):
                self.ne_orig = np.loadtxt(prof_path + 'PROFDEN.IN',skiprows=1)
                self.Te_orig = np.loadtxt(prof_path + 'PROFTE.IN', skiprows=1)
                self.Ti_orig = np.loadtxt(prof_path + 'PROFTI.IN', skiprows=1)
                self.Vz_orig = np.loadtxt(prof_path + 'PROFROT.IN', skiprows=1)

                #rescale Vz since it is given in m/s instead of cm/s
                self.Vz_orig[:,1] = self.Vz_orig[:,1] * 1e2

                break
            elif re.search(kilca_pattern, filename):
                self.ne_orig = np.loadtxt(prof_path + 'n.dat')
                self.Te_orig = np.loadtxt(prof_path + 'Te.dat')
                self.Ti_orig = np.loadtxt(prof_path + 'Ti.dat')
                self.Vz_orig = np.loadtxt(prof_path + 'Vz.dat')
                break
            else:
                raise ValueError('Other input profiles not yet implemented')
        if int(np.floor(np.log10(np.abs(self.ne_orig[0,1])))) > 16:
            # rescale since density is given in cubic meter instead of centimeter
            self.ne_orig[:,1] = self.ne_orig[:,1] * 1e-6


        #self.ne = np.interp(self.s, self.ne_orig[:,0], self.ne_orig[:,1])
        spline_ne = interp1d(self.ne_orig[:,0], self.ne_orig[:,1], kind='cubic')
        self.ne = spline_ne(self.s)

        #self.Te = np.interp(self.s, self.Te_orig[:,0], self.Te_orig[:,1])
        spline_Te = interp1d(self.Te_orig[:,0], self.Te_orig[:,1], kind='cubic')
        self.Te = spline_Te(self.s)

        spline_Ti = interp1d(self.Ti_orig[:,0], self.Ti_orig[:,1], kind='cubic')
        #self.Ti = np.interp(self.s, self.Ti_orig[:,0], self.Ti_orig[:,1])
        self.Ti = spline_Ti(self.s)

        #self.Vz = np.interp(self.s, self.Vz_orig[:,0], self.Vz_orig[:,1])
        spline_Vz = interp1d(self.Vz_orig[:,0], self.Vz_orig[:,1], kind='cubic')
        self.Vz = spline_Vz(self.s)

        np.savetxt(save_path + 'n.dat', np.column_stack((self.r_eff, self.ne)))
        np.savetxt(save_path + 'Te.dat', np.column_stack((self.r_eff, self.Te)))
        np.savetxt(save_path + 'Ti.dat', np.column_stack((self.r_eff, self.Ti)))
        np.savetxt(save_path + 'Vz.dat', np.column_stack((self.r_eff, self.Vz)))
        np.savetxt(save_path + 'q.dat', np.column_stack((self.r_eff, self.q)))
        print('Profiles written to ' + save_path)

        if plot==True:
            plt.figure()
            plt.plot(self.r_eff, self.ne/np.max(self.ne), label=r'n$_e$')
            plt.plot(self.r_eff, self.Te/np.max(self.Te), label=r'T$_e$')
            plt.plot(self.r_eff, self.Ti/np.max(self.Ti), label=r'T$_i$')
            plt.plot(self.r_eff, self.Vz/np.max(self.Vz), label=r'V$_z$')
            plt.plot(self.r_eff, np.abs(self.q)/np.max(np.abs(self.q)), label=r'|q|')
            plt.legend(bbox_to_anchor=(1.0,1.0))
            plt.xlabel('r [cm]')
            plt.ylabel('normalized profiles')
            plt.minorticks_on()
            plt.grid(which='major')
            plt.tight_layout()
            plt.show()

    def determine_anomalous_diff_coeff(self):
        
        self.Da = np.ones_like(self.r_eff) * 1e4

        np.savetxt(self.save_path + 'Da.dat', np.column_stack((self.r_eff, self.Da)))

    def calc_Er_prof(self, recalc=False):
        '''Caluclate Er profile with NEO-2.'''

        dat = np.loadtxt(self.flux_data + 'btor_rbig.dat')
        self.Btor = dat[0]
        self.R0 = dat[1]

        self.solve_cyl_equilibrium()

        if not os.path.exists(self.save_path + 'kprof/k.dat') or recalc:
            # don't calculate k if it exists and if the recalculation is not needed
            self.neo2 = neo2_for_Er(self.save_path, self.flux_data + 'equil_r_q_psi.dat')
            self.neo2.run_neo2(self.gfile, self.convex_wall, self.flux_data)
            self.collect_k_profile(self.save_path + 'kprof/', self.save_path + 'kprof/')
        else:
            # read k profile
            k_dat = np.loadtxt(self.save_path + '/kprof/k.dat')
            self.k_prof = k_dat[:,1]
            self.k_prof_r = k_dat[:,0]

        self.press_ion = self.ne * self.Ti * self.kB * self.eVK
        self.dpress_ion = np.gradient(self.press_ion, self.r_eff)

        self.dTi = np.gradient(self.Ti, self.r_eff) * self.eV_to_erg
        self.dne = np.gradient(self.ne, self.r_eff)

        #self.v_hat = self.c * self.Bz * self.dTi / (self.echarge * self.B**2)

        self.k_prof_r = np.append(self.k_prof_r, self.r_eff[-1])
        self.k_prof = np.append(self.k_prof, 0.5 * self.k_prof[-1])
        #self.k = np.interp(self.r_eff, self.k_prof_r, self.k_prof)
        spline_k = interp1d(self.k_prof_r, self.k_prof, kind='cubic', fill_value='extrapolate')
        self.k = spline_k(self.r_eff)

        #self.vth = self.k * self.v_hat

        #self.Vpol = self.k * self.c * self.dTi / self.echarge/ self.ne

        self.Er = self.Ti * self.eV_to_erg * self.dne / (self.echarge * self.ne) + (1.0 - self.k) * self.dTi / self.echarge + self.r_eff * self.B * self.Vz / (self.c * self.q * self.R0)

        #self.Er = 1.0 / (self.echarge * self.ne) * self.dpress_ion - self.k * self.dTi / self.echarge
        #self.Er = self.Er + self.r_eff * self.Vz * self.B / (self.R0 * self.c * self.q)
        np.savetxt(self.save_path + 'Er.dat', np.column_stack((self.r_eff, self.Er)))



    def collect_k_profile(self, kpath, fname):
        """Collect the k profile from the NEO-2 output."""
        wd = os.getcwd()
        os.chdir(kpath)

        content = os.listdir(kpath)
        content = content[2:-1]

        k = np.array([])
        s = np.array([])

        abort = 0

        for j, con in enumerate(content):
            #print(con)
            if os.path.exists(con) and os.path.isdir(con) and not os.path.basename(con) == 'TEMPLATE_DIR':
                try:
                    fulltransp = h5py.File(con + '/fulltransp.h5')
                    neo2config = h5py.File(con + '/neo2_config.h5')

                    try:
                        if np.isnan(np.array(fulltransp['k_cof'])):
                            continue
                        k = np.append(k, np.array(fulltransp['k_cof']))
                    except:
                        ntvout = np.loadtxt(con + '/ntv_out.dat')
                        k = ntvout[6]

                    s = np.append(s, np.array(neo2config['settings']['boozer_s']))
                except:
                    abort = abort + 1
        sort_ind = np.argsort(s)
        s = np.sort(s)
        k = k[sort_ind]

        M = np.loadtxt('surfaces.dat')
        r = np.interp(s, M[:,0], M[:,1])

        np.savetxt(fname + 'k.dat', np.column_stack((r,k)))
        os.chdir(wd)

        self.k_prof = k
        self.k_prof_r = r


    def solve_cyl_equilibrium(self):

        r_ode = self.r_eff
        q_ode = - self.q
        p_tot = self.ne * (self.Te + self.Ti) * self.kB * self.eVK

        g_ode = 1.0 + (r_ode / -q_ode / self.R0)**2

        dp_tot = np.gradient(p_tot, r_ode)

        fq = lambda x: np.interp(x, r_ode, q_ode)
        fg = lambda x: np.interp(x, r_ode, g_ode)
        fdp = lambda x: np.interp(x, r_ode, dp_tot)

        odefun = lambda x,y: -2.0 * x * y / fq(x)**2 / fg(x) / self.R0**2 - 8.0 * np.pi * fdp(x)

        u0 = self.Btor**2 * g_ode[0]

        u_ode = solve_ivp(odefun, (r_ode[0], r_ode[-1]), [u0], method='RK45', t_eval=r_ode)

        Bz_ode = np.sign(self.Btor) * np.sqrt(u_ode.y[0] / g_ode)
        Bth_ode = r_ode * Bz_ode / q_ode / self.R0

        self.Bz = Bz_ode
        self.Bth = Bth_ode
        self.B = np.sqrt(Bz_ode**2 + Bth_ode **2)
        print(f'Cylinder B = {self.B[0]}')


    def get_resonant_radii(self, m_mode, n_mode):
        """Get the effective radius of the resonant surfaces for which the mode numbers are given.
            input:
                - m_mode ... numpy array with poloidal mode numbers
                - n_mode ... single number of toroidal mode number 
            """
        r_res = []
        for i, m in enumerate(m_mode):
            r_res.append(np.interp(m/n_mode, self.q, self.r_eff))
        return r_res


    def plot_Er_profile(self, save=False):
        plt.figure()
        plt.plot(self.r_eff, self.Er, lw=2)

        plt.xlabel(r'r$_{\mathrm{eff}}$ [cm]')
        plt.ylabel(r'E$_r$ [statV cm$^{-1}$]')

        res = []
        for i, m in enumerate(self.m_mode):
            res.append(np.interp(m/self.n_mode[0], np.abs(self.q), self.r_eff))
            plt.axvline(res[i], ls=':', c='dimgrey')
            plt.text(res[i], plt.ylim()[1]- 0.03 - np.mod(i,3) * 0.025, f'{self.m_mode[i]}', horizontalalignment='center')
        ax = plt.gca()
        ax.tick_params(axis='both', which='both', direction='out', top=True, right=True)
        plt.grid()
        plt.minorticks_on()
        plt.axhline(0.0, ls='--', c='k', lw=1.5)
        plt.tight_layout()
        if save:
            plt.savefig(self.plot_dir + 'Er_of_reff.pdf', transparent=False, bbox_inches='tight')
        plt.show()

    def plot_profiles(self, save=True):
        
        os.makedirs(self.plot_dir, exist_ok=True)

        plt.rc('font', size=12)
        self.plot_Er_profile(True)
        self.plot_kinetic_profs(True)


    def plot_kinetic_profs(self, save=False):
        fig,ax = plt.subplots(2,2, figsize=(8,6))

        ax[0,0].plot(self.r_eff, self.ne, lw=2)
        ax[0,0].set_ylabel(r'n$_e$ [cm$^{-3}$]')
        ax[0,0].set_xlabel(r'r$_{\mathrm{eff}}$ [cm]')

        ax[0,1].plot(self.r_eff, self.Te, lw=2)
        ax[0,1].set_ylabel(r'T$_e$ [eV]')
        ax[0,1].set_xlabel(r'r$_{\mathrm{eff}}$ [cm]')

        ax[1,0].plot(self.r_eff, self.Ti, lw=2)
        ax[1,0].set_ylabel(r'T$_i$ [eV]')
        ax[1,0].set_xlabel(r'r$_{\mathrm{eff}}$ [cm]')

        ax[1,1].plot(self.r_eff, self.Vz, lw=2)
        ax[1,1].set_ylabel(r'V$_{\mathrm{tor}}$ [cm s$^{-1}$]')
        ax[1,1].set_xlabel(r'r$_{\mathrm{eff}}$ [cm]')

        list(map(lambda x: x.tick_params(axis='both', which='both', direction='out', top=True, right=True), ax[:,0]))
        list(map(lambda x: x.tick_params(axis='both', which='both', direction='out', top=True, right=True), ax[:,1]))
        
        list(map(lambda x: x.grid(), ax[:,0]))
        list(map(lambda x: x.grid(), ax[:,1]))
        list(map(lambda x: x.minorticks_on(), ax[:,0]))
        list(map(lambda x: x.minorticks_on(), ax[:,1]))


        res = []
        for i, m in enumerate(self.m_mode):
            res.append(np.interp(m/self.n_mode[0], np.abs(self.q), self.r_eff))
            list(map(lambda x: x.axvline(res[i], ls=':', c='dimgrey'), ax[:,0]))
            list(map(lambda x: x.axvline(res[i], ls=':', c='dimgrey'), ax[:,1]))
            list(map(lambda x: x.text(res[i], x.get_ylim()[1]*(1.0 - np.mod(i,3) * 0.25), f'{self.m_mode[i]}', horizontalalignment='center'), ax[:,0]))
            list(map(lambda x: x.text(res[i], x.get_ylim()[1]*(1.0 - np.mod(i,3) * 0.25), f'{self.m_mode[i]}', horizontalalignment='center'), ax[:,1]))

        plt.tight_layout()
        if save:
            plt.savefig(self.plot_dir + 'kinetic_of_reff.pdf', transparent=False, bbox_inches='tight')
        plt.show()

