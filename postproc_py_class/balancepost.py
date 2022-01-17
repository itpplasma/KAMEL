##########################################################
# balancepost.py
##########################################################
# Description:
# Class for the post processing of the Balance circle.
# --------------------------------------------------------
# Methods:
# def loadfile(self):
# def shotinfo(self):
# def plotprofiles(self, export = 0):
# def showcontent(self):
# def plotdata(self, key, dataname, export = 0):
# def plotbifurcation(self, export = 0):
# --------------------------------------------------------
# To do:
# - Add second x axis to plt_shielding_fac_over_ant_fac that shows time
# - Adjust x label in plt_shielding_fac_wo_diag_all_modes_one_fig
# - Adjust y label in plt_antenna_ramp_up
##########################################################
# Author: Markus Markl
# Created: 31.03.2021
##########################################################

import os
from h5py._hl import group
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import h5py
from scipy.interpolate import CubicSpline
import re

class postproc:

    prof_t_plot_r_offset = 0.5  # this is in centimeters
    prof_t_plot_y_offset_lower = 1.0 # this is in percent of the profile
    prof_t_plot_y_offset_upper = 1.0 # this is in percent of the profile

    path2inp = ''
    path2out = ''

    profile_types = {'n': r'/$cm^{-3}$', 'Te': r'/$eV$', 'Ti': r'/$eV$', 'Vz': r'/$cm/s$', 'Er': r'/$statV cm^{-1}$'}

    def __init__(self, path2inp, path2out):
        """ Constructor of class.
        input: path2inp ... path to the input hdf5 file containing also the file name
               path2out ... path to the output hdf5 file containing also the file name"""

        self.path2inp = path2inp
        self.path2out = path2out
        self.loadfile()
        self.loadquantities()


    def loadfile(self):
        """Open the hdf5 input and output files to read."""
        self.h5inp = h5py.File(self.path2inp, 'r')
        self.h5out = h5py.File(self.path2out, 'r')


    def loadquantities(self):
        """ Load specific initial quantities, e.g. shot number, time, ..."""
        self.fac_n = np.array(self.h5out['factors/fac_n'])[0]
        self.fac_Ti = np.array(self.h5out['factors/fac_Ti'])[0]
        self.fac_Te = np.array(self.h5out['factors/fac_Te'])[0]
        self.fac_vz = np.array(self.h5out['factors/fac_vz'])[0]
        #self.balance = self.h5file['balance']
        self.shot = self.h5out['shot'][0][0]
        self.time = self.h5out['time'][0][0]
        self.m = self.h5out['m'][:][:,0]
        self.n = self.h5out['n'][:][:,0]
        self.r_res = self.h5inp['output/r_res'][:]
        self.da_res = self.h5inp['output/Da_res'][:]
        self.scalefactors_sq = self.h5inp['output/scalefactors_sq'][:]


    def get_scan_names(self):
        """ Get the names of the parameter scan groups"""
        self.scans_list = []
        for key in self.h5out.keys():
            if re.search('\An', key):
                self.scans_list.append(key)
        return self.scans_list

    def get_scan_names_from_factors(self):
        """ Use the factor values to generate the group names for the parameter scan."""
        #fac_Te = np.array(Tnscan.h5out['factors/fac_Te'])[0]
        #fac_Ti = np.array(Tnscan.h5out['factors/fac_Ti'])[0]
        #fac_n = np.array(Tnscan.h5out['factors/fac_n'])[0]
        #fac_vz = np.array(Tnscan.h5out['factors/fac_vz'])[0]

        self.scanlist = []
        for nfac in self.fac_n:
            for Tefac in self.fac_Te:
                for Tifac in self.fac_Ti:
                    for vzfac in self.fac_vz:
                        nform = "{:0.3f}".format(nfac)
                        if nform[0] == '0':
                            nform = nform[1:]
                        Teform = "{:0.3f}".format(Tefac)
                        if Teform[0] == '0':
                            Teform = Teform[1:]
                        Tiform = "{:0.3f}".format(Tifac)
                        if Tiform[0] == '0':
                            Tiform = Tiform[1:]
                        vzform = "{:0.3f}".format(vzfac)
                        if vzform[0] == '0':
                            vzform = vzform[1:]

                        self.scanlist.append('n'+nform+'Te'+Teform+'Ti'+Tiform+'vz'+vzform)
        return self.scanlist


    def get_mode_names(self, scanname):
        """Get the mode names, e.g. f_5_2, for a specific scanname. If no parameter scan was done, use '/' to indicate the root group."""
        self.modes_list = []
        for key in self.h5out[scanname].keys():
            if re.search('\Af_', key):
                self.modes_list.append(key)
        return self.modes_list


    def group_string(self, scanid, group, mode='f_5_2', dataset = ''):
        """Concatenate the strings of the scanid, mode, an arbitrary group and a specific dataset."""
        if scanid:
            scanid = scanid + '/'
        return scanid + mode + '/' + group + '/' + dataset


    def deql22(self, scanid, modes, time = 0):
        """Returns the profile of deql22"""
        return self.h5out[scanid][modes]['fort.5000'][str(5000+time)][3]

    def tor_resc(self, dqle22, mode='f_5_2'):
        """Toroidal rescaling of the diffusion coefficient with scalefactor from input hdf5"""
        self.check_mode(mode)
        self.get_mode_names('/')
        self.ind = self.modes_list.index(mode)
        self.scalefactor_sq = np.array(self.h5inp['/output/scalefactors_sq'])[self.ind]
        return dqle22*self.scalefactor_sq


    def get_deq22_res(self,scanid = '', time=0):
        """ Get De_22 at the resonant surfaces. This is done with the profile of deq22. Since its value at the resonant surface is now saved separately, this method is not used anymore."""
        self.scanid = scanid
        self.deq22_res_vec = np.empty([0])
        # iterate over given m values
        if (self.scanid == ''):
             for i in range(0,len(self.m)):
                self.modestring = 'f_'+ str(int(self.m[i])) + '_'+ str(int(self.n[i]))
                # get r and De22 values for given mode number and scanid (=factor setting)
                self.r = self.h5out[self.modestring]['fort.5000'][str(5000+time)]['r']
                self.deq22 = self.h5out[self.modestring]['fort.5000'][str(5000+time)]['dqle22']
                self.deq22 = self.deq22 * self.scalefactors_sq[i]
                self.fun = CubicSpline(self.r, self.deq22)
                self.deq22_res_vec = np.append(self.deq22_res_vec, self.fun(self.r_res[i]))

        else:
            for i in range(0,len(self.m)):
                self.modestring = 'f_'+ str(int(self.m[i][0])) + '_'+ str(int(self.n[i][0]))
                # get r and De22 values for given mode number and scanid (=factor setting)
                self.r = self.h5out[self.scanid][self.modestring]['fort.5000'][str(5000+time)]['r']
                self.deq22 = self.h5out[self.scanid][self.modestring]['fort.5000'][str(5000+time)]['dqle22']
                self.deq22 = self.deq22 * self.scalefactors_sq[i]
                self.fun = CubicSpline(self.r, self.deq22)
                self.deq22_res_vec = np.append(self.deq22_res_vec, self.fun(self.r_res[i]))
        return self.deq22_res_vec

    def get_deq22_res_read(self,scanid='/', time=0):
        self.scanid = scanid
        self.deq22_res_vec = np.empty([0])
        if (self.scanid=='/'):
            for i in range (0,len(self.m)):
                self.modestring = 'f_'+ str(int(self.m[i])) + '_'+ str(int(self.n[i]))
                self.deq22_res_vec = np.append(self.deq22_res_vec, np.array(self.h5out[self.scanid][self.modestring]['dqle22_res'])[0])
        return self.deq22_res_vec


    def nscan_plotdeq22_res(self,time=0, figuresize=(8,5)):
        """ Plot De_22_res """
        self.get_scan_names()
        self.deq22_mat = np.empty((0,3))
        # get all deq22_res values
        for scan in self.scans_list:
            self.deq22_mat = np.append(self.deq22_mat, [self.get_deq22_res(scan,time)], axis=0)
        self.mode_num = self.deq22_mat.shape[1] # get the number of modes from the dimension of the deq22 matrix
        fig, axs = plt.subplots(self.mode_num, sharex=True, figsize=figuresize)
        for i in range(0,self.mode_num):
            axs[i].set_yscale('log')
            axs[i].set_xscale('log')
            axs[i].set_xlabel('n factor/1')
            axs[i].set_ylabel('$D^{QL}_{e22}$/$D_a$')
            axs[i].axhline(y=1, color='grey')
            axs[i].grid(True, which="both")
            axs[i].scatter(self.fac_n[0], self.deq22_mat[:,i]/self.da_res[i])
        plt.tight_layout()


    def nscan_plotdeq22_res_interp(self, time=0, figuresize=(8,5)):
        """ Plot De_22_res """
        self.get_scan_names()
        self.deq22_mat = np.empty((0,3))
        # get all deq22_res values
        for scan in self.scans_list:
            self.deq22_mat = np.append(self.deq22_mat, [self.get_deq22_res(scan,time)], axis=0)
        self.mode_num = self.deq22_mat.shape[1] # get the number of modes from the dimension of the deq22 matrix

        # linear interpolation of logarithmic data interpolation
        # i.e. y = k*x + d, where y = log(Deq22) and x = log(n_fac)
        self.k = np.empty((self.mode_num))
        self.d = np.empty((self.mode_num))
        self.interp = []    # list that contains the interpolation functions
        self.x = np.linspace(np.log(self.fac_n[0][0]), np.log(self.fac_n[0][-1]), 100)

        fig, axs = plt.subplots(self.mode_num, sharex=True, figsize=figuresize)
        for i in range(0,self.mode_num):
            # create the interpolations
            self.interp.append(CubicSpline(np.log(self.fac_n[0]), np.log(self.deq22_mat[:,i])))
            self.k[i] = np.average(np.gradient(self.interp[i](self.x), self.x))
            self.d[i] = np.average(self.interp[i](self.x) - self.k[i] * self.x)
            axs[i].set_yscale('log')
            axs[i].set_xscale('log')
            axs[i].xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))

            axs[i].set_ylabel('$D^{QL}_{e22}$')
            axs[i].axhline(y=self.da_res[i], color='grey', linestyle='dashed',label='Threshold: $D_a^{res}$')
            axs[i].axvline(x=1, color='black', linestyle='dotted', label='Unscaled')
            #axs[i].grid(True, which="both")
            axs[i].scatter(self.fac_n[0], self.deq22_mat[:,i])
            axs[i].plot(np.exp(self.x), np.exp(self.k[i]*self.x + self.d[i]),label='$D^{QL}_{e22}$ ~$n^{'+"{:5.3f}".format(self.k[i])+'}$')
            axs[i].legend()

        axs[self.mode_num-1].set_xlabel('n factor / 1')
        plt.tight_layout()


    def plotprofiles(self, mode, time=0, scan=''):
        """ WIP!!!"""
        self.get_scan_names() # saves possible scan names in list
        self.profs = {'r' : np.empty((0)), 'n' : np.empty((0)), 'vz' : np.empty((0)), 'Te' : np.empty((0)),
                      'Ti' : np.empty((0)), 'Er' : np.empty((0)), 'sqrtg_Btheta_over_c' : np.empty((0))}
        if (scan == ''):
            count = 0
            for key in self.profs.keys():
                self.profs[key] = self.h5out[mode]['fort1000'][str(1000+time)][count]
                count = count + 1
        else:
            pass
            ### here the case if parameterscan was done has to be covered

    def getbalanceprof(self, time):
        """ Get the profiles from the balance code run (contained in fort.1000 group) after the specified time """
        # order of f.1000 data
        # r, n, vz, Te, Ti, Er, sqrtg_Btheta_over_c
        pass

    def closefile(self):
        """ Close the files. Necessary to free the file for other applications. """
        self.h5out.close()
        self.h5inp.close()


    def plt_simple_criterion(self, scanid = '/', time=0):
        """ Create a plot that illustrates the simple bifurcation criterion, i.e. the values of dqle22 at the resonant surfaces. Also plots the plasma pressure as a reference."""
        self.da_res = np.array(self.h5inp['output/Da_res']).transpose()[0]
        #if (scanid == ''):
            #self.deq22_res = self.get_deq22_res(scanid, time)
        if (scanid == 'all'):
            self.get_scan_names()

            for scan in self.scans_list:
                self.calc_pres(scan)
                self.deq22_res = self.get_deq22_res(scan, time)
                self.bifurcfactors = self.deq22_res/self.da_res
                fig, ax = plt.subplots(figsize=(10,4))
                ax_twin = plt.twinx()
                ax.set_title(scan)
                ax.plot(self.r_res.transpose()[0], self.bifurcfactors,marker='o', label='bifurcation factors')
                ax.axhline(y=1.0,linestyle = ':', color = 'grey' ,label='threshold')
                ax.set_yscale('log')
                ax.set_ylabel('$D^{QL}_{e22}/D_a$ @ res')
                ax.set_xlabel('r/cm')
                ax.set_xlim((40,70))
                ax_twin.plot(self.prof_rc, self.prof_p*2, color='orange', label ='pressure')
                ax_twin.set_ylabel('$p_{tot} / Pa', color='orange')
                #ax_twin.set_ylim((0,12e4));

                for elem in self.h5inp['output/zero_veperp']:
                    if (not np.isnan(elem)):
                        ax.axvline(x=elem, color = 'grey', linestyle = '--', label='Fluid resonance')

                for elem in self.h5inp['output/zero_vExB']:
                    if (not np.isnan(elem)):
                        ax.axvline(x=elem, color = 'grey', linestyle = '-.', label='ExB resonance')
                ax.legend()
                #plt.show()
            return list(map(plt.figure, plt.get_fignums()))

        if (scanid == '/'):
            self.calc_pres(scanid=scanid)
            self.deq22_res = self.get_deq22_res_read(scanid, time)
            self.bifurcfactors = self.deq22_res/self.da_res
            fig, ax = plt.subplots(figsize=(10,4))
            ax_twin = plt.twinx()
            ax.set_title(scanid)
            ax.plot(self.r_res.transpose()[0], self.bifurcfactors,marker='o', label='bifurcation factors')
            ax.axhline(y=1.0,linestyle = ':', color = 'grey' ,label='threshold')
            ax.set_yscale('log')
            ax.set_ylabel('$D^{QL}_{e22}/D_a$ @ res')
            ax.set_xlabel('r/cm')
            ax.set_xlim((40,70))
            ax_twin.plot(self.prof_rc, self.prof_p*2, color='orange', label ='pressure')
            ax_twin.set_ylabel('$p_{tot}$ / Pa', color='orange')
            ax_twin.set_ylim((0,14e4))

            ax_twin.ticklabel_format(useMathText=True)

            fluid_resonance = np.amax(np.array(self.h5inp['output/zero_veperp']))
            ExB_resonance = np.amax(np.array(self.h5inp['output/zero_vExB']))

            ax.axvline(x=fluid_resonance, color = 'grey', linestyle = '--', label='Fluid resonance')
            ax.axvline(x=ExB_resonance, color = 'grey', linestyle = '-.', label='ExB resonance')

            #for elem in self.h5inp['output/zero_veperp']:
            #    if (not np.isnan(elem)):
            #        ax.axvline(x=elem, color = 'grey', linestyle = '--', label='Fluid resonance')

            #for elem in self.h5inp['output/zero_vExB']:
            #    if (not np.isnan(elem)):
            #        ax.axvline(x=elem, color = 'grey', linestyle = '-.', label='ExB resonance')
            ax.legend()
            #plt.show()
            return list(map(plt.figure, plt.get_fignums()))

        else:
            self.calc_pres(scanid=scanid)
            self.deq22_res = self.get_deq22_res(scanid, time)
            self.bifurcfactors = self.deq22_res/self.da_res
            fig, ax = plt.subplots(figsize=(10,4))
            ax_twin = plt.twinx()
            ax.set_title(scanid)
            ax.plot(self.r_res.transpose()[0], self.bifurcfactors,marker='o', label='bifurcation factors')
            ax.axhline(y=1.0,linestyle = ':', color = 'grey' ,label='threshold')
            ax.set_yscale('log')
            ax.set_ylabel('$D^{QL}_{e22}/D_a$ @ res')
            ax.set_xlabel('r/cm')
            ax.set_xlim((40,70))
            ax_twin.plot(self.prof_rc, self.prof_p*2, color='orange', label ='pressure')
            ax_twin.set_ylabel('$p_{tot}$ / Pa', color='orange')
            ax_twin.set_ylim((0,14e4))

            ax_twin.ticklabel_format(useMathText=True)

            for elem in self.h5inp['output/zero_veperp']:
                if (not np.isnan(elem)):
                    ax.axvline(x=elem, color = 'grey', linestyle = '--', label='Fluid resonance')

            for elem in self.h5inp['output/zero_vExB']:
                if (not np.isnan(elem)):
                    ax.axvline(x=elem, color = 'grey', linestyle = '-.', label='ExB resonance')
            ax.legend()
            #plt.show()
            return list(map(plt.figure, plt.get_fignums()))


    def load_profile(self, proftype='', scanid='', mode='f_5_2', time=0):
        """ Load specific profiles from hdf5 file.
            Possible proftypes: n, Te, Ti, Er, vz, all """
        self.prof_rc = np.array(self.h5out[self.group_string(scanid, mode=mode, group ='/fort.1000/'+str(1000+time)+'/rc')])
        if (proftype==''):
            self.prof_n = np.array(self.h5out[self.group_string(scanid, mode=mode, group= '/fort.1000/'+str(1000+time)+'/n')])
            self.prof_Te = np.array(self.h5out[self.group_string(scanid, mode = mode, group ='/fort.1000/'+str(1000+time)+'/Te')])
            self.prof_Ti = np.array(self.h5out[self.group_string(scanid, mode = mode, group = '/fort.1000/'+str(1000+time)+'/Ti')])
            self.prof_Er = np.array(self.h5out[self.group_string(scanid, mode = mode, group = '/fort.1000/'+str(1000+time)+'/Er')])
            self.prof_vz = np.array(self.h5out[self.group_string(scanid, mode = mode, group = '/fort.1000/'+str(1000+time)+'/Vz')])

        elif (proftype=='n'):
            self.prof_n = np.array(self.h5out[self.group_string(scanid, mode = mode, group= '/fort.1000/'+str(1000+time)+'/n')])
            return self.prof_n
        elif (proftype=='Te'):
            self.prof_Te = np.array(self.h5out[self.group_string(scanid, mode = mode, group ='/fort.1000/'+str(1000+time)+'/Te')])
            return self.prof_Te
        elif (proftype=='Ti'):
            self.prof_Ti = np.array(self.h5out[self.group_string(scanid, mode = mode, group = '/fort.1000/'+str(1000+time)+'/Ti')])
            return self.prof_Ti
        elif (proftype=='Er'):
            self.prof_Er = np.array(self.h5out[self.group_string(scanid, mode = mode, group = '/fort.1000/'+str(1000+time)+'/Er')])[0:-1]
            return self.prof_Er
        elif (proftype=='Vz'):
            self.prof_vz = np.array(self.h5out[self.group_string(scanid, mode = mode, group = '/fort.1000/'+str(1000+time)+'/Vz')])
            return self.prof_vz

        else:
            print('Type not supported')



    def calc_pres(self, scanid='', time=0):
        """ Calculates the pressure from the density, electron temperature and constants. """
        k_B = 1.3807e-16
        EVK = 1.1604e4
        self.load_profile(proftype='n', scanid=scanid)
        self.load_profile(proftype='Te', scanid=scanid)
        self.prof_p = self.prof_n * self.prof_Te * k_B * EVK


    def calc_total_pres(self, scanid='', time=0):
        """ Calculate the total pressure, i.e. p_i + p_e."""
        k_B = 1.3807e-16
        EVK = 1.1604e4
        self.load_profile(proftype='n', scanid=scanid)
        self.load_profile(proftype='Te', scanid=scanid)
        self.load_profile(proftype='Ti', scanid=scanid)
        self.prof_p_tot = self.prof_n * (self.prof_Te+self.prof_Ti) * k_B * EVK

######################################################
# Time evolution specific

    def plt_Br_abs_ant(self, scanid ='', mode='f_5_2'):
        """ Plot the radial magnetic field at the antenna over time."""
        iterates = list(self.h5out[self.group_string(scanid=scanid, group='/fort.5000')].keys())
        self.Br_abs_ant = np.empty((0))

        for it in iterates:
            self.Br_abs_ant = np.append(self.Br_abs_ant, self.h5out[scanid+mode+ '/fort.5000/'+it+'/Br_abs'][-1])
        self.Br_ant_time = np.array(self.h5out[scanid+mode+'/timstep_evol.dat']).transpose()[5]
        self.Br_ant_timesel = np.append(self.Br_ant_time[0::10], self.Br_ant_time[-1])
        self.Br_abs_ant_inter = CubicSpline(self.Br_ant_timesel, self.Br_abs_ant)
        self.Br_abs_ant_new = self.Br_abs_ant_inter(self.Br_ant_time)
        self.antenna_factor = np.array(self.h5out[scanid+mode+'/br_abs_res.dat']).transpose()[2]
        fig = plt.figure()
        plt.plot(self.Br_ant_time, self.Br_abs_ant_new*np.sqrt(self.antenna_factor[0:-1]), marker='o')
        plt.xlabel('t/s')
        plt.ylabel('$|B_r^{ant}|C_{52}$')
        return list(map(plt.figure, plt.get_fignums()))


    def plt_Br_abs_res(self, scanid='', mode='f_5_2', diag=False):
        """ Plot the radial magnetic field at the resonant surface for a certain mode over time."""
        #timedataindex = 2
        #if diag:
        #    timedataindex = 5

        self.Br_abs_res = np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_res.dat')]).transpose()[3]
        #self.Br_abs_res_time = np.array(self.h5out[self.group_string(scanid,mode=mode, group='/timstep_evol.dat')]).transpose()[timedataindex]
        self.Br_abs_res_time = np.array(self.h5out[self.group_string(scanid, mode=mode, group='br_abs_res.dat')]).transpose()[1]
        fig = plt.figure()
        plt.plot(self.Br_abs_res_time, self.Br_abs_res)
        plt.xlabel('t/s')
        plt.ylabel('$|B_r^{res}|C_{52}$')
        return list(map(plt.figure, plt.get_fignums()))

    def plt_Br_over_Br0_abs_res(self, scanid ='', mode = 'f_5_2'):
        """ Plot shielding factor at resonant surface over time."""
        self.brvac = np.array(self.h5out[self.group_string(scanid, mode=mode, group ='Brvac.dat')][1])
        self.rr = np.array(self.h5out[self.group_string(scanid, group= 'fort.5000/5000/r')])

        keys = list(self.h5out[self.group_string(scanid, group='fort.5000')].keys())
        self.br_frac = np.empty((0))
        for key in keys:
            fac = np.array(self.h5out[self.group_string(scanid,group='fort.5000/')+key+'/Br_abs'])/self.brvac
            ffun = CubicSpline(self.rr, fac)
            self.br_frac = np.append(self.br_frac, ffun(self.r_res[1]))

        self.Br_frac_time = np.array(self.h5out[self.group_string(scanid,group='/timstep_evol.dat')]).transpose()[5]
        self.Br_frac_timesel = np.append(self.Br_frac_time[0::10], self.Br_frac_time[-1])
        self.Br_frac_interp_fun = CubicSpline(self.Br_frac_timesel, self.br_frac)
        self.Br_frac_interp = self.Br_frac_interp_fun(self.Br_frac_time)

        fig = plt.figure()
        plt.scatter(self.Br_frac_timesel, self.br_frac,marker='o', label = 'Data')

        plt.xlabel('t/s')
        plt.ylabel('$|B_r^{res}(t)/B_0^{res}|$')
        plt.plot(self.Br_frac_time[100:], self.Br_frac_interp[100:], color='orange', label='Interpolated')
        plt.legend()
        return list(map(plt.figure, plt.get_fignums()))

    def plt_shielding_fac_wo_diag(self, scanid='', mode = 'f_5_2', save=False, title=False, out_type='pdf'):
        """ Plot shielding factor at resonant surface over time, when the balance run was done without the diagnostics output."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        self.brvac = np.array(self.h5out[self.group_string(scanid, mode=mode, group ='Brvac_res')])
        self.Br_abs_res = np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_res.dat')]).transpose()[3]
        self.Br_abs_res_time = np.array(self.h5out[self.group_string(scanid, mode=mode, group='br_abs_res.dat')]).transpose()[1]
        self.antenna_factor = np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_res.dat')]).transpose()[2]

        fig = plt.figure(figsize=(7.5,5))
        plt.tight_layout()
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = ' + mode )
        plt.plot(self.Br_abs_res_time, self.Br_abs_res/self.brvac/np.sqrt(self.antenna_factor))
        plt.xlabel('t/s')
        plt.ylabel('$|B_r^{res}(t)/B_0^{res}|$')
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_Br_over_Br0_abs' + mode+ '.'+ out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))



    def plt_antenna_ramp_up(self, scanid='', save=False, title=False, out_type='pdf', plt_type='plot', perc=True):
        """ Plot antenna factor over time for the case when the diagnostics output was turned off. If perc is True, it plots it as percentage of the experimental value."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')
        fig = plt.figure(figsize=(7.5,5))
        plt.rc('font', size=18)
        plt.tight_layout()
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms.')

        mode_names = self.get_mode_names('/'+scanid)

        for mode in mode_names:

            self.Br_abs_res_time = np.array(self.h5out[self.group_string(scanid, mode=mode, group='br_abs_time')])
            self.antenna_factor = np.sqrt(np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_antenna_factor')]))
            if perc:
                self.antenna_factor = 100*self.antenna_factor/np.sqrt(self.scalefactors_sq[1])

            if plt_type == 'scatter':
                plt.scatter(self.Br_abs_res_time, self.antenna_factor, label='m = ' + mode[2] + '; n = '+ mode[4])
            if plt_type == 'plot':
                plt.plot(self.Br_abs_res_time, self.antenna_factor, label='m = ' + mode[2] + '; n = '+ mode[4],lw=4)
            else:
                SystemError('Plot type not supported')

        plt.xlabel(r'$t \; / \; s$')
        if perc:
            plt.ylabel('Antenna factor / %')
        else:
            plt.ylabel('Antenna factor / 1')
        plt.legend()
        plt.grid()
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_antenna_factor.'+ out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))

    def plt_shielding_fac_over_ant_fac(self, scanid='', save=False, title=False, out_type = 'pdf', plt_type = 'plot', perc = True):
        """ Plot shielding factor over antenna factor for the case when the diagnostics output was turned off. If perc is True, it plots the antenna factor as percentage of the experimental value."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')
        fig = plt.figure(figsize=(7.5,5))
        plt.rc('font', size=18)
        plt.tight_layout()
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms.')

        mode_names = self.get_mode_names('/'+scanid)

        for mode in mode_names:
            self.brvac = np.array(self.h5out[self.group_string(scanid, mode=mode, group ='Brvac_res')])
            self.Br_abs_res = np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_res')])
            self.Br_abs_res_time = np.array(self.h5out[self.group_string(scanid, mode=mode, group='br_abs_time')])
            self.antenna_factor = np.sqrt(np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_antenna_factor')]))

            # If true, scale antenna factor with antenna factor max
            if perc:
                #max_ant = np.sqrt(np.array(self.h5inp['/output/scalefactors_sq'])[1])
                self.antenna_factor_x = 100*self.antenna_factor/np.sqrt(self.scalefactors_sq[1])

            if plt_type == 'scatter':
                plt.scatter(self.antenna_factor_x, self.Br_abs_res/self.brvac/self.antenna_factor, label='m = ' + mode[2] + '; n = '+ mode[4])
            if plt_type == 'plot':
                plt.plot(self.antenna_factor_x, self.Br_abs_res/self.brvac/self.antenna_factor, label='m = ' + mode[2] + '; n = '+ mode[4],lw=4)
            else:
                SystemError('Plot type not supported')

        plt.xlabel(r'antenna factor / %')
        plt.ylabel('$|B_r^{res}(t)/B_0^{res}|$')
        plt.legend()
        plt.grid()
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_shielding_fac_over_ant_fac.'+ out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))


    def plt_shielding_fac_wo_diag_all_modes_one_fig(self, scanid='', save=False, title=False, out_type='pdf', plt_type='scatter'):
        """ Plot Br over Br_vas absolute at the resonant surface, for the case when the diagnostics output was turned off.
        Plots the values for all modes in one plot."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        fig = plt.figure(figsize=(7.5,5))
        plt.rc('font', size=18)
        plt.tight_layout()
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms.')

        mode_names = self.get_mode_names('/'+scanid)

        for mode in mode_names:

            self.brvac = np.array(self.h5out[self.group_string(scanid, mode=mode, group ='Brvac_res')])
            self.Br_abs_res = np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_res')])
            self.Br_abs_res_time = np.array(self.h5out[self.group_string(scanid, mode=mode, group='br_abs_time')])
            self.antenna_factor = np.sqrt(np.array(self.h5out[self.group_string(scanid,mode=mode, group='br_abs_antenna_factor')]))

            if plt_type == 'scatter':
                plt.scatter(self.Br_abs_res_time, self.Br_abs_res/self.brvac/self.antenna_factor, label='m = ' + mode[2] + '; n = '+ mode[4])
            if plt_type == 'plot':
                plt.plot(self.Br_abs_res_time, self.Br_abs_res/self.brvac/self.antenna_factor, label='m = ' + mode[2] + '; n = '+ mode[4],lw=4)
            else:
                SystemError('Plot type not supported')

        plt.xlabel('t/s')
        plt.ylabel('$|B_r^{res}(t)/B_0^{res}|$')
        plt.legend()
        plt.grid()
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_shielding_fac_t_all_modes.'+ out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))




    def plt_n_t(self, time='last', scanid='', mode='f_5_2', save=False, title=False, out_type='pdf', y_range=0.1, r_range=0.5):
        """ Plot the density profile for t = 0 and arbitrary second time.
        The second time is defaulted to the last time step.
        y_range and r_range are used to control the zoom window. The former is given in percent (of the value of the profile at the resonant surface) and the latter in cm."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        mode = self.check_mode(mode)

        if time =='last':
            # get the last entry of the groups, i.e. the last time step
            self.f1000tlist = list(self.h5out[self.group_string(scanid, group='fort.1000', mode=mode)].keys())
            time = self.f1000tlist[-1]
            time = int(time) - 1000
            #print(time)

        fig = plt.figure(figsize=(15,5))
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Density time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = ' + mode )
        plt.tight_layout()
        plt.subplot(121)
        self.load_profile(proftype='n', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_n, label='t=0')
        self.load_profile(proftype='n', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_n, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('n/$cm^{-3}$')
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        # zoom in on the resonance
        plt.subplot(122)
        self.load_profile(proftype='n', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_n, label='t=0s')
        interp_n_t0 = CubicSpline(self.prof_rc, self.prof_n)
        self.load_profile(proftype='n', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_n, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('n/$cm^{-3}$')
        plt.xlabel('r/cm')
        plt.legend()
        # set x and y limits by averaging the values of the lines
        plt.xlim(left=self.r_res[0]- self.prof_t_plot_r_offset, right = self.r_res[0] + self.prof_t_plot_r_offset)
        interp_n = CubicSpline(self.prof_rc, self.prof_n)
        #avg_n_lower = (interp_n(self.r_res[0] + self.prof_t_plot_r_offset) + interp_n_t0(self.r_res[0] + self.prof_t_plot_r_offset))/2
        #avg_n_upper = (interp_n(self.r_res[0] - self.prof_t_plot_r_offset) + interp_n_t0(self.r_res[0] - self.prof_t_plot_r_offset))/2

        plt.ylim((interp_n(self.r_res[0]) *(1- y_range), interp_n(self.r_res[0]) *(1+ y_range)))
        #plt.ylim((interp_n(self.r_res[0]+0.5), interp_n(self.r_res[0]-0.5)))
        plt.axvline(self.r_res[0], c='grey', ls=':')
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_'+mode+ \
                '_nprof_t.' + out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))


    def plt_Te_t(self, time='last', scanid='', mode='f_5_2', save=False, title=False, out_type ='pdf'):
        """ Plot the electron temperature profile for t = 0 and arbitrary
        second time. The second time is defaulted to the last time step."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        mode = self.check_mode(mode)

        # define a scale for the plot
        keVscale = 1e3

        if time =='last':
            # get the last entry of the groups, i.e. the last time step
            self.f1000tlist = list(self.h5out[self.group_string(scanid, group='fort.1000', mode=mode)].keys())
            time = self.f1000tlist[-1]
            time = int(time) - 1000
        fig = plt.figure(figsize=(15,5))
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Electron temperature time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = ' + mode )
        plt.tight_layout()
        plt.subplot(121)
        self.load_profile(proftype='Te', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_Te / keVscale, label='t=0')
        self.load_profile(proftype='Te', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_Te / keVscale, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('Te/$keV$')
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        # zoom in on the resonance
        plt.subplot(122)
        self.load_profile(proftype='Te', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_Te / keVscale, label='t=0s')
        interp_Te_t0 = CubicSpline(self.prof_rc, self.prof_Te)
        self.load_profile(proftype='Te', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_Te / keVscale, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('Te/$keV$')
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        plt.xlim(left=self.r_res[0] - self.prof_t_plot_r_offset, right = self.r_res[0] + self.prof_t_plot_r_offset)
        interp_Te = CubicSpline(self.prof_rc, self.prof_Te)
        avg_n_lower = (interp_Te(self.r_res[0] + self.prof_t_plot_r_offset) + interp_Te_t0(self.r_res[0] + self.prof_t_plot_r_offset))/2 / keVscale
        avg_n_upper = (interp_Te(self.r_res[0] - self.prof_t_plot_r_offset) + interp_Te_t0(self.r_res[0] - self.prof_t_plot_r_offset))/2 / keVscale

        plt.ylim((avg_n_lower * self.prof_t_plot_y_offset_lower, avg_n_upper * self.prof_t_plot_y_offset_upper))
        #plt.ylim((interp_Te(self.r_res[0]+0.5) / keVscale, interp_Te(self.r_res[0]-0.5) / keVscale))
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_'+mode+ \
                '_Teprof_t.'+ out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))

    def plt_Ti_t(self, time='last', scanid='', mode='f_5_2', save=False, title=False, out_type='pdf'):
        """ Plot the ion temperature profile for t = 0 and arbitrary
        second time. The second time is defaulted to the last time step."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        mode = self.check_mode(mode)

        # define a scale for the plot
        keVscale = 1e3

        if time =='last':
            # get the last entry of the groups, i.e. the last time step
            self.f1000tlist = list(self.h5out[self.group_string(scanid, group='fort.1000', mode=mode)].keys())
            time = self.f1000tlist[-1]
            time = int(time) - 1000
        fig = plt.figure(figsize=(15,5))
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Ion temperature time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = ' + mode )
        plt.tight_layout()
        plt.subplot(121)
        self.load_profile(proftype='Ti', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_Ti / keVscale, label='t=0')
        self.load_profile(proftype='Ti', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_Ti / keVscale, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('Ti/$keV$')
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        # zoom in on the resonance
        plt.subplot(122)
        self.load_profile(proftype='Ti', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_Ti / keVscale, label='t=0s')
        interp_Ti_t0 = CubicSpline(self.prof_rc, self.prof_Ti)
        self.load_profile(proftype='Ti', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_Ti / keVscale, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('Ti/$keV$')
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        plt.xlim(left=self.r_res[0] - self.prof_t_plot_r_offset, right = self.r_res[0] + self.prof_t_plot_r_offset)
        interp_Ti = CubicSpline(self.prof_rc, self.prof_Ti)
        avg_n_lower = (interp_Ti(self.r_res[0] + self.prof_t_plot_r_offset) + interp_Ti_t0(self.r_res[0] + self.prof_t_plot_r_offset))/2 / keVscale
        avg_n_upper = (interp_Ti(self.r_res[0] - self.prof_t_plot_r_offset) + interp_Ti_t0(self.r_res[0] - self.prof_t_plot_r_offset))/2 / keVscale

        plt.ylim((avg_n_lower * self.prof_t_plot_y_offset_lower, avg_n_upper * self.prof_t_plot_y_offset_upper))
        #plt.ylim((interp_Ti(self.r_res[0]+0.5) / keVscale, interp_Ti(self.r_res[0]-0.5) / keVscale))
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_'+mode+ \
                '_Tiprof_t.'+ out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))


    def plt_vz_t(self, time='last', scanid='', mode='f_5_2', save=False, title=False, out_type='pdf'):
        """ Plot the toroidal velocity profile for t = 0 and arbitrary
        second time. The second time is defaulted to the last time step."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        mode = self.check_mode(mode)

        if time =='last':
            # get the last entry of the groups, i.e. the last time step
            self.f1000tlist = list(self.h5out[self.group_string(scanid, group='fort.1000', mode=mode)].keys())
            time = self.f1000tlist[-1]
            time = int(time) - 1000
        fig = plt.figure(figsize=(15,5))
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Toroidal velocity time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = ' + mode )
        plt.tight_layout()
        plt.subplot(121)
        self.load_profile(proftype='Vz', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_vz, label='t=0')
        self.load_profile(proftype='Vz', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_vz, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('vz/$cm/s$')
        plt.xlabel('r/cm')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        # zoom in on the resonance
        plt.subplot(122)
        self.load_profile(proftype='Vz', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_vz, label='t=0s')
        interp_vz_t0 = CubicSpline(self.prof_rc, self.prof_vz)
        self.load_profile(proftype='Vz', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_vz, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('vz/$cm/s$')
        plt.xlabel('r/cm')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        plt.xlim(left=self.r_res[0] - self.prof_t_plot_r_offset, right = self.r_res[0] + self.prof_t_plot_r_offset)
        interp_vz = CubicSpline(self.prof_rc, self.prof_vz)
        avg_n_lower = (interp_vz(self.r_res[0] + self.prof_t_plot_r_offset) + interp_vz_t0(self.r_res[0] + self.prof_t_plot_r_offset))/2
        avg_n_upper = (interp_vz(self.r_res[0] - self.prof_t_plot_r_offset) + interp_vz_t0(self.r_res[0] - self.prof_t_plot_r_offset))/2

        plt.ylim((avg_n_lower * self.prof_t_plot_y_offset_lower, avg_n_upper * self.prof_t_plot_y_offset_upper))
        #plt.ylim((interp_vz(self.r_res[0]+0.5), interp_vz(self.r_res[0]-0.5)))
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_'+mode+ \
                '_vzprof_t.' + out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))


    def plt_Er_t(self, time='last', scanid='', mode='f_5_2', save=False, title=False, out_type='pdf'):
        """ Plot the radial electric field profile for t = 0 and arbitrary
        second time. The second time is defaulted to the last time step."""

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        mode = self.check_mode(mode)

        if time =='last':
            # get the last entry of the groups, i.e. the last time step
            self.f1000tlist = list(self.h5out[self.group_string(scanid, group='fort.1000', mode=mode)].keys())
            time = self.f1000tlist[-1]
            time = int(time) - 1000
        fig = plt.figure(figsize=(15,5))
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Radial electric field time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = '+mode )
        plt.tight_layout()
        plt.subplot(121)
        self.load_profile(proftype='Er', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_Er[0:-1], label='t=0')
        self.load_profile(proftype='Er', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_Er[0:-1], label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('Er/$statV cm^{-1}$')
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()
        # zoom in on the resonance
        plt.subplot(122)
        self.load_profile(proftype='Er', scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, self.prof_Er[0:-1], label='t=0s')
        interp_Er_t0 = CubicSpline(self.prof_rc, self.prof_Er[0:-1])
        self.load_profile(proftype='Er', scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, self.prof_Er[0:-1], label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel('Er/$statV cm^{-1}$')
        plt.xlabel('r/cm')
        plt.legend()
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.xlim(left=self.r_res[0] - self.prof_t_plot_r_offset, right = self.r_res[0] + self.prof_t_plot_r_offset)
        interp_Er = CubicSpline(self.prof_rc, self.prof_Er[0:-1])
        avg_n_lower = (interp_Er(self.r_res[0] + self.prof_t_plot_r_offset) + interp_Er_t0(self.r_res[0] + self.prof_t_plot_r_offset))/2
        avg_n_upper = (interp_Er(self.r_res[0] - self.prof_t_plot_r_offset) + interp_Er_t0(self.r_res[0] - self.prof_t_plot_r_offset))/2

        plt.ylim((avg_n_lower * self.prof_t_plot_y_offset_lower, avg_n_upper * self.prof_t_plot_y_offset_upper))
        #plt.ylim((interp_Er(self.r_res[0]+0.5), interp_Er(self.r_res[0]-0.5)))
        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_'+mode+ \
                '_Erprof_t.' + out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))


    def plt_all_profiles_t(self, time='last', scanid='', mode='f_5_2', save=False, title=False, out_type='pdf'):
        """ Plot all profiles over t. Contains zoomed in subplot."""

        mode = self.check_mode(mode)
        self.plt_Br_over_Br0_abs_res_wo_diag(scanid=scanid, mode=mode, save=save, title=title, out_type=out_type)
        self.plt_n_t(time=time, scanid=scanid, mode=mode, save=save, title=title, out_type=out_type)
        self.plt_Te_t(time=time, scanid=scanid, mode=mode, save=save, title=title, out_type=out_type)
        self.plt_Ti_t(time=time, scanid=scanid, mode=mode, save=save, title=title, out_type=out_type)
        self.plt_vz_t(time=time, scanid=scanid, mode=mode, save=save, title=title, out_type=out_type)
        self.plt_Er_t(time=time, scanid=scanid, mode=mode, save=save, title=title, out_type=out_type)

    def plt_all_profiles_t_all_modes(self, time='last', scanid='', save=False, title=False, out_type='pdf'):
        """ plt_all_profiles_t for all mode numbers."""
        for i in self.m:
            #print('mode = '+ 'f_'+"{:.0f}".format(i)+'_2')
            self.plt_all_profiles_t(time=time, scanid=scanid, mode='f_'+"{:.0f}".format(i)+'_2', save=save, title=title, out_type=out_type)


    def get_real_time(self, timenum, scanid = '', mode = 'f_5_2'):
        """Converts balance run step number to the corresponding real time value."""
        if type(timenum) == str:
            try:
                timenum = int(timenum)
            except:
                print('not a valid input')

        self.timedata = np.array(self.h5out[self.group_string(scanid, group= '/timstep_evol.dat', mode=mode)]).transpose()[2]
        if (timenum>999 and timenum < 4999):
            timenum = timenum - 1000
        elif (timenum>4999):
            timenum = timenum - 5000
        try:
            return self.timedata[timenum-2]
        except:
            print('Something went wrong')


    def get_ped_val(self):
        """ WIP!!!"""
        self.load_profile(proftype='n')
        #self.prof_n
        #self.prof_rc
        #interpfun = interpolate.interp1d(self.prof_rc, self.prof_n,kind='cubic')
        self.grad_n = np.gradient(self.prof_n, self.prof_rc)


    def plt_nTscan(self, mode='f_5_2'):
        """ Plot 2D parameter scan over density and electron temperature."""
        self.get_scan_names() # saves possible scan names in list
        Temesh, nmesh = np.meshgrid(self.fac_Te, self.fac_n)

        k_B = 1.3807e-16
        EVK = 1.1604e4
        self.prof_rc = np.array(self.h5out['init_params/r'])[0:-1]
        prof_n = np.array(self.h5out['init_params/n'])
        prof_Te = np.array(self.h5out['init_params/Te'])
        prof_p = prof_n * prof_Te * k_B * EVK

        # interpolate
        interp_grad_p = CubicSpline(self.prof_rc, np.gradient(prof_p))
        interp_grad_grad_p = CubicSpline(self.prof_rc, np.gradient(np.gradient(prof_p)))
        interp_p = CubicSpline(self.prof_rc, prof_p)


        #interp_n = CubicSpline(self.prof_rc, self.h5out['init_params/n'])
        interp_n = CubicSpline(self.prof_rc, prof_n)
        interp_grad_n = CubicSpline(self.prof_rc, np.gradient(prof_n))
        interp_Te = CubicSpline(self.prof_rc, np.array(prof_Te))#/(1.6022*10**-12))

        gradgradproots = interp_grad_grad_p.roots()
        gradproots = interp_grad_p.roots()
        gradnroots = interp_grad_n.roots()
        Teped = interp_Te(gradgradproots[2])
        nped = interp_n(gradgradproots[2])

        self.get_mode_names(self.scans_list[1])
        #modelist = ['/f_5_2/dqle22_res','/f_6_2/dqle22_res','/f_7_2/dqle22_res']
        #count1 = 0
        #count2 = 0
        #for scan in self.scans_list:
        #    for mode in self.modes_list:
        #        dqle22[count1, count2] = np.array(self.h5out[mode+'/dqle22_res'])
        #        count2 += 1
        #    count1 += 1
        #    count2 = 0
        dqle22 = np.array(self.h5out[mode+'/dqle22_res'])
        # toroidal rescaling of the diffusion coefficient
        dqle22 = self.tor_resc(dqle22, mode)

        self.da_res = np.array(self.h5inp['output/Da_res']).transpose()[0]

        nempval = 3.3e13
        tempval = 1e3

        nscaling = 1e13
        Tscaling = 1e3

        upper = 3.5
        lower = -4.5
        xlimright = 35
        nums = 15
        lvl = nums#np.linspace(lower,upper,nums)
        vcmap = 'coolwarm'


        fig = plt.figure(figsize=(16,8))
        # plotting
        for i in range(0,1):
            dql_re = np.transpose(np.reshape(dqle22/self.da_res[i],(np.size(self.fac_Te), np.size(self.fac_n))))

            plt.title(str(int(self.shot)) + ' @ '+ str(int(self.time))+ 'ms ; m = '+str(int(self.m[i])))
            # filled contour plot
            plot = plt.contourf(nmesh*nped/nscaling, Temesh*Teped/Tscaling, np.log10(dql_re), levels=lvl,cmap=vcmap)
            #plt.xlim(right=xlimright)

            # contours
            conpltgrey = plt.contour(nmesh*nped/nscaling, Temesh*Teped/Tscaling, np.log10(dql_re), linestyles='solid', levels=lvl)
            plt.clabel(conpltgrey, fmt='%2.1f', colors = 'k', fontsize=12)

            # threshold line
            conplt0 = plt.contour(nmesh*nped/nscaling, Temesh*Teped/Tscaling, np.log10(dql_re), levels=np.array([0.0]))
            plt.clabel(conplt0,fmt='%2.1f', colors='k', fontsize=12)
            h1,_ = conplt0.legend_elements() # handle for legend

            pltref = plt.scatter(1.0*nped/nscaling,1.0*Teped/Tscaling, color='k', marker='*')

            # color bar
            cbar = plt.colorbar(plot)
            cbar.set_label('$log_{10}(D^{ql}_{e22}/D_a)$')

            # calculated points
            pltcalc = plt.scatter(nmesh*nped/nscaling,Temesh*Teped/Tscaling, color = 'grey', marker = '.')

            #empirical bounds
            pltnemp = plt.axvline(x=nempval/nscaling, linestyle='--', color='w')
            #plt.text(nempval/nscaling + 0.2, plt.Axes.get_ylim()[1]*0.75 ,s='$n_{p}^{emp}$',color='w',rotation=90,fontsize=18)

            plttemp = plt.axhline(y=tempval/Tscaling, linestyle='--',color='w')
            #plt.text(plt.Axes.get_xlim()[1]*0.75, (tempval+300)/Tscaling,s='$T_{p}^{emp}$',color='w',fontsize=18)

            # legend
            plt.legend([h1[0], pltref,pltcalc] , ['Threshold','Reference point','Calculated values'])

            # label
            plt.xlabel('$n_{p} \, /10^{13} \, cm^{-3}$',fontsize = 14)
            plt.ylabel('$T_{p} \, /10^3 \, eV$', fontsize = 14)
            plt.show()
            #fig.savefig('plots/parscan_'+str(int(np.array(Tnscan.shot)[0][0])) + '.'+ str(int(np.array(Tnscan.time)[0][0]))+'_m'+str(int(Tnscan.m.transpose()[0][i]))+'.pdf')


    def plt_nTscan_all_modes(self, save=False, type_out='pdf', nums=15):
        """ Plot 2D parameter scan over density and electron temperature for all modes."""
        self.get_scan_names() # saves possible scan names in list
        Temesh, nmesh = np.meshgrid(self.fac_Te, self.fac_n)

        k_B = 1.3807e-16
        EVK = 1.1604e4
        self.prof_rc = np.array(self.h5out['init_params/r'])[0:-1]
        prof_n = np.array(self.h5out['init_params/n'])
        prof_Te = np.array(self.h5out['init_params/Te'])
        prof_p = prof_n * prof_Te * k_B * EVK

        # interpolate
        interp_grad_p = CubicSpline(self.prof_rc, np.gradient(prof_p))
        interp_grad_grad_p = CubicSpline(self.prof_rc, np.gradient(np.gradient(prof_p)))
        interp_p = CubicSpline(self.prof_rc, prof_p)


        #interp_n = CubicSpline(self.prof_rc, self.h5out['init_params/n'])
        interp_n = CubicSpline(self.prof_rc, prof_n)
        interp_grad_n = CubicSpline(self.prof_rc, np.gradient(prof_n))
        interp_Te = CubicSpline(self.prof_rc, np.array(prof_Te))#/(1.6022*10**-12))

        gradgradproots = interp_grad_grad_p.roots()
        gradproots = interp_grad_p.roots()
        gradnroots = interp_grad_n.roots()
        Teped = interp_Te(gradgradproots[2])
        nped = interp_n(gradgradproots[2])

        self.get_mode_names(self.scans_list[1])

        nempval = 3.3e13
        tempval = 1e3

        nscaling = 1e13
        Tscaling = 1e3

        upper = 3.2
        lower = -7
       # xlimright = 35
        lvl = np.linspace(lower,upper,nums)
        vcmap = 'coolwarm'

        i=0
        for mode in self.modes_list:
            dqle22 = np.array(self.h5out[mode+'/dqle22_res'])
            dqle22 = self.tor_resc(dqle22,mode)
            self.da_res = np.array(self.h5inp['output/Da_res']).transpose()[0]

            fig = plt.figure(figsize=(16,8))
            # plotting
            dql_re = np.reshape(dqle22/self.da_res[i],(np.size(self.fac_n), np.size(self.fac_Te))).transpose()

            plt.title(str(int(self.shot)) + ' @ '+ str(int(self.time))+ 'ms ; m = '+str(int(self.m[i])))
            # filled contour plot
            plot = plt.contourf(nmesh*nped/nscaling, Temesh*Teped/Tscaling, np.log10(dql_re), levels=lvl,cmap=vcmap)
            #plt.xlim(right=xlimright)

            # contours
            conpltgrey = plt.contour(nmesh*nped/nscaling, Temesh*Teped/Tscaling, np.log10(dql_re), linestyles='solid', levels=lvl)
            plt.clabel(conpltgrey, fmt='%2.1f', colors = 'k', fontsize=12)

            # threshold line
            conplt0 = plt.contour(nmesh*nped/nscaling, Temesh*Teped/Tscaling, np.log10(dql_re), levels=np.array([0.0]))
            plt.clabel(conplt0,fmt='%2.1f', colors='k', fontsize=12)
            h1,_ = conplt0.legend_elements() # handle for legend

            pltref = plt.scatter(1.0*nped/nscaling,1.0*Teped/Tscaling, color='k', marker='*')

            # color bar
            cbar = plt.colorbar(plot)
            cbar.set_label('$log_{10}(D^{ql}_{e22}/D_a)$')

            # calculated points
            pltcalc = plt.scatter(nmesh*nped/nscaling,Temesh*Teped/Tscaling, color = 'grey', marker = '.')

            #empirical bounds
            pltnemp = plt.axvline(x=nempval/nscaling, linestyle='--', color='w')
            #plt.text(nempval/nscaling + 0.2, plt.Axes.get_ylim()[1]*0.75 ,s='$n_{p}^{emp}$',color='w',rotation=90,fontsize=18)

            plttemp = plt.axhline(y=tempval/Tscaling, linestyle='--',color='w')
            #plt.text(plt.Axes.get_xlim()[1]*0.75, (tempval+300)/Tscaling,s='$T_{p}^{emp}$',color='w',fontsize=18)

            # legend
            plt.legend([h1[0], pltref,pltcalc] , ['Threshold','Reference point','Calculated values'])

            # label
            plt.xlabel('$n_{p} \, /10^{13} \, cm^{-3}$',fontsize = 14)
            plt.ylabel('$T_{p} \, /10^3 \, eV$', fontsize = 14)
            plt.show()


            if save:
                fig.savefig('plots/parscan_'+str(int(self.shot)) + '.'+ str(int(self.time))+'_m'+str(int(self.m[i]))+'.'+type_out, bbox_inches='tight', dpi=150)

            i = i+1

    def plt_profile(self, prof_type='n', time='last', mode='f_5_2', title=False, save=False, out_type='pdf', y_range=0.1, r_range=0.5, scanid=''):
        """ Plot a given profile for t = 0 and arbitrary second time.
        The second time is defaulted to the last time step.
        y_range and r_range are used to control the zoom window. The former is given in percent (of the value of the profile at the resonant surface) and the latter in cm."""
        if not (prof_type in self.profile_types.keys()):
            raise ValueError('Profile type does not exist. Supported types:'+ str(self.profile_types.keys()))

        mode = self.check_mode(mode)

        if not (out_type == 'pdf' or out_type == 'jpg'):
            raise ValueError('Wrong output type! Only jpg and pdf are supported.')

        if time =='last':
            # get the last entry of the groups, i.e. the last time step
            self.f1000tlist = list(self.h5out[self.group_string(scanid, group='fort.1000', mode=mode)].keys())
            time = self.f1000tlist[-1]
            time = int(time) - 1000
            #print(time)
        fig = plt.figure(figsize=(15,5))
        if title:
            plt.suptitle('ELM Suppression in Hydrogen. Density time evolution of shot ' + "{:.0f}".format(self.shot) \
                + ' at ' + "{:.0f}".format(self.time) + 'ms. Mode = ' + mode )
        plt.tight_layout()

        plt.subplot(121)
        profile = self.load_profile(proftype=prof_type, scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, profile, label='t=0')
        profile = self.load_profile(proftype=prof_type, scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, profile, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel(prof_type+self.profile_types[prof_type])
        plt.xlabel('r/cm')
        plt.axvline(self.r_res[0], c='grey', ls=':')
        plt.legend()

        # zoom in on the resonance
        plt.subplot(122)
        profile = self.load_profile(proftype=prof_type, scanid=scanid, mode=mode, time=0)
        plt.plot(self.prof_rc, profile, label='t=0s')
        interp_profile_t0 = CubicSpline(self.prof_rc, profile)
        profile = self.load_profile(proftype=prof_type, scanid=scanid, mode=mode, time=time)
        plt.plot(self.prof_rc, profile, label='t='+"{:5.3f}".format(self.get_real_time(time, mode=mode))+'s')
        plt.ylabel(prof_type+self.profile_types[prof_type])
        plt.xlabel('r/cm')
        plt.legend()

        plt.xlim(left=self.r_res[0]- r_range, right = self.r_res[0] + r_range)
        interp_profile = CubicSpline(self.prof_rc, profile)
        plt.ylim((interp_profile(self.r_res[0]) *(1- y_range), interp_profile(self.r_res[0]) *(1+ y_range)))
        # plot resonant line
        plt.axvline(self.r_res[0], c='grey', ls=':')

        if save==True:
            fig.savefig('plots/'+"{:.0f}".format(self.shot)+'_'+"{:.0f}".format(self.time)+'_'+mode+ \
                '_'+prof_type+'_prof_t.' + out_type, bbox_inches='tight', dpi=150)
        return list(map(plt.figure, plt.get_fignums()))

    def check_mode(self, mode):
        if type(mode) == int:
            mode = 'f_' + "{:.0f}".format(self.m[mode]) + '_' + "{:.0f}".format(self.n[mode])
        return mode
