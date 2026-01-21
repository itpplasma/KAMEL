import h5py
import numpy as np
import os

class Balance_Input_h5:
    """ Class to prepare the input h5 file for the balance code."""
    
    debug = False
    
    def __init__(self, h5_file_path, profile_path, debug=False):
        self.profile_path = profile_path
        self.debug = debug

    def get_required_data(self):
        self.n = np.loadtxt(self.profile_path + '/n.dat')
        self.Te = np.loadtxt(self.profile_path + '/Te.dat')
        self.Ti = np.loadtxt(self.profile_path + '/Ti.dat')
        self.Vz = np.loadtxt(self.profile_path + '/Vz.dat')
        self.Vth = np.loadtxt(self.profile_path + '/Vth.dat')
        self.Er = np.loadtxt(self.profile_path + '/Er.dat')
        self.q = np.loadtxt(self.profile_path + '/q.dat')

        if not os.path.isfile(self.profile_path + '/Da.dat'):
            self.Da = np.zeros_like(self.n)
            self.Da[:,0] = self.n[:,0]
            self.Da[:,1] = 1e4 # order of magnitude guess
    
    def write_data_to_h5(self, file_name, facs):
        if self.debug: print('Writing data to input h5 file: ', file_name)
        h5f = h5py.File(file_name, 'w')
        self.write_with_bound_info(h5f, '/da_estimation/Da', data=self.Da[:,1])
        self.write_with_bound_info(h5f, '/da_estimation/r', data=self.Da[:,0])

        self.write_with_bound_info(h5f, '/preprocprof/r_out', data=self.n[:,0])
        self.write_with_bound_info(h5f, '/preprocprof/n', data=self.n[:,1])
        self.write_with_bound_info(h5f, '/preprocprof/Te', data=self.Te[:,1])
        self.write_with_bound_info(h5f, '/preprocprof/Ti', data=self.Ti[:,1])
        self.write_with_bound_info(h5f, '/preprocprof/Vz', data=self.Vz[:,1])
        self.write_with_bound_info(h5f, '/preprocprof/Vth', data=self.Vth[:,1])

        self.write_with_bound_info(h5f, '/preprocprof/Er', data=self.Er[:,1])
        self.write_with_bound_info(h5f, '/preprocprof/q', data=self.q[:,1])

        self.write_fac_with_bound_info(h5f, '/factors/fac_n', data=[facs['fac_n']])
        self.write_fac_with_bound_info(h5f, '/factors/fac_Te', data=[facs['fac_Te']])
        self.write_fac_with_bound_info(h5f, '/factors/fac_Ti', data=[facs['fac_Ti']])
        self.write_fac_with_bound_info(h5f, '/factors/fac_vz', data=[facs['fac_vz']])
        
        h5f.close()
    
    def write_with_bound_info(self, file, dataset, data):
        ds = file.create_dataset(dataset, data=data)
        ds.attrs['lbounds'] = [1]
        ds.attrs['ubounds'] = [len(data)]

    def write_fac_with_bound_info(self, file, dataset, data):
        ds = file.create_dataset(dataset, data=data)
        ds.attrs['lbounds'] = [1]
        ds.attrs['ubounds'] = [len(data[0])]
        
