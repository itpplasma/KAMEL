import h5py
import numpy as np

class Balance_Input_h5:
    """ Class to prepare the input h5 file for the balance code."""
    
    def __init__(self, h5_file_path, profile_path):
        self.profile_path = profile_path

    def get_required_data(self):
        self.Da = np.loadtxt(self.profile_path + '/Da.dat')
        self.n = np.loadtxt(self.profile_path + '/n.dat')
        self.Te = np.loadtxt(self.profile_path + '/Te.dat')
        self.Ti = np.loadtxt(self.profile_path + '/Ti.dat')
        self.Vz = np.loadtxt(self.profile_path + '/Vz.dat')
        self.Vth = np.loadtxt(self.profile_path + '/Vth.dat')
        self.Er = np.loadtxt(self.profile_path + '/Er.dat')
        self.q = np.loadtxt(self.profile_path + '/q.dat')
    
    def write_data_to_h5(self, file_name):
        print('Writing data to h5 file: ', file_name)
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
        h5f.close()
    
    def write_with_bound_info(self, file, dataset, data):
        ds = file.create_dataset(dataset, data=data)
        ds.attrs['lbounds'] = 0
        ds.attrs['ubounds'] = len(data)
        