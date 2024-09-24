import numpy as np
from scipy.interpolate import CubicSpline

class Profile_Extender:
    
    d = 0.1
    dr_cut = 0.2
    
    cut_exp = lambda self, x, xc, d: 1.0/(1.0 + np.exp((x-xc)/d))
    cut_ee = lambda self, x, xc, d: np.exp(-np.exp((x-xc)/d))

    def __init__(self, file_in, output_file, factor):
        self.file_in = file_in
        self.output_file = output_file
        self.factor = factor

    def read(self):
        dat = np.loadtxt(self.file_in)
        self.psi_in = dat[:, 0]
        self.y_in = dat[:, 1]

    def process(self, r_eff, r_min, r_max, y_inf, type_of_cut):
        
        r_end = r_eff[-1]
        dr = r_eff[-1] - r_eff[-2]
        self.r_out = np.concatenate((r_eff, np.arange(r_end + dr, r_max, dr)))
        
        self.r_eff = np.append(r_eff, self.r_out[-1])
        self.y_in = np.append(self.y_in, y_inf)
        self.y_interp = CubicSpline(self.r_eff, self.y_in * self.factor, bc_type='natural')
        self.y_inf = y_inf * self.factor
        
        cut = r_end + self.dr_cut
        
        if type_of_cut == "" or type_of_cut == "exp":
            cut_fun = self.cut_exp
        elif type_of_cut == "ee":
            cut_fun = self.cut_ee
        
        self.y_out = self.y_interp(self.r_out) * cut_fun(self.r_out, cut, self.d) + self.y_inf * (1.0 - cut_fun(self.r_out, cut, self.d))

    def write(self):
        np.savetxt(self.output_file, np.column_stack((self.r_out, self.y_out)))