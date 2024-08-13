import sys
import os
sys.path.append(os.path.dirname(__file__) + '/../kilca_interface/')
from read_in import read_in
from change_opts import change_opts_fieldpy
from save_file import save_file
import subprocess
import os
import shutil

class fieldpy:

    field_divB0_inp = {}
    path_to_fourier_modes_exe = os.path.dirname(__file__) + '/../../matlab/fourier/'
    
    def __init__(self, gfile, pfile, convex_file, fluxdata_path):
        '''
        Constructor.
        input:
            gfile ... path to g-file (EQDSK standard)
            pfile ... path to p-file (field.dat)
            convex_file ... path to convexfile
            fluxdata_path ... path to fluxdata directory
        '''
        
        # default values of field_divB0.inp:
        self.field_divB0_inp['ipert']         = 0    #0=eq only, 1=vac, 2,3=vac+plas
        self.field_divB0_inp['iequil']        = 1    # 0=pert. alone, 1=with equil.
        self.field_divB0_inp['ampl']          = 1.00 # amplitude of perturbation, a.u.
        self.field_divB0_inp['ntor']          = 72   # number of toroidal harmonics
        self.field_divB0_inp['cutoff']        = 0.99 # inner cutoff in psi/psi_a units
        self.field_divB0_inp['icftype']       = 4    # type of coil file
        self.field_divB0_inp['gfile']         = ''   # equilibrium file
        self.field_divB0_inp['pfile']         = ''   # field.dat - coil file (field.dat) from Kisslinger_asdex or coil_tools
        self.field_divB0_inp['convexfile']    = ''   # convex file for stretchcoords
        self.field_divB0_inp['fluxdatapath']  = ''   # directory with data in flux coord.
        self.field_divB0_inp['winsize_Rfilt'] = 0    # window size for filtering of psi array over R
        self.field_divB0_inp['winsize_Zfilt'] = 0    # window size for filtering of psi array over Z

        self.gfile = gfile
        self.pfile = pfile
        self.convex_file = convex_file
        self.fluxdata_path = fluxdata_path
    
    def write_field_divB0_inp(self, path_from, path_to):

        self.existing = read_in(path_from)
        
        self.field_divB0_inp['gfile'] = "'" + self.gfile + "'" 
        self.field_divB0_inp['pfile'] = "'" + self.pfile + "'"
        self.field_divB0_inp['convexfile'] = "'" + self.convex_file + "'"
        self.field_divB0_inp['fluxdatapath'] = "'" + self.fluxdata_path + "'"

        self.ind = [6, 7, 8, 9]

        self.changed = change_opts_fieldpy(self.existing, self.ind, [self.field_divB0_inp['gfile'], self.field_divB0_inp['pfile'], self.field_divB0_inp['convexfile'], self.field_divB0_inp['fluxdatapath']], ' ')

        save_file(self.changed, path_to)

    def run_fourier_modes(self):
        
        try:
            wd = os.getcwd()
            os.chdir(self.path_to_fourier_modes_exe)
            process = subprocess.Popen('./fouriermodes.x', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

            for line in iter(process.stdout.readline, ''):
                print(line.strip())


            files = ['btor_rbig.dat', 'equil_r_q_psi.dat', 'axis.dat', 'box_size.dat', 'separ.dat', 'phinorm_arr.dat', 'thetabooz.dat', 'theta_of_theta_qt_flabel.dat', 'amn.dat']
            
            for el in files:
                if os.path.exists(self.fluxdata_path + el):
                    os.remove(self.fluxdata_path + el)
                shutil.move(el, self.fluxdata_path, copy_function=shutil.copy2)

            os.chdir(wd)
        except Exception as e:
            print(f'An error occurred: {e}')
