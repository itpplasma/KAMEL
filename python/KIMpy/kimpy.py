import subprocess
import os
from pathlib import Path
import re
import shutil
import numpy as np
import h5py

def get_base_path() -> Path:
    '''
        Check if the current KAMEL path is the one of the environment variable CODE.
        If not, assume this is a development repo and use this path instead.
    '''
    env_base = Path(os.environ['CODE']).resolve()
    code_base = Path(__file__).resolve().parents[2]

    # If code is not inside the env base, assume dev checkout
    if not code_base.is_relative_to(env_base):
        return code_base

    return env_base

CODE = get_base_path()

class KIMpy:

    kim_config_nml = str(CODE.absolute()) + '/KAMEL/KIM/nmls/KIM_config.nml'
    kim_exe_path = str(CODE.absolute()) + '/KAMEL/build/install/bin/KIM.x'
    omp_num_threads = 8

    def __init__(self, runpath):

        self.runpath = runpath

        self.command = './KIM.x'

        if not os.path.isfile(self.runpath + 'KIM.x'):
            os.system('ln -sf ' + self.kim_exe_path + ' '+ self.runpath + 'KIM.x')
        
        if not os.path.isfile(self.runpath + 'KIM_config.nml'):
            shutil.copy2(self.kim_config_nml, self.runpath + 'KIM_config.nml')
        self.kim_config_nml = self.runpath + 'KIM_config.nml'
        

    def run(self, no_out=False):
        try:
            cwd = os.getcwd()
            os.chdir(self.runpath)
            print(os.getcwd())
            # Run the command and capture the output in real-time
            process = subprocess.Popen(self.command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True, env=dict(os.environ, OMP_NUM_THREADS=str(self.omp_num_threads)))

            # Read and print the output line by line
            if (not no_out):
                for line in iter(process.stdout.readline, ''):
                    if not re.search(r'^Progress*', line):
                        print(line.strip())
                    else:
                        print(line.strip(), end='\r')

            # Wait for the process to complete
            process.wait()
            os.chdir(cwd)

        except Exception as e:
            print(f"Error: {str(e)}")


    def generate_constant_profiles(self, r_max=65.0, ne=1e13, Te=1e3, Ti=1e3, Vz=1e5, Er=0.5):
        self.profile_path = self.runpath + 'profiles/'
        r = np.linspace(0.1, r_max, 1000)
        ne = ne * np.ones_like(r)
        Te = Te * np.ones_like(r)
        Ti = Ti * np.ones_like(r)
        Vz = Vz * np.ones_like(r)
        Er = Er * np.ones_like(r)

        q = (r/np.max(r) * 2)**2 + 1.05
        
        np.savetxt(self.profile_path + 'n.dat', np.column_stack((r, ne)))
        np.savetxt(self.profile_path + 'Te.dat', np.column_stack((r, Te)))
        np.savetxt(self.profile_path + 'Ti.dat', np.column_stack((r, Ti)))
        np.savetxt(self.profile_path + 'Vz.dat', np.column_stack((r, Vz)))
        np.savetxt(self.profile_path + 'Er.dat', np.column_stack((r, Er)))
        np.savetxt(self.profile_path + 'q.dat', np.column_stack((r, q)))

    def take_existing_profiles(self, from_path):
        
        if os.path.islink(self.runpath + 'profiles/'):
            # Remove the symlink
            os.unlink(self.runpath + 'profiles/')
            print("Symlink removed successfully.")
        else:
            print("The specified path is either not a symlink or does not exist.")
        try:
            shutil.rmtree(self.runpath + 'profiles/')
        except:
            print("No profiles directory found in runpath.")
        os.symlink(from_path, self.runpath + 'profiles/') 

    def read_Eperp(self, m_mode=6, n_mode=2):
        try:
            E_perp = np.loadtxt(self.runpath + f'out/m{m_mode}_n{m_mode}/fields/E_perp.dat')
            E_perp_psi = np.loadtxt(self.runpath + f'out/m{m_mode}_n{m_mode}/fields/E_perp_psi.dat')
            E_perp_MA = np.loadtxt(self.runpath + f'out/m{m_mode}_n{m_mode}/fields/E_perp_MA.dat')
            return E_perp, E_perp_psi, E_perp_MA
        except FileNotFoundError:
            print("Eperp.dat not found in the runpath.")
            return None
    
    def transfer_solution_to_h5_reduced(self, h5file, m_mode=6, n_mode=2):

        h5f = h5py.File(h5file, 'w')
        phi_sol = np.loadtxt(self.runpath + f'out/m{m_mode}_n{n_mode}/fields/phi_sol.dat')
        E_perp = np.loadtxt(self.runpath + f'out/m{m_mode}_n{n_mode}/fields/E_perp.dat')
        E_perp_psi = np.loadtxt(self.runpath + f'out/m{m_mode}_n{n_mode}/fields/E_perp_psi.dat')
        E_perp_MA = np.loadtxt(self.runpath + f'out/m{m_mode}_n{n_mode}/fields/E_perp_MA.dat')

        h5f.create_dataset('r', data=phi_sol[:,0])
        h5f.create_dataset('phi_sol', data=phi_sol[:,1:3])
        h5f.create_dataset('E_perp', data=E_perp[:,1:3])
        h5f.create_dataset('E_perp_psi', data=E_perp_psi[:,1:3])
        h5f.create_dataset('E_perp_MA', data=E_perp_MA[:,1:3])
        h5f.close()
