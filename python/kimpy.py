import subprocess
import os
import re
import shutil
import numpy as np

class kimpy:

    kim_config_nml = os.path.dirname(__file__) + '/../nmls/KIM_config.nml'
    kim_exe_path = os.path.dirname(__file__) + '/../KIM_exe'

    def __init__(self, runpath):

        self.runpath = runpath

        self.command = './KIM_exe'
        if os.path.exists(self.runpath + 'KIM_exe'):
            os.remove(self.runpath + 'KIM_exe')
        shutil.copy2(self.kim_exe_path, self.runpath + 'KIM_exe')
        
        shutil.copy2(self.kim_config_nml, self.runpath + 'KIM_config.nml')
        self.kim_config_nml = self.runpath + 'KIM_config.nml'
        

    def run(self, no_out=False):
        try:
            cwd = os.getcwd()
            os.chdir(self.runpath)
            print(os.getcwd())
            # Run the command and capture the output in real-time
            process = subprocess.Popen(self.command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

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
        r = np.linspace(0, r_max, 1000)
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
       os.symlink(from_path, self.runpath + 'profiles/') 