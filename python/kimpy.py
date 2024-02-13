import subprocess
import os
import re
import shutil

class kimpy:

    kim_config_nml = os.path.dirname(__file__) + '/../nmls/KIM_config.nml'
    kim_exe_path = os.path.dirname(__file__) + '/../KIM_exe'

    def __init__(self, runpath):

        self.runpath = runpath

        self.command = './KIM_exe'
        shutil.copyfile(self.kim_exe_path, self.runpath + 'KIM_exe')

        shutil.copyfile(self.kim_config_nml, self.runpath + 'KIM_config.nml')
        self.kim_config_nml = self.runpath + 'KIM_config.nml'
        

    def run(self):
        try:
            cwd = os.getcwd()
            os.chdir(self.runpath)
            print(os.getcwd())
            # Run the command and capture the output in real-time
            process = subprocess.Popen(self.command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

            # Read and print the output line by line
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
