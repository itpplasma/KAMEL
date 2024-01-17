import subprocess
import os
import re

class kimpy:
    def __init__(self, runpath):

        self.runpath = runpath
        self.command = './KIM_exe'
        

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
