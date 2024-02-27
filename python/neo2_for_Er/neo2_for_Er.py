import os
import numpy as np
import shutil
import subprocess
import h5py
import re
import time
from datetime import datetime
import sys
sys.path.append(os.path.dirname(__file__) + '/../fieldpy/')
from fieldpy import fieldpy

class neo2_for_Er():

    pmass = 1.6726e-24
    echarge = 4.8031e-10
    eV = 1.6022e-12

    neo2_path = os.path.normpath(os.path.dirname(__file__)+'/../../matlab/neo2/')+'/'

    def __init__(self, profilepath, equil_r_q_psi):


        def remove_non_monotonous_tail(data):
            non_mon_index = next((i for i, (x,y,) in enumerate(zip(data[:-1], data[1:])) if x > y), None)
            print(non_mon_index)
            if non_mon_index is not None:
                return [data[:non_mon_index+1], non_mon_index]
            else:
                return [data, len(data)-1]

        self.kpath = profilepath + 'kprof/'
        self.kprof = self.kpath + 'k.dat'

        self.equil_r_q_psi = np.loadtxt(equil_r_q_psi, skiprows=3)
        self.r_eff = self.equil_r_q_psi[:,0]
        [self.r_eff, ind] = remove_non_monotonous_tail(self.r_eff)

        self.S = self.equil_r_q_psi[:ind+1,2] / self.equil_r_q_psi[ind,2]
        self.R_beg = self.equil_r_q_psi[:ind+1,7]
        self.Z_beg = self.equil_r_q_psi[:ind+1,8]

        self.ne = np.loadtxt(profilepath + 'n.dat')[:,1]
        self.Te = np.loadtxt(profilepath + 'Te.dat')[:,1]
        self.Ti = np.loadtxt(profilepath + 'Ti.dat')[:,1]
        self.Vz = np.loadtxt(profilepath + 'Vz.dat')[:,1]



    def run_neo2(self, gfile, convex_wall, flux_data):
        '''Run NEO-2 to get the Er profile.'''

        os.makedirs(self.kpath, exist_ok=True)
        os.system('rm -r ' + self.kpath + '*')

        # prepare profiles for neo-2
        # set points for neo-2

        #self.r_eff = self.equil_r_q_psi[:,0]
        self.r_min = self.r_eff[0] + 1.0

        n_points = [5,3,9] # points between borders

        self.r_border = [[self.r_min, 10], [11, 53], [58, np.max(self.r_eff)-1]]

        r_neo = np.zeros(np.sum(n_points))
        k2 = 0

        for l in np.arange(0,len(self.r_border)):
            k1 = k2
            k2 = k1 + n_points[l]
            r_neo[k1:k2] = np.linspace(self.r_border[l][0], self.r_border[l][1], n_points[l])

        s_neo = np.interp(r_neo, self.r_eff, self.S)
        r_beg = np.interp(r_neo, self.r_eff, self.R_beg)
        Z_beg = np.interp(r_neo, self.r_eff, self.Z_beg)
        Ti_neo = np.interp(r_neo, self.r_eff, self.Ti)
        ne_neo = np.interp(r_neo, self.r_eff, self.ne)        

        # calculate kappa from ions
        collog = 39.1 - 1.15 * np.log10(1e6 * ne_neo) + 2.3 * np.log10(1e-3 * Ti_neo)
        v_th = np.sqrt(2 * Ti_neo * self.eV / self.pmass)
        tau = 3 * self.pmass**2 * v_th**3 / (16 * np.sqrt(np.pi) * ne_neo * self.echarge**4 * collog)
        kappa = - 2 / (v_th * tau)

        # prepare run for neo-2

        # copy content
        if os.path.exists(self.kpath):
            shutil.rmtree(self.kpath)

        shutil.copytree(self.neo2_path, self.kpath)
        #for filename in os.listdir(self.neo2_path):
        #    source_file = os.path.join(self.neo2_path, filename)
        #    if os.path.isfile(source_file):
        #        destination_file = os.path.join(self.kpath, filename)
        #        shutil.copy2(source_file, destination_file)

        shutil.copy2(convex_wall, os.path.join(self.kpath, 'convexwall.dat'))

        # output matrix
        M = np.column_stack((s_neo, r_neo, r_beg, Z_beg, Ti_neo, ne_neo, kappa))
        np.savetxt(self.kpath + 'surfaces.dat', M)

        self.fp = fieldpy(gfile, '',  convex_wall, flux_data)
        self.fp.field_divB0_inp['ipert'] = 0
        self.fp.write_field_divB0_inp(self.fp.path_to_fourier_modes_exe + 'template_field_divB0.inp', self.kpath + 'TEMPLATE_DIR/field_divB0.inp')

        wd = os.getcwd()
        os.chdir(self.kpath)

        subprocess.run(['python', 'create_surf_realspace.py'])

        self.run_remote_neo2(self.neo2_path + 'remote_run.conf', self.kpath)


    def run_remote_neo2(self, configfile, kpath):
        """Run NEO-2 on computers specified in the config file.
            input:
                configfile ... path to config file
                kpath .. path of run 
        """

        minram = 30000.0 # minimum RAM per NEO-2 run
        mincpu = 1.0 # minimum number of CPUs per NEO-2 run

        host = []
        cpu = []
        mem = []
        maxjob = []

        # get the computers, cpu and ram from the config file
        with open(configfile, 'r') as file:
            for line in file:
                if line[0] == '#':
                    continue
                # split the line
                values = line.split()
                host.append(values[0])
                cpu.append(int(values[1]))
                mem.append(int(values[2]))
                maxjob.append(int(np.min([np.floor(int(values[2]) / minram), np.floor(int(values[1]) / mincpu)])))


        queue = []
        for k in np.arange(0, len(maxjob)):
            queue.append(np.tile(host[k], maxjob[k]).tolist())

        wd = os.getcwd()
        os.chdir(self.kpath)

        jobs = [d for d in os.listdir(self.kpath) if os.path.isdir(os.path.join(self.kpath, d))]
        jobs = [d for d in jobs if d != 'TEMPLATE_DIR']
        totaljobs = len(jobs)

        jobsrun = []
        queuerun = []
        jobsfinished = []
        queuewait = queue.copy()
        jobswait = jobs.copy()

        start_time = time.time()
        print('Start NEO-2 at ', datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        finished = []
        for k in np.arange(0,len(jobswait)):
            if os.path.exists(jobswait[k] + '/done.out'):
                print('here in done')
                finished.append(k)
                jobsfinished.append(jobswait[k])
                jobswait.remove(k)

        while True:
            if len(jobsrun) > 0:
                finished = []
                jobsrun_count = len(jobsrun)
                jobsrun_copy = jobsrun.copy()
                queuerun_copy = queuerun.copy()

                for k, job in enumerate(jobsrun):
                    if os.path.exists(job + '/log.txt'):
                        with open(job + '/log.txt') as file:
                            for line in file:
                                if re.compile(r'I give up').search(line):
                                    os.system(f'touch {job}/done.out')
                                    print(f'Got I give up in: {job}')
                    if os.path.exists(job + '/done.out'):
                        print(f'Remove due to done.out in {job}')
                        finished.append(k)
                        jobsfinished.append(job)
                        queuewait.append(queuerun_copy[k])
                        jobsrun.remove(job)
                        queuerun.remove(queuerun_copy[k])
            if len(queuewait) > 0 and len(jobswait) > 0:
                start = []
                for k in np.arange(np.min([len(queuewait), len(jobswait)])):
                    #command = f"ssh {queuewait[k][0]} ''cd {self.kpath}; cd {jobswait[k]}/ ; hostname > log.txt ; ./neo_2.x >> log.txt 2>&1 ; touch done.out &''"
                    command =f"ssh {queuewait[k][0]} ''cd {self.kpath}; cd {jobswait[k]}/; hostname > log.txt; ./neo_2.x >> log.txt 2>&1; touch done.out'' &"
                    jobnum = len(jobsfinished) + len(jobsrun) + 1
                    print(f'Start job {jobnum}/{totaljobs} : {jobswait[k]} @ {queuewait[k][0]}')
                    #os.system(command)
                    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    #print(f"Started at {queuewait[k][0]}")
                    jobsrun.append(jobswait[k])
                    queuerun.append(queuewait[k])
                    start.append(k)
                    stdout, stderr = proc.communicate()
                    if stderr is not None:
                        print("SSH stderr:", stderr)
                queuewait = [value for index, value in enumerate(queuewait) if index not in start]
                jobswait = [value for index, value in enumerate(jobswait) if index not in start]

            print(f'Time elapsed: {time.time() - start_time:.2f} s; Jobs: {len(jobsfinished)} / {totaljobs} finished')

            if len(jobsfinished) == totaljobs:
                break

            time.sleep(10)
        print(f'Finished NEO-2 at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
        os.chdir(wd)
