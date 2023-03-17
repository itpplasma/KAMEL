import os
import f90nml
import subprocess
from datetime import datetime
from AUG_coil import *
import logging

class GPEC_interface:

    gpec_home = '/proj/plasma/CODE/GPEC/' # path to GPEC source files
    path_to_run = ''
    shot = 0
    time_slice = 0
    name = ''

    lib_GPEC_interface = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
    template_input = os.path.abspath(os.path.join(lib_GPEC_interface, '../matlab/balance/GPEC_interface/template/')) + '/'
    

    coil_nml  = []   # input file for coil
    equil_nml = []  # input file for equi
    dcon_nml  = []   # input file for dcon
    gpec_nml  = []   # input file for gpec
    vac_nml   = []    # input file for vac

    def __init__(self, runpath, shot, time, name):
        '''Constructor of the GPEC python interface.
        Input:
        - runpath ... path of run of the balance code
        - shot    ... shot number
        - time    ... time in ms
        - name    ... name of the run'''

        self.path_run = runpath
        self.shot = shot
        self.time = time
        self.name = name

    def setEqui(self, gfile):
        '''Sets equilibrium data for GPEC.
        Input:
        - gfile ... path of equilibrium file.'''

        if (not os.path.exists(gfile)):
            raise ValueError('Equi file not found for gpec.')
        self.equi_nml = f90nml.read(self.template_input + 'equil.in')
        self.equi_nml['EQUIL_CONTROL']['eq_filename'] = gfile

    def setCoil(self, coil_file):
        '''Sets coil data for GPEC
        Input:
        - coil_file ... path to coil file'''

        if (not os.path.exists(coil_file)):
            raise ValueError('Coil file not found for gpec.')

        # load coil file using AUG_coil class
        rawcoil = AUG_coil(coil_file)
        rawcoil.read()

        # rewrite coil data to list
        self.current = np.reshape([rawcoil.Iu, rawcoil.Il],(len(rawcoil.Iu),2)).tolist()

        self.coil_nml = f90nml.read(self.template_input + 'coil.in')
        self.coil_nml['COIL_CONTROL']['coil_cur'] = self.current

    def loadDGV(self):
        '''Description:
        loads dcon, gpec and vac namelists from template.'''

        self.dcon_nml = f90nml.read(self.template_input + 'dcon.in')
        self.gpec_nml = f90nml.read(self.template_input + 'gpec.in')
        self.vac_nml = f90nml.read(self.template_input + 'vac.in')


    def write(self):
        '''Description:
        Creates directory structure and copies all input files.'''

        # load namelists if not done yet
        if not self.dcon_nml or not self.gpec_nml or not self.vac_nml:
            self.loadDGV()

        # create run path
        os.system('mkdir -p ' + self.path_run)

        # delete old outputs
        os.system('rm ' + self.path_run + '*.nc')
        os.system('rm ' + self.path_run + '*.bin')
        os.system('rm ' + self.path_run + '*.out')
        os.system('rm ' + self.path_run + '*.log')

        # write namelist files
        if self.coil_nml:
            self.coil_nml.write(self.path_run + 'coil.in', force=True)
        if self.equi_nml:
            self.equi_nml.write(self.path_run + 'equil.in', force=True)
        
        self.dcon_nml.write(self.path_run + 'dcon.in', force=True)
        self.gpec_nml.write(self.path_run + 'gpec.in', force=True)
        self.vac_nml.write(self.path_run + 'vac.in', force=True)

    def run(self):
        '''Runs dcon and gpec.'''

        self.run_code('dcon')
        self.run_code('gpec')

    def run_code(self, code):
        """Either dcon or gpec"""
        pypath = os.getcwd()
        os.chdir(self.path_run)

        # create sotlink to executable
        sourcepath = self.gpec_home + code + '/' + code
        os.system('ln -sfT ' + sourcepath + ' ./' + code)

        # get time and create log file
        start_time = datetime.now()
        logfile = code + '_' + str(start_time).replace(' ', '_') + '.log'
        fid = open(logfile, 'w')
        fid.writelines('Start of ' + code + ' at ' + str(start_time) + '\n')
        fid.writelines('Shot: ' + str(self.shot) + ', Time: ' + str(self.time) + 'ms, Name: ' + self.name + '\n\n')

        print('Start of ' + code + ' at ' + str(start_time))
        print('Shot: ' + str(self.shot) + ', Time: ' + str(self.time) + 'ms, Name: ' + self.name)

        # execute code
        self.out_code = subprocess.Popen("".join(self.getConf()) + ' ; ./'+ code, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
        self.out_code.wait()

        #logging.info(self.out_code.communicate())
        with self.out_code as proc:
            fid.write(proc.stdout.read().decode('UTF-8'))
        fid.close()

        print('Finished ' + code + ' at ' + str(datetime.now()))
        print('Total runtime was ' + str(datetime.now() - start_time))

        os.chdir(pypath)

    def getConf(self):
        '''return config lines'''

        conf = ['source /afs/itp.tugraz.at/opt/intel/2018.1/bin/compilervars.sh intel64 ; ',
                    'export FC=ifort ; ',
                    'export F77=ifort ; ',
                    'export F90=ifort ; ',
                    'export CC=icc ; ',
                    'export CXX=icc ; ',
                    'export PROJLIBS=/proj/plasma/Libs/intel ; ',
                    'export NETCDFHOME=$PROJLIBS/NetCDF ; ',
                    'export CFLAGS=-I$NETCDFHOME/include/ ; ',
                    'export FFLAGS=-I$NETCDFHOME/include/ ; ',
                    'export CXXFLAGS=-I$NETCDFHOME/include/ ; ',
                    'export LDFLAGS=-L$NETCDFHOME/lib/ ; ',
                    'export LD_LIBRARY_PATH=$NETCDFHOME/lib/:$LD_LIBRARY_PATH ; ',
                    'export F90HOME=/afs/itp.tugraz.at/opt/intel/2018.1 ; ',
                    'export GPECHOME=', self.gpec_home, ' ; ',
                    'ulimit -s unlimited']

        return conf