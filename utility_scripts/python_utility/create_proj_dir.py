#!/usr/bin/env python3
###############################################################################
# create_proj_dir.py
# ----------------------------------------------------------------------------
# description:
# ------------
# This python script creates the directory structure for a project running 
# the balance code. Also, it copies the template scripts to the SCRIPTS dir.
# ----------------------------------------------------------------------------
# author: Markus Markl
# created: 22.12.21
###############################################################################

import os
import shutil

#parentdir = '/temp/markl_m/' + projectname
projectname = 'test_proj'
parentdir = './'+projectname

datadir = os.path.join(parentdir, 'DATA')
coildir = 'COIL'
equidir = 'EQUI'
indexdir = 'INDEX'
profdir = 'PROF'

os.makedirs(os.path.join(datadir, coildir))
os.makedirs(os.path.join(datadir, equidir))
os.makedirs(os.path.join(datadir, indexdir))
os.makedirs(os.path.join(datadir, profdir))


runsdir = os.path.join(parentdir, 'RUNS')
os.makedirs(runsdir)

runspredir = os.path.join(parentdir, 'RUNS_PRE')
os.makedirs(runspredir)

scriptsdir = os.path.join(parentdir, 'SCRIPTS')

# copy template scripts
source_folder = '../../template_scripts'
destination_folder = scriptsdir
shutil.copytree(source_folder, destination_folder)

postprocdir = 'postprocessing'
os.makedirs(os.path.join(scriptsdir, postprocdir))


prerundatadir = os.path.join(parentdir, 'PRERUNDATA')
os.makedirs(prerundatadir)


