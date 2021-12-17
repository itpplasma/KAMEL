#!/usr/bin/env python3
#############################################################################
# fetch_aug_data_from_list.py
#############################################################################
# This script fetches data from the IPP network that was previously generated
# by the use of augped and the script_collect.py script. It reads the shot 
# and time combinations from a text file.
#
# NOTE: For usage, change the paths according to your needs.
#############################################################################
# author: Markus Markl
# created: 11.2021
#############################################################################

import paramiko
from scp import SCPClient
import sys
import os
import getpass
import numpy as np

    # progress function, to show the progress of the copying
def progress4(f, size, sent, p):
    prog = sent / size * 100.
    sys.stdout.write(f"({p[0]}:{p[1]}) {f}\'s progress: {prog}\r\n")
    sys.stdout.write(" - - - \n")

policy = paramiko.client.AutoAddPolicy

# ask for username and password
user = input("Username: ")
passwd = getpass.getpass("Password: ")

listpath = '/itp/MooseFS/markl_m/BALANCE_4CASES/DATA/shottimelist.txt'
data = np.loadtxt(listpath, dtype='int')

shots = data[:,0]
shots = np.unique(shots)
time_slices = data[:,1]
profile_offset = data[:,2]

prof_time = time_slices + profile_offset

#shot = input("shot number: ")
#time_slice = input("at time (in ms): ")
# also check for the time of the profiles, it may be different by a few ms
#profile_offset = input("profile offset: ")

#prof_time = str(int(time_slice)+int(profile_offset))
profile_string = "_PED_MMARKL_rho_pol.dat"

local_project_location = "/temp/markl_m/BALANCE_4CASES/DATA/"

for shot in shots:
    local_path_coil = local_project_location + "COIL/" + str(int(shot)) + "/"
    local_path_equi = local_project_location + "EQUI/" + str(int(shot)) + "/"
    local_path_prof = local_project_location + "PROF/" + str(int(shot)) + "/"

    # check if the paths exist

    if not os.path.exists(local_path_coil):
        print("coil dir for "+str(int(shot))+ " does not exist, will be created")
        try:
            os.mkdir(local_path_coil)
        except:
            raise OSError("Can't create directory " + local_path_coil)
    else:
        print('coil dir exists')

    if not os.path.exists(local_path_equi):
        print("equi dir does not exist, will be created")
        try:
            os.mkdir(local_path_equi)
        except:
            raise OSError("Can't create directory " + local_path_equi)
    else:
        print("equi dir exists")

    if not os.path.exists(local_path_prof):
        print("prof dir does not exist, will be created")
        try:
            os.mkdir(local_path_prof)
        except:
            raise OSError("Can't crate directory " + local_path_prof)
    else:
        print("prof dir exists")

# establish ssh connection
    with paramiko.SSHClient() as client:
        client.set_missing_host_key_policy(policy)
    try:
        client.connect(hostname='gate2.aug.ipp.mpg.de', username=user, password=passwd)
        print("Connection established")
    except paramiko.ssh_exception.NoValidConnectionsError:
            print("Connection failed")

    # paths in the ipp network

    count = 0
    for t in time_slices:
        remote_coil = "~/pyout/coil/" + str(int(shot)) + "/" + str(int(shot)) + "." + str(t) + "_coil.dat"
        remote_equi = "~/pyout/equi/" + str(int(shot)) + "/g" + str(int(shot)) + "." + str(t) + "_EQH"
        remote_prof = "~/pyout/profiles/" + str(int(shot)) + "/" + str(int(shot)) + "." + str(int(prof_time[count])) + "_"
        remote_prof_ne = remote_prof + "ne" + profile_string
        remote_prof_Te = remote_prof + "Te" + profile_string
        remote_prof_Ti = remote_prof + "Ti" + profile_string
        remote_prof_vt = remote_prof + "vt" + profile_string
            
        # get the data
        with SCPClient(client.get_transport(), progress4=progress4) as scp:
            scp.get(local_path=local_path_coil, remote_path=remote_coil)
            scp.get(local_path=local_path_equi, remote_path=remote_equi)
            scp.get(local_path=local_path_prof, remote_path=remote_prof_ne)
            scp.get(local_path=local_path_prof, remote_path=remote_prof_Te)
            scp.get(local_path=local_path_prof, remote_path=remote_prof_Ti)
            scp.get(local_path=local_path_prof, remote_path=remote_prof_vt)

        count = count+1



