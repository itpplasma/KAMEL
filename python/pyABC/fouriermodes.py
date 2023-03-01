import os
import subprocess
import datetime


def fouriermodes(fourierpath, fluxdatapath):
    ###########################################################################
    #function Fouriermodes(fourierpath, fluxdatapath, fdb0)
    ###########################################################################
    # description:
    #--------------------------------------------------------------------------
    # runs fouriermodes.x located in fourierpath for already written
    # field_divB0.inp file(MATLAB object fdb0) and saves output to
    # fluxdatapath.
    ###########################################################################
    # input:
    #--------------------------------------------------------------------------
    # fourierpath  ... location of fouriermodes.x
    # fluxdatapath ... location to save data
    ###########################################################################

    curr_path = os.getcwd()

    os.chdir(fourierpath)
    time_start = datetime.datetime.now()
    print('Start of Fouriermodes at ' + str(time_start))
    # execute fouriermodes
    res = subprocess.Popen('./fouriermodes.x', stdout=subprocess.PIPE)
    time_end = datetime.datetime.now()
    logname = fourierpath + 'fouriermodes_' + str(time_end) + '.log'
    fid = open(logname, 'w')
    fid.writelines('Start of Fouriermodes at ' + str(time_start) + '\n')
    fid.writelines('End of Fouriermodes at ' + str(time_end) + '\n\n')
    #fid.writelines(res)
    with res as proc:
        fid.write(str(proc.stdout.read()))
    fid.close()

    os.system('mkdir -p ' + fluxdatapath)
    outfiles = ['amn.dat', 'btor_rbig.dat', 'equil_r_q_psi.dat', 'axis.dat', 'box_size.dat', 'separ.dat', 'phinorm_arr.dat', 'thetabooz.dat', 'theta_of_theta_qt_flabel.dat', logname]

    for file in outfiles:
        os.system('mv -f ' + file + ' ' + fluxdatapath)

    os.system('cp -f field_divB0.inp ' + fluxdatapath)
    os.system('cp -f fouriermodes.inp ' + fluxdatapath)

    print('Finished Fouriermodes at ' + str(time_end))
    os.chdir(curr_path)
