#!/usr/bin/env python3
#####################################################################
# extended_4Cases_make_tex.py
#####################################################################
# This script creates a tex file, using the plot files contained in 
# a directory. 
#####################################################################

import os
import csv

indexpath = '/temp/markl_m/BALANCE_4CASES/DATA/INDEX/'
indexfiles = filter(lambda x: x.endswith('.index'), os.listdir(indexpath))
texpath = '/temp/markl_m/BALANCE_4CASES/slides_velocities/inserts/'

timeslices = []
types = []

for ind in indexfiles:
    print(ind)
    shot, _ = ind.split('.')
    texfile = 'plots'+str(shot)+'.tex'
    
    # read index file to get times and types
    with open(indexpath+ind, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            timeslices.append(row[0])
            types.append(row[1])
    #print(timeslices)
    #print(types)
    # sort times
    timeslices, types = (list(t) for t in zip(*sorted(zip(timeslices, types))))

    # create tex file
    count = 0
    currenttype =''

    with open(texpath+texfile, 'w') as f:
        pass

    for t in timeslices:
        if 'MMARKL' in types[count]:
            currenttype = 'MMARKL'
        elif 'ULBLP' in types[count]:
            currenttype = 'ULBLP'
        else:
            print('Unknown profile type')
        

        with open(texpath+texfile, 'a') as f:
            # velocity plot
            f.write('\\begin{frame}\n')
            f.write('\t\\frametitle{'+str(shot)+'/'+t+' '+currenttype+'}\n')
            f.write('\t\\centering\n')
            f.write('\t\\includegraphics[width=1.0\\textwidth]{plots/'+types[count]+'_NRef_'+str(shot)+'_'+t+'_vel.png}\n')
            f.write('\\end{frame}\n\n')
            # pressure plot
            f.write('\\begin{frame}\n')
            f.write('\t\\frametitle{'+str(shot)+'/'+t+' '+currenttype+', pressure profiles}\n')
            f.write('\t\\centering\n')
            f.write('\t\\includegraphics[width=1.0\\textwidth]{plots/'+str(shot)+'_'+t+'_pressure.png}\n')
            f.write('\\end{frame}\n\n')
        
        count = count + 1


    timeslices =[]
    types = []
    
    print(str(shot) + ' done')






