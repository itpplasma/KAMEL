#! /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
from os.path import exists

sys.path.insert(0,__file__ + '/../../postproc_py_class')
#from utility_class import utility
#from balancepost import balance
#import utility

try:
    path = sys.argv[1]
except:
    path = './'

print("Path: " + path)

file_exists = exists(path + 'EB.dat')

if not file_exists:
    raise ValueError("EB.dat does not exist in " + path)


data = np.loadtxt(path+ 'EB.dat')

#ut = utility()

fig,ax = plt.subplots()
ax.plot(data[:,0], data[:,7], label=r'Re($B^r$)', lw=2)
ax.plot(data[:,0], data[:,8], label=r'Im($B^r$)', lw=2)
ax.legend()
ax.set_xlabel('r / cm')
ax.set_ylabel('B / G')
#ut.add_grid_to_axis(ax)
plt.tight_layout()
plt.show()