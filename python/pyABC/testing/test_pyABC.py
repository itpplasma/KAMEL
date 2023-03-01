
import sys
sys.path.append('../')
from pyABC import *

fluxpath = './fluxdata/'
gfilepath = './g33353.2900_EQH'
coilpath = '/proj/plasma/DATA/BALANCE/COIL/33353/33353.2900_coil.dat'

abc = pyABC(33353, 2900, './run/')
abc.set_coil(coilpath)
abc.set_equi(gfilepath, fluxpath)
