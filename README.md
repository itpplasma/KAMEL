#Balance Suite
This directory contains everything concerning the balance code in the hdf5 version. When using the template scripts, make sure to change the code according to your project, before using it.
Most importantly, the paths need to be changed. Also, the ion mass can be changed.

 - postproc_py_class
	contains a python class used for the post processing of the hdf5 version of the balance code

 - template_scripts
	contains the scripts. script_prerun.m is the Matlab script that must be run first. It creates an "input" hdf5 file for the balance code. This input file is the same for all types of runs (time evolution, parameter scan, ...). So, it has to be created only once.

 - ql-balance
	contains the balance code itself. The matlab class is contained in libneo.

 - utility_scripts
	contains matlab and python scripts that make certain tasks easier.

 - matlab
	contains the matlab interface classes for the ql-balance, KiLCA and GPEC code, as well as things like NEO-2 and the kisslinger code.
