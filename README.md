#Balance Suite
This directory contains everything concerning the balance code in the hdf5 version. When using the template scripts, make sure to change the code according to your project before using it.
Most of the time, paths need to be changed. 

 - postproc_py_class
	Contains a python class used for the post processing of the hdf5 version of the balance code
	Note however, that this class is currently a mess and needs some changes.

 - template_scripts
	Contains template scripts.
		- script_prerun.m is the Matlab script that does the prerun. It creates an "input" hdf5 file for the balance code. This input file is the same for all types of runs (time evolution, parameter scan, ...). So, it has to be created only once. script_prerun_multiple.m runs the prerun for multiple times for a single shot.
		- /linearrun/ contains scripts to run the balance code in a linear fashion. This is used to determine the quasilinear diffusion coefficients for the full RMP coil current. There are templates for a single time/shot combination and multiple times (single shot).
		- /timeevol/ contains scripts to do time evolution runs of the balance code.
		- /parameterscan/ contains scripts to do parameter scans with the balance code.

 - ql-balance
	Contains the balance Fortran code itself. Compile with "make -f Balance.mk_mpi".

 - utility_scripts
	Contains matlab and python scripts that make certain tasks easier.

 - matlab
	Contains the matlab interface classes for the ql-balance, KiLCA and GPEC code, as well as things like NEO-2 and the kisslinger code. Also, blueprints for e.g. balance_conf.nml can be found there.

- Documentation
	Self explanatory.
