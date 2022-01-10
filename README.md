#Balance Suite
This repository contains everything concerning the balance code in the hdf5 version. To do runs, a compiled version of KiLCA needs to be available.
When using the template scripts, make sure to change the code, most importantly the paths, according to your project before using it.

## Content
### postproc_py_class
	Contains a python class used for the post processing of the hdf5 version of the balance code
	Note however, that this class is currently a mess and needs some changes.

### template_scripts
	Contains matlab scripts that can be used as templates for certain balance code runs.
	
### ql-balance
	Contains the balance Fortran code itself. Compile with 

		make -f Balance.mk_mpi

### utility_scripts
	Contains matlab and python scripts that make life easier.

### matlab
	Contains the matlab interface classes for the ql-balance, KiLCA and GPEC code, as well as things like NEO-2 and the kisslinger code. Also, blueprints for e.g. balance_conf.nml can be found there.

### Documentation
	- Short introduction to the Balance Suite.
	- List of variables contained in the balance configuration namelist balance_conf.nml.
