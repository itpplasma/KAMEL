# Balance Code Framework
This repository contains the KiLCA/QL-Balance framework. 
When using the template scripts, make sure to change the code, most importantly the paths, according to your project before using it.

## Content
### postproc_py_class
Contains a python class used for the post processing of the hdf5 version of the balance code
Note however, that this class is currently a mess and needs some changes.

### template_scripts
Contains matlab scripts that can be used as templates for certain balance code runs.
	
### ql-balance
Contains the balance Fortran code itself. Compile with 
```
make -f Balance.mk_mpi
```
### KiLCA
Contains the source code of KiLCA. In the directory, compile with

```
mkdir build
cd build
cmake ..
make
```

Note that certain libraries (lapack-3.2.1, sundials-5.7.0 and gsl-2.4) have to be available. In the code framework of the ITP plasma group, use the setup_kilca.sh shell script in the scripts directory.
So far, the compilation and execution of the (Normal, Release, NOMD, FPGEN) version of the code was tested on Linux and MacOS. 

### utility_scripts
Contains matlab and python scripts that make life easier.

### matlab
Contains the matlab interface classes for the ql-balance, KiLCA and GPEC code, as well as things like NEO-2 and the kisslinger code. Also, blueprints for e.g. balance_conf.nml can be found there.

### Documentation
- Short introduction to the balance code framework.
- List of variables contained in the balance configuration namelist balance_conf.nml.
