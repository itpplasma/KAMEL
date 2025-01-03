# KAMEL - Kinetic plAsma response ModEL
This repository contains the kinetic plasma response framework containing the linear plasma response codes KiLCA and KIM, as well as the quasilinear transport code QL-Balance. KiLCA and KIM are cylindrical linear plasma response solvers based on a finite Larmor radius and an integral formalism, respectively. QL-Balance is a quasilinear 1D radial transport code. In combination, they are used to model the plasma response to external magnetic perturbation in toroidally confined fusion plasmas. 

Note, when using the template scripts, make sure to change the code, most importantly the paths, according to your project before using it.

## Compilation
For the **initial** compilation of the whole framework, use the top-level shell script

```
./kamel.sh
```

This invokes the sub-level shell scripts of the individual codes, which can also be called individually (see below). Note that in the compilation process the folder 'external' is created which contains libraries required for the compilation. It is assumed that the libraries are not present on the system and the libraries are therefore downloaded during the execution of the script. For more details, see the respective scripts.

Generally, for Apple Silicon the clang/gfortran (version 16.0 and 14.2., respectively) compiler combination is tested. On debian, the gnu compiler version 12.2.0 is tested.

After code changes, use the usual cmake and make commands in the respective colde folders for building.

## Codes

### KiLCA
Contains the source code of KiLCA. Compile in top level directory KAMEL with the shell script

```
./KiLCA/kilca.sh
```

This also installs the necessary external libraries like gsl, sundials, etc in the 'external' directory.
So far, the compilation and execution of the (Normal, Release, NOMD, FPGEN) version of the code was tested on Linux and MacOS. 

### KIM
Contains the source code of KiLCA and additional python code (e.g. calculation of the dispersion relation). Compile in top level directory KAMEL with the shell script

```
./KIM/kim.sh
```

### QL-Balance
Quasilinear transport code based on KiLCA. Requires the prior compilation of KiLCA. Compile in QL-Balance folder with

```
make
```

### PreProc
PreProc contains the fouriermodes code used to calculate r_eff, q, and the toroidal and poloidal fluxes. Also, it contains the neo-2 templates used to run NEO-2 on the ITP machines with condor. This requires the NEO-2 code (see github.com/itpplasma/neo-2).

## python
Contains python classes and functions to use the code.

## template_scripts
Contains matlab scripts that can be used as templates for certain balance code runs.

## utility_scripts
Contains matlab and python scripts that make life easier.

## matlab
Contains the matlab interface classes for the ql-balance, KiLCA and GPEC code, as well as things like NEO-2 and the kisslinger code. Also, blueprints for e.g. balance_conf.nml can be found there.

## Documentation
- Short introduction to the balance code framework.
- List of variables contained in the balance configuration namelist balance_conf.nml.
