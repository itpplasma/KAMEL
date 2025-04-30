# KAMEL - Kinetic plAsma response ModEL
This repository contains the kinetic plasma response framework containing the linear plasma response codes KiLCA and KIM, as well as the quasilinear transport code QL-Balance. KiLCA and KIM are cylindrical linear plasma response solvers based on a finite Larmor radius and an integral formalism, respectively. QL-Balance is a quasilinear 1D radial transport code. In combination, they are used to model the plasma response to external magnetic perturbation in toroidally confined fusion plasmas. 

Note, when using the template scripts, make sure to change the code, most importantly the paths, according to your project before using it.

## Compilation
For compilation of the codes use in the repositories root directory either
```
make all
```
to compile all three codes (QL-Balance needs a compiled version of KiLCA in any case), or 
```
make $CODE_NAME
```
to compile each code individually. Each specific code directory comes also with its own Makefile.

Generally, for Apple Silicon the clang/gfortran (version 16.0 and 14.2., respectively) compiler combination is tested. On debian, the gnu compiler version 12.2.0 is tested.

## Codes

### KiLCA
Contains the source code of KiLCA.

So far, the compilation and execution of the (Normal, Release, NOMD, FPGEN) version of the code was tested on Linux and MacOS. 

### KIM
Contains the source code of KIM (KiLCA Integral Model).

### QL-Balance
Quasilinear transport code based on KiLCA. Requires the prior compilation of KiLCA.

### PreProc
PreProc contains the fouriermodes code used to calculate r_eff, q, and the toroidal and poloidal fluxes. Also, it contains the neo-2 templates used to run NEO-2 on the ITP machines with condor. This requires the NEO-2 code (see github.com/itpplasma/neo-2).

## python
Contains python classes and functions to use the code. Comes with its own Makefile.

## template_scripts
Contains matlab scripts that can be used as templates for certain balance code runs.

## utility_scripts
Contains matlab and python scripts that make life easier.

## matlab
Contains the matlab interface classes for the ql-balance, KiLCA and GPEC code, as well as things like NEO-2 and the kisslinger code. Also, blueprints for e.g. balance_conf.nml can be found there.

## Documentation
- Short introduction to the balance code framework.
- List of variables contained in the balance configuration namelist balance_conf.nml.
