# KAMEL - Kinetic plAsma response ModEL
This repository contains the kinetic plasma response framework containing the linear plasma response codes KiLCA and KIM, as well as the quasilinear transport code QL-Balance. KiLCA and KIM are cylindrical linear plasma response solvers based on a finite Larmor radius and an integral formalism, respectively. QL-Balance is a quasilinear 1D radial transport code. In combination, they are used to model the plasma response to external magnetic perturbation in toroidally confined fusion plasmas. 

Note, when using the template scripts, make sure to change the code, most importantly the paths, according to your project before using it.

## Dependencies

### Core Dependencies

The KAMEL framework requires the following external libraries and tools:

#### Build Tools
- **CMake** (>= 3.24): Modern build system generator
- **Make**: Build automation tool
- **C/C++ Compiler**: 
  - Apple Silicon: clang 16.0 / clang++
  - Linux: gcc/g++ 12.2.0
- **Fortran Compiler**: 
  - Apple Silicon: gfortran 14.2
  - Linux: gfortran 12.2.0
- **MPI**: Message Passing Interface (mpif90 wrapper required)

#### Mathematical Libraries
- **LAPACK**: Linear Algebra Package (for matrix operations)
- **BLAS**: Basic Linear Algebra Subprograms
- **GSL** (2.4): GNU Scientific Library (for special functions, integration, etc.)
- **SuiteSparse**: Suite of sparse matrix algorithms (specifically UMFPACK for sparse solvers)
- **SUNDIALS**: SUite of Nonlinear and DIfferential/ALgebraic equation Solvers (for ODE/DAE integration)

#### Data Handling
- **HDF5**: Hierarchical Data Format (for data storage and exchange)
- **NetCDF**: Network Common Data Form (optional, for certain data formats)

#### Additional Libraries
- **OpenMP**: For parallel computing support
- **FFTW3**: Fastest Fourier Transform in the West (optional, for spectral methods)
- **Zeal**: (Internal dependency)

### Python Dependencies

For the Python interface (KAMELpy), the following packages are required:

```
numpy       # Numerical arrays and operations
scipy       # Scientific computing tools
matplotlib  # Plotting and visualization
h5py        # HDF5 file interface
f90nml      # Fortran namelist file parser
```

Install Python dependencies via:
```bash
cd python
pip install -r requirements.txt
# or
pip install .
```

### MATLAB Dependencies

The MATLAB interface requires:
- MATLAB R2020a or newer
- No additional toolboxes required for basic functionality

### External Libraries (Bundled)

The following libraries are included in the `external/` directory:
- GSL 2.4
- LAPACK 3.2.1
- SuiteSparse (including UMFPACK, AMD, CHOLMOD, etc.)
- SUNDIALS

### Platform-Specific Notes

#### macOS (Apple Silicon)
- Accelerate framework is used for optimized BLAS/LAPACK
- Homebrew recommended for package management:
  ```bash
  brew install cmake gfortran open-mpi gsl hdf5
  ```

#### Linux
- System packages typically required:
  ```bash
  # Debian/Ubuntu
  sudo apt-get install cmake gfortran libopenmpi-dev libgsl-dev \
                       liblapack-dev libblas-dev libhdf5-dev
  
  # Red Hat/CentOS
  sudo yum install cmake gcc-gfortran openmpi-devel gsl-devel \
                   lapack-devel blas-devel hdf5-devel
  ```

## Compilation
For compilation of the codes invoke `make` in this directory.  
Generally, for Apple Silicon the clang/gfortran (version 16.0 and 14.2., respectively) compiler combination is tested. On debian, the gnu compiler version 12.2.0 is tested.

The build system will automatically fetch and build missing dependencies from the `external/` directory if system versions are not found.

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
