# KiLCA Integral Model (KIM)
This is the integral plasma response model based on the code and model Kinetic Linear Cylindrical Approximation (KiLCA).


## Configuration
The code is configured with the namelist file */nmls/KIM_config.nml*. The various variables are:

**Namelist KIM_CONFIG**
- profile_location ... string, path to the profiles
- output_path ... string, path to the output directory
- hdf5_input ... boolean, if true reads from a hdf5 file
- hdf5_output ... boolean, if true writes to a hdf5 output file
- fdebug ... integer, flag for debugging
- fstatus ... integer, flag for giving status updates during the execution of the code
- ispecies ... integer, number of ion species

** Namelist KIM_SETUP **
- btor ... double precision, toroidal magnetic field at the magnetic axis in Gauss
- R0 ... double precision, major radius of the magnetic axis in cm
- m_mode ... integer, poloidal mode number of the resonant surface
- n_mode ... integer, toroidal mode number of the resonant surface
- Zi ... integer, ion charge number; if ispecies > 1 write multiple separated by spaces, i.e. Zi = 1 2 4
- Ai ... integer, ion mass number; if ispeceis > 1 write multiple separated by spaces, i.e. Ai = 2 4 8
- k_space_dim ... integer, dimension of the k space
- reduce_r ... boolean, if true reduces r dim of the input profiles
- reduced_r_dim ... integer, dimension of reduced r dim
- omega ... double precision, perturbation frequency
- 
## Compilation
To compile the code:
```
mkdir build
cd build
cmake ..
make
```