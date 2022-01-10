# Template Scripts
Contains matlab scripts that can be used as templates for certain balance code runs.

## Content
- script_prerun.m is the Matlab script that does the prerun. It creates an "input" hdf5 file for the balance code. This input file is the same for all types of runs (time evolution, parameter scan, ...). It has to be created only once.
- script_prerun_multiple.m runs the prerun for multiple time slices for a single shot.
- linearrun/ contains scripts to run the balance code in linear mode. This is used to determine the quasilinear diffusion coefficients for the full RMP coil current. There are templates for single time/shot combination and multiple times/single shot combination.
- timeevol/ contains scripts to do time evolution runs of the balance code.
- parameterscan/ contains scripts to do parameter scans with the balance code.
