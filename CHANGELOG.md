# Changelog

## 2022-01-03 Markus Markl markl@tugraz.at

- Added ramp up namelist variable "faster_ramp_up" that enables faster ramp up if set to 1
- Added ramp up namelist variable "t_max_ramp_up" which is the time value at which the antenna factor should reach the experimental value if faster ramp up is used
- Added option to run KiLCA in "write_kilca" method on Balance.m.
- Changed "reset_profiles" method in Balance.m. Now, the profiles written to /profiles/*.dat are interpolated on the "KiLCA r" grid, which is taken from the hdf5 input file group KiLCA_vac/output/background/R.

## 2022-01-10 Markus Markl markl@tugraz.at

- Added Documentation: short introduction to Balance Suite, Quick Start Guide
- Added README.md in template_scripts

## 2022-01-11 Markus Markl markl@tugraz.at

- Added method plt_antenna_ramp_up to postproc python class, which plots the antenna factor over time, also as percentage of experimental value.
- Changed some names of methods in postproc python class.
