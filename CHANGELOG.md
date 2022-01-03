# Changelog

2021-01-03 Markus Markl markl@tugraz.at

	* Added ramp up namelist variable "faster_ramp_up" that enables faster ramp up if set to 1
	* Added ramp up namelist variable "t_max_ramp_up" which is the time value at which the antenna factor should reach the experimental value if faster ramp up is used

	* Added option to run KiLCA in "write_kilca" method on Balance.m.
	* Changed "reset_profiles" method in Balance.m. Now, the profiles written to /profiles/*.dat are interpolated on the "KiLCA r" grid, which is taken from the hdf5 input file group KiLCA_vac/output/background/R.
