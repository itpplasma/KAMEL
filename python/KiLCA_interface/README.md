# KiLCA_interface

Python interface for KiLCA that simplifies configuration and execution. KiLCA uses plain text configuration files; this interface automates reading and writing these files, allowing parameters to be changed by modifying class variables.

# KiLCA_postprocessor

This is a postprocessing class that handles the output of KiLCA. Most importantly, it reads the `EB.dat` files in `linear-data/m_<m>_n_<n>_flab_[<re>,<im>]/` that contain the electromagnetic perturbation fields of each mode. `EB.dat` has 13 columns: the radius `r` followed by the real and imaginary parts of `Er`, `Etheta`, `Ez`, `Br`, `Btheta` and `Bz` (CGS units; see [Data.md](../../Data.md)). The class also includes plotting methods.
