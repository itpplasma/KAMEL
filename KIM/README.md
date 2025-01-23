# KiLCA Integral Model (KIM)
This is the integral plasma response model based on the code and model Kinetic Linear Cylindrical Approximation (KiLCA).

## Configuration
The code is configured with the namelist file */nmls/KIM_config.nml*.

## Compilation
To compile the code:
```
make
```

The build process downloads the Zeal package (complex root finder) which requires a working LAPACK installation.

## KIM dispersion relation in WKB approximation
The dispersion relation for the kinetic integral kernels is determined in Python with the complex root solver module cxroots. The source code is contained in python/WKB-dispersion. The main class and a couple of tutorial functions is contained in wkb.py.