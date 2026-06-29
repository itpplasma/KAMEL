# KiLCA Integral Model (KIM)
This is the integral plasma response model based on the code and model Kinetic Linear Cylindrical Approximation (KiLCA).

## Dependencies

### Required Libraries
- **LAPACK/BLAS** - Linear algebra
- **GSL** - GNU Scientific Library
- **HDF5** - Data storage (with Fortran bindings)
- **SuiteSparse** - Sparse matrix operations (UMFPACK)
- **SuperLU** - Sparse direct solver

### Auto-fetched Dependencies
The following are automatically downloaded during build:
- **QUADPACK** - Adaptive quadrature (from Netlib)
- **libcerf** - Complex error function
- **SLATEC** - Mathematical library
- **ddeabm** - ODE solver

### Build Dependencies
- **KiLCA** - Must be built before KIM (provides core library)

## Configuration
The code is configured with the namelist file */nmls/KIM_config.nml*.

## Compilation
To compile the code:
```
make
```

## WKB Dispersion Solver

Solves the kinetic dispersion relation D(k_r) = 0 for complex radial wavenumber k_r at each radial grid point.

### Configuration (KIM_config.nml)

```fortran
&WKB_dispersion
  WKB_dispersion_solver = 'Muller'  ! 'Muller'
  WKB_dispersion_mode = 'KIM'       ! 'KIM' (full Bessel) or 'FLRE' (finite Larmor radius expansion)
/
```

### Solvers

**Muller**: Iterative root finder using previous root as initial guess. Reliable branch tracking.

### Output

Results written to `out/m*_n*/dispersion/` directory:
- `muller_branches_*.dat` - Branch-tracked roots
- HDF5 files with same data for post-processing

### Python Alternative

Python implementation using cxroots: `python/KIMpy/WKB-dispersion/wkb.py`
