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
- **ZEAL** - Complex root finder ([Kravanja et al., 2000](#references))
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

The build process downloads the ZEAL package (complex root finder) which requires a working LAPACK installation.

## WKB Dispersion Solver

Solves the kinetic dispersion relation D(k_r) = 0 for complex radial wavenumber k_r at each radial grid point.

### Configuration (KIM_config.nml)

```fortran
&WKB_dispersion
  WKB_dispersion_solver = 'Muller'  ! 'Muller' or 'ZEAL'
  WKB_dispersion_mode = 'KIM'       ! 'KIM' (full Bessel) or 'FLRE' (finite Larmor radius expansion)
/
```

### Solvers

**Muller** (recommended): Iterative root finder using previous root as initial guess. Reliable branch tracking.

**ZEAL**: Contour integration + Newton refinement. Per-branch tracking with configurable parameters in `zeal_input.f90`:
- `MAX_TRACKED_BRANCHES = 4` - Maximum simultaneous branches
- `BRANCH_SEARCH_HALFWIDTH = 1.5` - Search window per branch
- `BROAD_SEARCH_HALFWIDTH = 5.0` - Initial discovery window

### Output

Results written to `out/m*_n*/dispersion/` directory:
- `muller_branches_*.dat` or `zeal_branches_*.dat` - Branch-tracked roots
- HDF5 files with same data for post-processing

### Python Alternative

Python implementation using cxroots: `python/KIMpy/WKB-dispersion/wkb.py`

## References

P. Kravanja, M. Van Barel, O. Ragos, M.N. Vrahatis, F.A. Zafiropoulos,
*ZEAL: A mathematical software package for computing zeros of analytic functions*,
Computer Physics Communications **124** (2000) 212-232.
[doi:10.1016/S0010-4655(99)00429-4](https://doi.org/10.1016/S0010-4655(99)00429-4)