# KiLCA Integral Model (KIM)
This is the integral plasma response model based on the code and model Kinetic Linear Cylindrical Approximation (KiLCA).

## Dependencies

### Required Libraries
- **LAPACK/BLAS** - Linear algebra
- **fortnum** - Special functions (Bessel) and error function
- **HDF5** - Data storage (with Fortran bindings)
- **SuiteSparse** - Sparse matrix operations (UMFPACK)
- **SuperLU** - Sparse direct solver

### Auto-fetched Dependencies
The following are automatically downloaded during build:
- **fortnum** - Numerical core: complex region-root finder, quadrature, ODE, special functions
- **QUADPACK** - Adaptive quadrature (from Netlib)
- **libcerf** - Complex error function
- **SLATEC** - Mathematical library
- **ddeabm** - ODE solver

### Build Dependencies
- **KiLCA** - Must be built before KIM (provides core library)

## Configuration
The code is configured with the namelist file */nmls/KIM_config.nml*.

### Fokker-Planck equilibrium parallel flow

`Vz_file` is the cylindrical `z` velocity (the toroidal component in KIM's
straight-cylinder model), in cm/s. A finite file must be declared explicitly:

```fortran
&KIM_PROFILES
  Vz_file = 'Vz.dat'
  parallel_flow_convention = 'toroidal'
/
```

KIM treats this as the axial component of a field-aligned ion flow and applies
the cylindrical projection once,
`V_parallel = Vz/h_z`, where
`h_z = sign(B_z)/sqrt(1 + (r/(q R0))**2)`. All configured ion Maxwellians use
that profile; electrons remain stationary. Radial force balance continues to
use the original toroidal component, while Fokker-Planck susceptibilities use

```text
x2(s,m_phi) = -(omega_E + k_parallel V_parallel(s)
                 + m_phi omega_c(s) - omega) / nu(s).
```

No `Vz.dat`, or an identically zero profile, preserves the historical zero-flow
path bit-for-bit. A finite profile with no declared convention or insufficient
radial coverage is rejected before susceptibility assembly. In-memory callers
may pass the explicitly field-parallel `Vpar_in` argument to
`set_profiles_from_arrays`; omission means zero flow. Poloidal neoclassical flow
is not added by this model.

### Fokker-Planck Fourier-kernel FLR model

The two-wavenumber FP kernels take their finite-Larmor-radius arguments from a
named model selected in `&KIM_CONFIG`:

```fortran
&KIM_CONFIG
  fp_ks_model_name = 'kperp_full'   ! or 'kr_only'
/
```

`'kperp_full'` (the default) keeps the helical wavenumber in the arguments,
`b_+ = rho_L^2 (2 k_s^2 + k_r^2 + k_r'^2)/2` and
`b_x = rho_L^2 sqrt(k_s^2+k_r^2) sqrt(k_s^2+k_r'^2)`, as required by the
oracle-verified derivation. `'kr_only'` is the named historical radial-only
approximation (`k_s` dropped), retained for compatibility comparisons with the
old diagonal routine; it is an approximation, not the physics default. Note
that earlier revisions dropped `k_s^2` by default — results at finite `k_s`
change under the new default. The FLR arguments use the stored per-species
Larmor radius, so `ion_flr_scale_factor` applies to this kernel exactly as it
does to the diagonal path.

## Compilation
To compile the code:
```
make
```

The complex region-root finder behind `WKB_dispersion_solver='region_roots'` is provided by fortnum and requires a working LAPACK installation.

## WKB Dispersion Solver

Solves the kinetic dispersion relation D(k_r) = 0 for complex radial wavenumber k_r at each radial grid point.

### Configuration (KIM_config.nml)

```fortran
&WKB_dispersion
  WKB_dispersion_solver = 'Muller'  ! 'Muller' or 'region_roots'
  WKB_dispersion_mode = 'KIM'       ! 'KIM' (full Bessel) or 'FLRE' (finite Larmor radius expansion)
/
```

### Solvers

**Muller** (recommended): Iterative root finder using previous root as initial guess. Reliable branch tracking.

**region_roots**: Contour integration + Newton refinement via fortnum `complex_region_roots`. Per-branch tracking with configurable parameters in the `&wkb_dispersion` namelist:
- `WKB_max_tracked_branches = 4` - Maximum simultaneous branches
- `WKB_branch_search_halfwidth = 1.5` - Search window per branch
- `WKB_broad_search_halfwidth = 5.0` - Initial discovery window
- `WKB_broad_search_interval = 0` - Broad discovery every N grid points; 0 means only when starting or all branches are lost

### Output

Results written to `out/m*_n*/dispersion/` directory:
- `muller_branches_*.dat` or `region_roots_branches_*.dat` - Branch-tracked roots
- HDF5 files with same data for post-processing

### Python Alternative

Python implementation using cxroots: `python/KIMpy/WKB-dispersion/wkb.py`

## References

P. Kravanja, M. Van Barel, O. Ragos, M.N. Vrahatis, F.A. Zafiropoulos,
*ZEAL: A mathematical software package for computing zeros of analytic functions*,
Computer Physics Communications **124** (2000) 212-232.
[doi:10.1016/S0010-4655(99)00429-4](https://doi.org/10.1016/S0010-4655(99)00429-4)
