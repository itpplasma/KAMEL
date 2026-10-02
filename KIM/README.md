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
- **libcerf** - Complex error function

### Build Dependencies
- **KiLCA** - Must be built before KIM (provides core library)

## Configuration
The code is configured with the namelist file */nmls/KIM_config.nml*.

For `grid_spacing_rg` or `grid_spacing_xl` set to `adaptive` or
`non-equidistant`, the generated physical nodes include `r_min` and `r_plas`
exactly. Interior nodes retain the Gaussian refinement around the signed
rational surface; the final step is shortened to the requested outer
boundary. A terminal remainder within an accumulated coordinate-roundoff
allowance is merged into that final step; a resolved short cell is retained.
The allowance accounts for the number of coordinate additions, not only the
last addition. It does not bound errors in the nonlinear monitor. A merge
is rejected if its allowance is not small relative to the local step.
Cell centers lie strictly inside their accepted cells.
Earlier adaptive outputs could extend beyond `r_plas`; those outputs used
a different physical boundary and must be interpreted with their saved
coordinates.

The nominal `rg_space_dim` and `l_space_dim` set the spacing scale, rather
than an exact adaptive node count. The two meshes have separate spacing
selectors and share `ampl_res`, `width_res` and `hrmax_scaling`. Invalid
nonfinite controls, nonpositive width or spacing scale, and spacing factors
that cannot stay positive are rejected. Endpoint-preserving refinement
does not by itself establish response convergence. The inherited
gyroaveraged virtual-hat extension at field boundaries is unchanged.

For `type_of_run = 'electrostatic_periodic'`, the HDF5 `/fields` group contains
the total parallel current `jpar`, the electron contribution `jpar_e`, the
summed ion contribution `jpar_i`, and one dataset per configured ion species:
`jpar_i1`, `jpar_i2`, and so on. Ion dataset indices follow the one-based ion
order in `&KIM_species`; consequently, `jpar = jpar_e + jpar_i` and
`jpar_i = sum(jpar_i1, jpar_i2, ...)`.

## FLR2 run type

Set `type_of_run = 'flr2'` to solve the second-order finite-Larmor-radius
response imported from KiLCA-FLR2 through the normal KIM lifecycle:

```fortran
&KIM_CONFIG
  type_of_run = 'flr2'
  collision_model = 'FokkerPlanck'
  number_of_ion_species = 1
/

&KIM_SETUP
  m_mode = -6
  n_mode = 2
  R0 = 165.0
  type_br_field = 12
  Br_boundary_re = 1.0
  Br_boundary_im = 0.0
/
```

This mode does not use KiLCA-FLR2's profile or equilibrium readers. KIM
calculates the radial grid, cylindrical equilibrium, thermodynamic forces,
collision frequencies, Larmor radii, and Fokker–Planck susceptibilities once;
the FLR2 response consumes those prepared arrays directly. The result is
written through the standard KIM field outputs:

- `fields/Br`: radial magnetic-field perturbation in G
- `fields/Phi`: electrostatic potential perturbation in statV
- `fields/jpar`: parallel current density in statA/cm²

`type_br_field = 12` applies the constant complex field configured by
`Br_boundary_re` and `Br_boundary_im`. `type_br_field = 11` reads
`./inp/Br_in.dat`. The current implementation requires Fokker–Planck
collisions, `omega = 0`, and exactly one positively charged ion species.

The optional `&KIM_FLR2` group enables or disables electron and ion FLR,
potential, and current terms. All switches default to `.true.`; see the
[namelist reference](nmls/README.md).

The existing `type_of_run = 'flr2_benchmark'` is different: it exercises the
FLR2 approximation of KIM's non-local integral kernel and global Poisson
solver.

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

Python implementation using cxroots: `python/KIMpy/wkb.py`

## References

P. Kravanja, M. Van Barel, O. Ragos, M.N. Vrahatis, F.A. Zafiropoulos,
*ZEAL: A mathematical software package for computing zeros of analytic functions*,
Computer Physics Communications **124** (2000) 212-232.
[doi:10.1016/S0010-4655(99)00429-4](https://doi.org/10.1016/S0010-4655(99)00429-4)
