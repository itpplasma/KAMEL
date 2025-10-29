# Fourier Pre-Processor

This directory contains the Fortran pre-processing utilities that convert an
axisymmetric equilibrium into the effective radius, magnetic fluxes, etc.  The build generates the executable
`fouriermodes.x`, which

1. reconstructs the magnetic axis and separatrix from the supplied equilibrium,
2. traces magnetic field lines to tabulate safety factor, flux, and Boozer
   angle information, and
3. projects the perturbation magnetic field into toroidal/poloidal Fourier
   coefficients that are emitted in machine-readable (`*.dat`) form.

The implementation relies on several modules from the open-source
[`libneo`](https://github.com/itpplasma/libneo) project; the required source
files are downloaded automatically during CMake configuration.

## Building

```bash
make
```

This yields `build/fouriermodes.x` and the static helper library
`build/lib/libfouriermodes.a`.  All modules write their unformatted Fortran
records in the current working directory, so run the executable from a folder
that contains the input files described below.

## Required input files

Two plain-text control files configure the run.  Representative templates live
under `inp/`.

### `field_divB0.inp`

Lines are read sequentially (free format) by `field_divB0.f90`:

1. `ipert` – perturbation mode selector  
   `0` (equilibrium only), `1` (vacuum coils), `2` (vacuum + plasma shielding,
   field only), `3` (vacuum + plasma shielding with derivatives – slowest).
2. `iequil` – whether to add the background equilibrium to the output (`0`
   disables it, which is useful if only the perturbation is needed).
3. `ampl` – scalar amplitude applied to all perturbation harmonics.
4. `ntor` – number of toroidal harmonics expected in the perturbation data.
5. `cutoff` – inner radial cutoff in units of ψ/ψₐ used when stitching the
   perturbation solution.
6. `icftype` – coil file format selector (see your dataset documentation).
7. `gfile` – path to the axisymmetric equilibrium (GEQDSK or compatible).
8. `pfile` – path to the coil data file (may be blank when `ipert=0`).
9. `convexfile` – convex hull description used by the coordinate stretching
   routines.
10. `fluxdatapath` – directory containing flux-coordinate data required when
    `ipert ≥ 2`.
11. `nwindow_r` – window size (in radial grid points) for smoothing ψ(R).
12. `nwindow_z` – window size (in vertical grid points) for smoothing ψ(Z).

Leave the string entries quoted; set them to empty quotes (`''`) when the data
are not needed for the selected `ipert`/`iequil` combination.  The default
values in `inp/template_field_divB0.inp` illustrate the expected formatting and
include inline comments describing each switch.

### `fouriermodes.inp`

The main driver reads six integers that control the resolution of the field-line
tracing and the radial/poloidal grids:

1. `mpol` – number of poloidal Fourier modes (|m| ≤ `mpol`).
2. `nstep` – number of integration steps per field-line circuit.
3. `nsqpsi` – radial grid size for 1-D profiles (uniform in √ψ).
4. `nlabel` – number of surfaces retained for 2-D/3-D quantities.
5. `ntheta` – number of points in the Boozer poloidal angle grid.
6. `nsurfmax` – trial surfaces launched when searching for the separatrix.

Adjust these values to balance accuracy against runtime; the defaults in
`inp/fouriermodes.inp` target high-resolution studies.

## Running

After tailoring the input files (and placing all referenced equilibrium/coil
data in accessible locations), execute

```bash
./fouriermodes.x
```
In a suitable directory. Soft-linking the executable in the run directory is usually the easiest way.

`fouriermodes.x` produces several diagnostic and data products in the working
directory:

| File | Description |
| --- | --- |
| `box_size.dat` | R–Z bounds of the computational box inferred from the equilibrium grid. |
| `axis.dat` | Magnetic-axis position (`raxis`, `zaxis`). |
| `btor_rbig.dat` | Toroidal-field strength and magnetic-axis major radius. |
| `separ.dat` | ψ values and R–Z curve describing the last closed flux surface. |
| `phinorm_arr.dat` | Tabulation of normalized toroidal flux vs. ψ on an evenly spaced √ψ grid. |
| `equil_r_q_psi.dat` | Derived 1-D profiles: minor radius, safety factor q, ψ, toroidal flux, volume, etc. |
| `theta_of_theta_qt_flabel.dat` (unformatted) | Mapping between geometric and Boozer poloidal angles on the labeled flux surfaces; header contains grid metadata. |
| `thetabooz.dat` | ASCII version of the previous mapping (one surface per line). |
| `amn.dat` (unformatted) | Complex Fourier coefficients of the vector potential components (A<sub>ψ</sub>, A<sub>θ</sub>) for every (m,n) mode retained. |

The unformatted files are written using Fortran stream I/O; consume them with
matching Fortran or Python readers that understand the binary record layout.

## Notes

- CMake re-downloads the required `libneo` sources whenever you delete
  `src/*_libneo.f90`; avoid running `make clean` if you are offline.
- All computations assume right-handed cylindrical coordinates (R, φ, Z) and
  the Boozer-angle conventions implemented in `libneo`.  Ensure your inputs are
  consistent with those assumptions.
- When only the background equilibrium is needed, set `ipert = 0` and leave the
  coil and flux-coordinate paths blank; the code will skip the perturbation
  branches automatically.

Refer to the inline comments in `template_field_divB0.inp` and
`fouriermodes.f90` for additional implementation details.
