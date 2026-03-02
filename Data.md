# KAMEL Data Flow and Input/Output Documentation

This document describes the actual data pipeline for KAMEL (Kinetic plAsma response ModEL) framework, based on examination of source code and example configurations.

## Overview

KAMEL consists of three main plasma response codes:
1. **KiLCA** - Cylindrical finite Larmor radius (FLR) solver
2. **KIM** - Integral model solver
3. **QL-Balance** - Quasilinear transport evolution code

## 1. KiLCA (Kinetic Linear Cylindrical Approximation)

### 1.1 Input Files

KiLCA uses plain text input files with comment-based configuration (values before `#` are read).

#### Main Configuration Files

Located in the run directory:

| File | Description |
|------|-------------|
| `background.in` | Machine parameters and background plasma settings |
| `antenna.in` | RMP coil antenna configuration |
| `modes.in` | List of (m,n) Fourier modes to compute |
| `output.in` | Output control flags |
| `zone_*_<type>.in` | Zone definitions (flre, vacuum, homogeneous medium) |
| `eigmode.in` | Eigenmode solver settings (if used) |

#### background.in Structure

```
170.69              # Major radius R0 (cm)
70.0                # Plasma minor radius (cm)
23176.46            # Toroidal B-field at axis (Gauss)
../profiles/        # Path to profile files
1                   # Recalculate background flag (1=yes, 0=no)
f                   # Background type: f=full, w=WKB, h=homogeneous
9                   # Spline degree (must be odd)
1.0e9               # V_gal_sys: moving frame velocity (cm/c)
1.0e0               # V_scale: velocity profile scaling
2.0                 # Ion mass (proton mass units)
1.0e-0              # Electron collision coefficient
1.0e-0              # Ion collision coefficient
0                   # Debug flag
```

#### antenna.in Structure

```
90.0             # Antenna radius (cm)
0.0              # Current density layer width (0=delta function)
1.0e13           # Coil current (statamp)
(1.0e0, 0.0e0)   # Complex frequency (1/c) in lab frame
5                # Number of antenna modes from modes.in
1                # Debug flag
0                # Eigenmode problem flag (0=no, 1=yes)
```

#### modes.in Structure

List of (m,n) mode pairs, one per line:
```
(-4,2)
(-5,2)
(-6,2)
...
```

#### output.in Structure

```
2                   # Background data: 0=skip, 1=calc, 2=calc+save
2                   # Linear data: 0=skip, 1=calc, 2=calc+save
2                   # Various quantities: 0=skip, 1=calc, 2=calc+save
0                   # Dispersion: 0=skip, 1=calc, 2=calc+save
0                   # Debug flag
1                   # Current density perturbation
1                   # Absorbed power
1                   # Dissipated power
1                   # Kinetic flux (r-component)
1                   # Poynting flux (r-component)
1                   # Total flux (r-component)
1                   # Number density perturbation
0                   # Lorentz torque
```

#### Profile Files

Located in the directory specified in `background.in` (typically `./profiles/` or `../profiles/`).

**Format**: Two-column ASCII files (`.dat` and `.dim` pairs)
- `.dim` file: Single integer specifying number of data points
- `.dat` file: Two columns (radius, value)

**Required profiles when recalculate flag = 1** (7 files):
| Filename | Description | Units |
|----------|-------------|-------|
| `q.dat`/`q.dim` | Safety factor profile | dimensionless |
| `n.dat`/`n.dim` | Electron density | cm^-3 |
| `Te.dat`/`Te.dim` | Electron temperature | eV |
| `Ti.dat`/`Ti.dim` | Ion temperature | eV |
| `Vth.dat`/`Vth.dim` | Poloidal velocity | cm/c |
| `Vz.dat`/`Vz.dim` | Parallel velocity | cm/c |
| `dPhi0.dat`/`dPhi0.dim` | Electric potential | Gauss units |

**Required profiles when recalculate flag = 0** (11 files):
Above 7 files plus:
| Filename | Description |
|----------|-------------|
| `nui.dat`/`nui.dim` | Ion collision frequency |
| `nue.dat`/`nue.dim` | Electron collision frequency |
| `B0t.dat`/`B0t.dim` | Toroidal magnetic field component |
| `B0z.dat`/`B0z.dim` | Parallel magnetic field component |

### 1.2 Output Files

KiLCA creates output in a directory structure under the run path:

#### Main Output Files

| File | Type | Description |
|------|------|-------------|
| `linear-data/m_<m>_n_<n>_flab_[<re>,<im>]/EB.dat` | Binary | Electromagnetic fields for each (m,n) mode |
| `linear-data/m_<m>_n_<n>_flab_[<re>,<im>]/zone_0_current_dens_p_0_t_lab.dat` | Text | Current density perturbations |
| `formfactors.dat` | Binary (unformatted Fortran) | Transport form factors for QL-Balance |
| `equil_r_q_psi.dat` | Text | Equilibrium profiles on computational grid |

#### EB.dat File Structure

**Format**: Binary/unformatted file written by C++ I/O routines

Contains for each (m,n) mode on (r, theta) grid:
- E_r, E_theta, E_z (complex, Gauss units)
- B_r, B_theta, B_z (complex, Gauss units)

Read using `load_data_file()` from [inout.cpp](KiLCA/io/inout.cpp).

#### formfactors.dat File Structure

**Format**: Fortran unformatted binary

Written by `save_form_facs_fortran()` in [save_fortran.f90](KiLCA/post_processing/save_fortran.f90:3):

```fortran
write(10) modes_dim, dimg, m_min, m_max, n_min, n_max
write(10) sqr_psi_grid(1), sqr_psi_grid(dimg)
write(10) modes_index(m_min:m_max, n_min:n_max)
write(10) transpose(ffpsi(1:dimg, 1:modes_dim))  ! Form factors
```

Used by QL-Balance to compute quasilinear transport coefficients.

#### Additional Output Files

When enabled in `output.in`:
- Various diagnostic quantities in subdirectories
- `label_grid.dat`, `radial_grid.dat` - Grid information
- Power absorption/dissipation profiles
- Flux diagnostics

### 1.3 Data Flow for KiLCA

```
Input Profiles        Configuration Files
(*.dat, *.dim)        (*.in)
     │                     │
     └──────────┬──────────┘
                │
                ▼
         ┌──────────────┐
         │   KiLCA.x    │
         │  Executable  │
         └──────┬───────┘
                │
        ┌───────┴────────┐
        │                │
        ▼                ▼
  EB.dat files    formfactors.dat
  (per mode)      (for QL-Balance)
        │                │
        └───────┬────────┘
                │
                ▼
    Used by QL-Balance or
    post-processing tools
```

---

## 2. KIM (KiLCA Integral Model)

### 2.1 Input Files

KIM uses a single Fortran namelist file for configuration.

#### Configuration File: KIM_config.nml

**Location**: Specified as command-line argument or defaults to `./KIM_config.nml`

**Namelist Groups**:

##### KIM_CONFIG Namelist

```fortran
&kim_config
    number_of_ion_species = 1
    read_species_from_namelist = .false.
    type_of_run = 'electrostatic'          ! 'electrostatic', 'electromagnetic', 'WKB_dispersion', 'flr2_benchmark'
    collision_model = 'FokkerPlanck'    ! or 'Krook'
    artificial_debye_case = 0           ! 0=full, 1=Debye, 2=exclude Debye
    turn_off_ions = .false.
    turn_off_electrons = .false.
    plasma_type = 'D'                   ! 'D' or 'H'
    rescale_density = .false.
    number_density_rescale = 1d-6
/
```

##### KIM_IO Namelist

```fortran
&kim_io
    profile_location = './profiles/'
    output_path = './out/'
    hdf5_input = .false.
    hdf5_output = .false.
    fdebug = 1                          ! 0=off, 1=basic, 2=detailed
    fstatus = 1
    calculate_asymptotics = .true.
/
```

##### KIM_SETUP Namelist

```fortran
&kim_setup
    btor = -17977.413                   ! Toroidal field at axis (Gauss)
    r0 = 165.0                          ! Major radius (cm)
    m_mode = -6                         ! Poloidal mode number
    n_mode = 2                          ! Toroidal mode number
    omega = 0.0                         ! Perturbation frequency (rad/s)
    spline_base = 1                     ! 1=hat functions
    type_br_field = 12                  ! Field type selector
    collisions_off = .false.
    eps_reg = 0.01
    set_profiles_constant = 0
    bc_type = 3                         ! Boundary condition type
    Br_boundary_re = 1.0                ! Real part of Br at right boundary (electromagnetic run type)
    Br_boundary_im = 0.0                ! Imaginary part of Br at right boundary (electromagnetic run type)
/
```

##### KIM_GRID Namelist

```fortran
&kim_grid
    r_min = 3.0                         ! Minimum radius (cm)
    r_plas = 67.0                       ! Plasma radius (cm)
    width_res = 0.5                     ! Resonance layer width factor
    ampl_res = 15.0                     ! Resonance layer amplitude factor
    hrmax_scaling = 1.0
    l_space_dim = 512                   ! Spline grid points
    rg_space_dim = 512                  ! Cell boundary points
    grid_spacing_rg = "equidistant"     ! "equidistant", "non-equidistant", "adaptive"
    grid_spacing_xl = "equidistant"
    theta_integration = "GaussLegendre" ! "GaussLegendre", "RKF45", "QUADPACK"
    rkf45_atol = 1.0d-9                 ! RKF45 absolute tolerance
    rkf45_rtol = 1.0d-6                 ! RKF45 relative tolerance
    quadpack_algorithm = 'QAG'          ! 'QAG' or 'QAGS'
    quadpack_key = 6                    ! Gauss-Kronrod rule (1-6)
    quadpack_limit = 500                ! Max subdivisions
    quadpack_epsabs = 1.0d-10
    quadpack_epsrel = 1.0d-10
    quadpack_use_u_substitution = .true.
    kernel_taper_skip_threshold = 1.0d-6
    Larmor_skip_factor = 5
    gauss_int_nodes_Nx = 31
    gauss_int_nodes_Nxp = 30
    gauss_int_nodes_Ntheta = 17
/
```

##### KIM_SPECIES Namelist

```fortran
&kim_species
    zi = 1     ! Charge numbers (array)
    ai = 2     ! Mass numbers (array)
/
```

#### Profile Files for KIM

Located in directory specified by `profile_location` (default `./profiles/`).

**Format**: Two-column ASCII `.dat` files (no `.dim` files needed)

**Profile Files**:
| Filename | Description | Units |
|----------|-------------|-------|
| `n.dat` | Density profile | cm^-3 |
| `Te.dat` | Electron temperature | eV |
| `Ti.dat` | Ion temperature | eV |
| `Er.dat` | Radial electric field | statvolt/cm |
| `q.dat` | Safety factor | dimensionless |

**File format**: Two columns per line: `radius  value`

Files are read in [species_mod.f90](KIM/src/background_equilibrium/species_mod.f90) using standard Fortran `open()` and `read()`.

### 2.2 Output Files

KIM writes output to `output_path/m<m>_n<n>/` directory structure.

#### Main Output Files

**Text Output Files** (when `hdf5_output=false`):

Created in [IO_collection.f90](KIM/src/util/IO_collection.f90):

| Directory | File | Description |
|-----------|------|-------------|
| `fields/` | `phi_*.dat` | Electrostatic potential profiles |
| `fields/` | `E_r_*.dat` | Radial electric field |
| `kernel/` | `A_sparse_check_*.dat` | Sparse kernel matrix elements |
| `backs/` | `B0z.dat` | Parallel magnetic field |
| `backs/` | `B0th.dat` | Poloidal magnetic field |
| `backs/` | `B0.dat` | Total magnetic field |
| `backs/` | `hz.dat`, `hth.dat` | Magnetic field unit vectors |
| `backs/` | `dpress.dat` | Pressure gradient |
| `backs/` | `p_tot.dat` | Total pressure |
| `grid/` | `rg_xb.dat`, `rg_xc.dat` | Grid boundary and center points |
| `grid/` | `xl_xb.dat`, `xl_xc.dat` | Spline grid points |

**Output format**: ASCII text, typically three columns: `x_coord  real_part  imag_part` for complex quantities.

**HDF5 Output** (when `hdf5_output=true`):
- Not yet fully implemented in current version
- Future: Complete run data in single HDF5 file

### 2.3 Data Flow for KIM

```
Profile Files          KIM_config.nml
(*.dat)               (namelists)
    │                      │
    └──────────┬───────────┘
               │
               ▼
         ┌─────────┐
         │  KIM.x  │
         └────┬────┘
              │
     ┌────────┼────────┐
     ▼        ▼        ▼
  fields/  kernel/  backs/
  (*.dat)  (*.dat)  (*.dat)
     │        │        │
     └────────┼────────┘
              │
              ▼
    Post-processing
   (Python: kim_data.py)
```

---

## 3. QL-Balance (Quasilinear Transport)

### 3.1 Input Files

QL-Balance requires both its own configuration and outputs from KiLCA runs.

#### Configuration File: balance_conf.nml

**Location**: Run directory

**BALANCENML Namelist**:

```fortran
&BALANCENML
    ! Paths to KiLCA runs
    flre_path = '/path/to/kilca/flre/'
    vac_path = '/path/to/kilca/vacuum/'

    ! Machine parameters
    btor = -18009.167                   ! Toroidal field (Gauss)
    rtor = 165                          ! Major radius (cm)
    rmin = 3                            ! Minimum radius (cm)
    rmax = 69.9                         ! Maximum radius (cm)
    rsepar = 64.37086                   ! Separatrix radius (cm)

    ! Grid parameters
    npoimin = 3000                      ! Minimum grid points
    gg_factor = 50                      ! Grid generation factor
    gg_width = 2
    gg_r_res = 95.34                    ! Resonant surface radius

    ! Time evolution
    Nstorage = 1000                     ! Storage array size
    tmax_factor = 50                    ! Max time factor
    stop_time_step = 1e-06              ! Stopping time step
    save_prof_time_step = 10            ! Save every Nth step

    ! Physics parameters
    antenna_factor = 5.4264272          ! RMP amplitude factor
    dperp = 10000                       ! Anomalous diffusion (cm^2/s)
    Z_i = 1                             ! Ion charge
    am = 2                              ! Ion mass number

    ! Boundary conditions
    iboutype = 1
    rb_cut_in = 20                      ! Inner boundary (begin cutoff)
    re_cut_in = 25                      ! Inner boundary (end cutoff)
    rb_cut_out = 67.5                   ! Outer boundary (begin cutoff)
    re_cut_out = 68                     ! Outer boundary (end cutoff)

    ! I/O paths
    path2inp = '/path/to/input.hdf5'   ! Input HDF5 file
    path2out = '/path/to/output.hdf5'  ! Output HDF5 file

    ! Run type and options
    type_of_run = 'SingleStep'          ! 'SingleStep', 'TimeEvolution', 'ParameterScan'
    ihdf5IO = 1
    iwrite = 0                          ! Verbose output flag
    eps = 1e-06                         ! Numerical tolerance
    write_formfactors = .false.
    paramscan = .false.
    diagnostics_output = .false.
    br_stopping = .true.                ! Br stopping criterion
    suppression_mode = .true.
    debug_mode = .false.
    readfromtimestep = 0

    ! Ramp-up mode
    ramp_up_mode = 0                    ! 0=no ramp, 1=ramp up
    t_max_ramp_up = 1e-2
    antenna_max_stopping = 2.0

    ! Additional physics
    temperature_limit = 20.0            ! Lower T limit (eV)
    gyro_current_study = 0
    viscosity_factor = 1
    misalign_diffusion = .false.

    ! Equilibrium data path
    equil_path = '/path/to/equil_r_q_psi.dat'
/
```

#### Required Input Files from KiLCA

From `flre_path` directory:
- `EB.dat` - Electromagnetic fields from plasma response calculation
- `formfactors.dat` - Quasilinear form factors

From `vac_path` directory:
- `EB.dat` - Electromagnetic fields from vacuum calculation

From `equil_path`:
- `equil_r_q_psi.dat` - Equilibrium data (r, q, psi profiles)

#### Input HDF5 File Structure

File specified by `path2inp` contains:
- Initial plasma profiles (ne, Te, Ti, vt)
- Equilibrium data
- Mode information
- Grid specifications

(Created by preprocessing scripts in MATLAB/Python)

### 3.2 Output Files

QL-Balance writes results to HDF5 file specified by `path2out`.

#### Output HDF5 File Structure

**Main output file**: `<shot>_<time>_<name>.hdf5`

**HDF5 Groups and Datasets** (typical structure):

```
/time_evolution/
    t_array                 - Time points array
    ne(r,t)                 - Density evolution
    Te(r,t)                 - Electron temperature evolution
    Ti(r,t)                 - Ion temperature evolution
    vt(r,t)                 - Toroidal velocity evolution

/transport/
    D_ql(r,t)               - Quasilinear diffusion coefficient
    D_anom(r,t)             - Anomalous diffusion

/diagnostics/
    - Various diagnostic quantities per time step

/metadata/
    - Run parameters
    - Git version
    - Timestamp
```

(Note: Exact HDF5 structure depends on `type_of_run` and flags set)

#### Profile Snapshots

When `save_prof_time_step > 0`:
- Profile snapshots saved periodically
- Format depends on output flags

### 3.3 Data Flow for QL-Balance

```
KiLCA outputs           Input HDF5          balance_conf.nml
(EB.dat, formfactors)   (profiles+equil)    (namelist)
        │                      │                  │
        └──────────────────────┼──────────────────┘
                               │
                               ▼
                        ┌──────────────┐
                        │ ql-balance.x │
                        └──────┬───────┘
                               │
                               ▼
                        Output HDF5 file
                     (time evolution data)
                               │
                               ▼
                        Post-processing
                   (balancepost.py, MATLAB)
```

---

## 4. Complete Workflow Pipeline

### 4.1 Typical KAMEL Workflow

```
1. PREPROCESSING
   ├─ Experimental data (GEQDSK, profiles)
   ├─ fouriermodes.x (Fourier preprocessor)
   │  └─ Outputs: amn.dat, equil_r_q_psi.dat
   └─ Profile processor (Python/MATLAB)
      └─ Outputs: Extended profiles, INPUT.hdf5

2. KiLCA CALCULATIONS
   ├─ Vacuum run
   │  ├─ Inputs: background.in, antenna.in, modes.in, profiles/
   │  └─ Outputs: vac/EB.dat
   └─ Plasma (flre) run
      ├─ Inputs: background.in (with plasma), antenna.in, modes.in, profiles/
      └─ Outputs: flre/EB.dat, formfactors.dat

3. QL-BALANCE TRANSPORT
   ├─ Inputs: KiLCA outputs + balance_conf.nml + INPUT.hdf5
   └─ Outputs: balance_out.hdf5 (time evolution)

4. POST-PROCESSING
   └─ Python/MATLAB analysis of HDF5 outputs
```

### 4.2 Alternative: KIM Standalone

```
1. PROFILE PREPARATION
   └─ Create profile .dat files (n, Te, Ti, Er, q)

2. KIM CALCULATION
   ├─ Input: KIM_config.nml + profiles/*.dat
   └─ Outputs: fields/*.dat, kernel/*.dat, backs/*.dat

3. POST-PROCESSING
   └─ Python analysis (kim_data.py)
```

---

## 5. File Format Summary

| Format | Usage | Tools |
|--------|-------|-------|
| `.in` files | KiLCA configuration | Text editor |
| `.nml` files | Fortran namelists (KIM, QL-Balance) | Text editor |
| `.dat` + `.dim` | KiLCA profiles | Custom I/O in inout.cpp |
| `.dat` (2-column) | KIM profiles | Standard Fortran I/O |
| `.hdf5` | QL-Balance I/O, data exchange | h5py, HDF5 tools |
| `EB.dat` | Binary EM fields | Custom binary format |
| `formfactors.dat` | Unformatted Fortran | Fortran unformatted I/O |

---

## 6. Python Interface Notes

### Python Interfaces

Located in `/python/`:
- **KiLCA_interface**: Manages KiLCA input generation and runs
- **KIMpy**: KIM dispersion relation calculations
- **QL_Balance_interface**: QL-Balance workflow management
- **postproc_class**: Post-processing for all codes
- **KIMgui**: Graphical interface for KIM configuration and execution

---

## 7. Units Convention

All codes use **CGS/Gaussian units**:
- Length: cm
- Magnetic field: Gauss
- Density: cm^-3
- Temperature: eV
- Velocity: cm/s or cm/c (velocity/speed of light)
- Electric field: statvolt/cm
- Current: statamp

---

**Last Updated**: 2026-01-17
**Based on**: Actual source code examination of KAMEL repository
