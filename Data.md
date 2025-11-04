# KAMEL Data Flow and Input/Output Documentation

This document describes the complete data pipeline for the KAMEL (Kinetic plAsma response ModEL) framework, including all input and output files for the three main codes: KiLCA, KIM, and QL-Balance.

## Overview

KAMEL uses **HDF5 format** as the standard for data exchange between all components. This ensures:
- Standardized data structure across workflow stages
- Metadata tracking (git version, timestamps) for reproducibility
- Efficient storage and access for large datasets

## Data Flow Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         PREPROCESSING STAGE                              │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                    ┌───────────────┼───────────────┐
                    ▼               ▼               ▼
            ┌──────────────┐ ┌──────────┐ ┌──────────────────┐
            │ Experimental │ │ Magnetic │ │  RMP Coil Data   │
            │   Profiles   │ │   Equil  │ │  (field.dat)     │
            │  (n,Te,Ti,v) │ │ (GEQDSK) │ │                  │
            | fun. of flux | |          | |   (optional)     |
            └──────┬───────┘ └────┬─────┘ └────────┬─────────┘
                   │              │                 │
                   └──────────────┼─────────────────┘
                                  ▼
                    ┌─────────────────────────────┐
                    │  fouriermodes.x             │
                    │  (Fourier Preprocessor)     │
                    │  - Flux coordinates         │
                    │  - Safety factor q(r)       │
                    │  - Fourier modes of Br      │
                    |    (if fields.dat given)    |
                    └──────────────┬──────────────┘
                                   │
                    ┌──────────────┼──────────────┐
                    │              │              │
                    ▼              ▼              ▼
              amn.dat      equil_r_q_psi.dat   theta_*.dat
                    │              │              │
                    └──────────────┼──────────────┘
                                   │
                                   ▼
                    ┌─────────────────────────────┐
                    │  Profile Processor          │
                    │  (Python/MATLAB)            │
                    │  - Extend profiles          │
                    │  - Map to effective radius  │
                    │  - Calculate derivatives    │
                    └──────────────┬──────────────┘
                                   │
                                   ▼
                         INPUT_[shot]_[time].hdf5
                                   │
┌──────────────────────────────────┼──────────────────────────────────────┐
│                        MAIN COMPUTATION STAGE                            │
└──────────────────────────────────────────────────────────────────────────┘
                                   │
                    ┌──────────────┼──────────────┐
                    ▼              ▼              ▼
        ┌─────────────────┐  ┌──────────┐  ┌─────────────────┐
        │  KiLCA          │  │   KIM    │  │  QL-Balance     │
        │  (FLR solver)   │  │ (Integ.) │  │  (Transport)    │
        └────────┬────────┘  └────┬─────┘  └────────┬────────┘
                 │                │                  │
                 ▼                ▼                  ▼
            EB.dat           kim_out.hdf5    balance_out.hdf5
         (EM fields)         (EM fields)     (Time evolution)
                 │                │                  │
┌────────────────┼────────────────┼──────────────────┼───────────────────┐
│                         POST-PROCESSING STAGE                           │
└─────────────────────────────────────────────────────────────────────────┘
                 │                │                  │
                 └────────────────┼──────────────────┘
                                  ▼
                    ┌─────────────────────────────┐
                    │  Python/MATLAB Analysis     │
                    │  - KiLCA_postprocessor      │
                    │  - kim_data                 │
                    │  - balancepost              │
                    └──────────────┬──────────────┘
                                   │
                                   ▼
                    Plots, Reports, Physical Quantities
```

## 1. Preprocessing Stage

### 1.1 Fourier Modes Preprocessor (fouriermodes.x)

**Location**: `/PreProc/fourier/`

**Purpose**: Converts axisymmetric equilibrium and coil data into Fourier modes and flux coordinates for KiLCA/QL-Balance.

#### Input Files

| File | Type | Description |
|------|------|-------------|
| `field_divB0.inp` | Text | Configuration for perturbation calculation (mode selector, amplitude, etc.) |
| `fouriermodes.inp` | Text | Grid resolution parameters (mpol, nstep, nsqpsi, nlabel, ntheta) |
| `g[shot].[time]` | Binary | GEQDSK equilibrium file (magnetic flux surfaces) |
| `field.dat` | Text/Binary | RMP coil geometry and current distribution, only used when vacuum perturbation is to be calculated. Can be omitted otherwise |

#### Output Files

| File | Type | Description |
|------|------|-------------|
| `amn.dat` | Binary (unformatted) | Complex Fourier coefficients A_psi, A_theta for (m,n) modes (if vacuum perturbation is calculated) |
| `equil_r_q_psi.dat` | Text | 1D radial profiles: minor radius r, safety factor q, flux psi, volume. Mapping of effective radius to flux labels. |
| `theta_of_theta_qt_flabel.dat` | Binary | Boozer angle mapping (geometric ↔ Boozer coordinates) |
| `thetabooz.dat` | Text | ASCII version of Boozer angle mapping |
| `phinorm_arr.dat` | Text | Normalized toroidal flux vs. poloidal flux |
| `axis.dat` | Text | Magnetic axis position (R_axis, Z_axis) |
| `separ.dat` | Text | Last closed flux surface (LCFS) coordinates |
| `btor_rbig.dat` | Text | Toroidal field strength and major radius |
| `box_size.dat` | Text | R-Z computational domain bounds |

### 1.2 Profile Processor

**Location**: `/python/profile_processor/`, `/utility_scripts/`

**Purpose**: Prepare experimental plasma profiles for KiLCA/QL-Balance input.

#### Input Files

| File | Format | Description |
|------|--------|-------------|
| `[shot].[time]_ne*.dat` | Text | Electron density profile ne(r) [cm^-3] |
| `[shot].[time]_Te*.dat` | Text | Electron temperature profile Te(r) [eV] |
| `[shot].[time]_Ti*.dat` | Text | Ion temperature profile Ti(r) [eV] |
| `[shot].[time]_vt*.dat` | Text | Toroidal rotation profile vt(r) [cm/s] |

#### Output Files

| File | Format | Description |
|------|--------|-------------|
| `INPUT_[shot]_[time].hdf5` | HDF5 | Complete input dataset with profiles, equilibrium, modes |

**HDF5 Structure** (INPUT file, **not yet implemented this way**): (TODO)
```
/profiles/
    ne(r)           - Electron density [cm^-3]
    Te(r)           - Electron temperature [eV]
    Ti(r)           - Ion temperature [eV]
    vt(r)           - Toroidal velocity [cm/s]
    dne_dr(r)       - Density gradient
    dTe_dr(r)       - Electron temperature gradient
    dTi_dr(r)       - Ion temperature gradient
    dvt_dr(r)       - Velocity gradient
/equilibrium/
    r_grid          - Minor radius grid [cm]
    q(r)            - Safety factor profile
    psi(r)          - Poloidal flux [Wb]
    Btor            - Toroidal field [Gauss]
    R0              - Major radius [cm]
/modes/
    m_modes         - Poloidal mode numbers (array)
    n_modes         - Toroidal mode numbers (array)
/metadata/
    shot            - Shot number
    time            - Time slice [ms]
    git_version     - Git hash for reproducibility
    timestamp       - Creation timestamp
```

## 2. Main Computation Stage

### 2.1 KiLCA (Finite Larmor Radius Solver)

**Executable**: `build/KiLCA/KiLCA.x`

**Purpose**: Solves plasma response using finite Larmor radius (FLR) formalism in cylindrical geometry.

#### Configuration Files

| File | Type | Description |
|------|------|-------------|
| `antenna.nml` | Namelist | Antenna configuration (radius, modes, boundary conditions) |
| `background.nml` | Namelist | Background plasma parameters (Btor, R0, ion mass, profiles) |
| `eigmode.nml` | Namelist | Eigenmode search settings (frequency range, convergence) |
| `output.nml` | Namelist | Output control flags |

**Blueprint Location**: `/python/KiLCA_interface/blueprints/`

#### Input Files

| File | Format | Description |
|------|--------|-------------|
| `profiles/ne.dat` | Text | Electron density profile |
| `profiles/Te.dat` | Text | Electron temperature profile |
| `profiles/Ti.dat` | Text | Ion temperature profile |
| `profiles/Vz.dat` | Text | Parallel velocity profile |
| `profiles/Er.dat` | Text | Radial electric field profile |
| `equil_r_q_psi.dat` | Text | Equilibrium data from fouriermodes |

#### Output Files

| File | Format | Description |
|------|------|-------------|
| `EB.dat` | Binary | Electromagnetic fields: E_r, E_theta, E_z, B_r, B_theta, B_z (r,theta grid) |

**EB.dat Structure**:
- Header: nr (radial points), ntheta (poloidal points), nmodes
- For each (m,n) mode:
  - Complex E_r(r,theta), E_theta(r,theta), E_z(r,theta)
  - Complex B_r(r,theta), B_theta(r,theta), B_z(r,theta)

### 2.2 KIM (KiLCA Integral Model)

**Executable**: `build/KIM/KIM.x`

**Purpose**: Solves plasma response using integral formalism with non-local kinetic effects.

#### Configuration File

| File | Type | Description |
|------|------|-------------|
| `KIM_config.nml` | Namelist | Complete KIM configuration (all namelists below) |

**Location**: `/KIM/nmls/KIM_config.nml`

**Namelist Groups**:

1. `kim_config`: Run type and collision model
   - `number_of_ion_species`: Number of ion species (default: 1)
   - `type_of_run`: 'electrostatic' or 'electromagnetic'
   - `collision_model`: 'Krook' or 'FokkerPlanck'
   - `artificial_debye_case`: 0 (full), 1 (Debye), 2 (exclude Debye)

2. `kim_io`: I/O paths and flags
   - `profile_location`: Path to profile .dat files
   - `output_path`: Output directory
   - `hdf5_input`: Use HDF5 input (default: false)
   - `hdf5_output`: Use HDF5 output (default: false)
   - `fdebug`: Debug level (0-2)
   - `calculate_asymptotics`: Enable WKB calculations

3. `kim_setup`: Physical parameters
   - `btor`: Toroidal field at axis [Gauss]
   - `r0`: Major radius [cm]
   - `m_mode`: Poloidal mode number
   - `n_mode`: Toroidal mode number
   - `omega`: Perturbation frequency [rad/s]
   - `type_br_field`: Br field type selector

4. `kim_grid`: Grid resolution and integration
   - `r_plas`: Plasma minor radius [cm]
   - `r_min`: Minimum computational radius [cm]
   - `l_space_dim`: Number of spline grid points (~512)
   - `rg_space_dim`: Number of cell boundaries (~512)
   - `grid_spacing_rg`: "equidistant", "non-equidistant", or "adaptive"
   - `theta_integration`: "GaussLegendre", "RKF45", or "QUADPACK"
   - `gauss_int_nodes_Nx`: x-integration nodes (default: 31)
   - `gauss_int_nodes_Ntheta`: theta nodes (default: 17)

5. `kim_species`: Species properties
   - `zi`: Charge numbers (array)
   - `ai`: Mass numbers (array)

#### Input Files

| File | Format | Description |
|------|------|-------------|
| `profiles/ne.dat` | Text | Electron density [cm^-3] |
| `profiles/Te.dat` | Text | Electron temperature [eV] |
| `profiles/Ti.dat` | Text | Ion temperature [eV] |
| `profiles/q.dat` | Text | Safety factor profile |
| `profiles/dn_dr.dat` | Text | Density gradient |
| `profiles/dTe_dr.dat` | Text | Temperature gradients |

#### Output Files

| File | Format | Description |
|------|------|-------------|
| `kim_out.txt` | Text | Main output: dispersion relation roots, growth rates, frequencies |
| `kim_out.hdf5` | HDF5 | (if hdf5_output=true) Complete output dataset |
| `kim_kernel.dat` | Binary | Non-local kernel matrix elements |
| `kim_phi.dat` | Text | Electrostatic potential phi(r) |
| `kim_asymptotic.dat` | Text | WKB/FLR2 asymptotic quantities |
| `kim_log.txt` | Text | Execution log |

**kim_out.hdf5 Structure** (when enabled):
```
/dispersion/
    omega_real      - Real part of frequency [rad/s]
    omega_imag      - Imaginary part (growth rate) [rad/s]
    k_r             - Radial wave vector [cm^-1]
/fields/
    phi(r)          - Electrostatic potential
    E_r(r)          - Radial electric field
/kernel/
    K_matrix        - Non-local kernel matrix
/asymptotics/
    phi_FLR2(r)     - FLR2 asymptotic potential
    k_fourier       - Fourier space kernel
```

### 2.3 QL-Balance (Quasilinear Transport)

**Executable**: `build/QL-Balance/ql-balance.x`

**Purpose**: Time evolution of plasma profiles including quasilinear and anomalous transport.

#### Configuration File

| File | Type | Description |
|------|------|-------------|
| `balance_conf.nml` | Namelist | Complete QL-Balance configuration |

**Location**: `run_directory/balance_conf.nml`

**Key Parameters**:
- `flre_path`: Path to KiLCA flre run directory
- `vac_path`: Path to KiLCA vacuum run directory
- `path2inp`: Path to input HDF5 file
- `path2out`: Path to output HDF5 file
- `type_of_run`: 'SingleStep', 'TimeEvolution', or 'ParameterScan'
- `btor`: Toroidal field [Gauss]
- `rtor`: Major radius [cm]
- `rmin`, `rmax`: Radial computational domain [cm]
- `antenna_factor`: RMP coil current amplitude factor
- `tmax_factor`: Maximum simulation time factor
- `dperp`: Anomalous diffusion coefficient [cm^2/s]
- `Z_i`: Ion charge
- `am`: Ion mass number
- `stop_time_step`: Time step for stopping criterion
- `save_prof_time_step`: Save profiles every Nth step
- `br_stopping`: Use Br stopping criterion (true/false)
- `suppression_mode`: Track RMP suppression (true/false)

#### Input Files

| File | Format | Description |
|------|------|-------------|
| `INPUT_*.hdf5` | HDF5 | Complete input from preprocessing |
| `flre/EB.dat` | Binary | EM fields from KiLCA flre run |
| `vacuum/EB.dat` | Binary | EM fields from KiLCA vacuum run |
| `equil_r_q_psi.dat` | Text | Equilibrium data |

#### Output Files

| File | Format | Description |
|------|------|-------------|
| `balance_out.hdf5` | HDF5 | Complete time evolution dataset |
| `profiles_*.dat` | Text | Profile snapshots at saved time steps |
| `diagnostics.dat` | Text | Time series of global quantities |
| `balance_log.txt` | Text | Execution log |

**balance_out.hdf5 Structure**:
```
/time_evolution/
    t_array                 - Time points [s]
    ne(r,t)                 - Electron density evolution
    Te(r,t)                 - Electron temperature evolution
    Ti(r,t)                 - Ion temperature evolution
    vt(r,t)                 - Rotation evolution
/transport_coefficients/
    D_ql(r,t)               - Quasilinear diffusion [cm^2/s]
    D_anom(r,t)             - Anomalous diffusion [cm^2/s]
    chi_i(r,t)              - Ion thermal diffusivity [cm^2/s]
    chi_e(r,t)              - Electron thermal diffusivity [cm^2/s]
/diagnostics/
    energy_confinement(t)   - Energy confinement time [s]
    particle_flux(r,t)      - Particle flux [cm^-2 s^-1]
    heat_flux(r,t)          - Heat flux [erg cm^-2 s^-1]
    antenna_amplitude(t)    - RMP amplitude vs time
/resonances/
    r_res(m,n)              - Resonant surface locations [cm]
    island_width(m,n,t)     - Magnetic island widths [cm]
```

## 3. Workflow Templates

**Location**: `/template_scripts/`

### 3.1 Prerun (Data Preparation)

**Script**: `script_prerun.m` (MATLAB) or Python equivalent

**Purpose**: Creates complete HDF5 input file for QL-Balance.

**Workflow**:
1. Read experimental profiles (ne, Te, Ti, vt)
2. Read equilibrium (GEQDSK)
3. Read coil data
4. Run fouriermodes.x to generate Fourier modes
5. Run profile processor to extend/map profiles
6. Run KiLCA vacuum calculation
7. Run KiLCA flre calculation
8. Package everything into `INPUT_[shot]_[time].hdf5`

**Output**: Single HDF5 file ready for QL-Balance

### 3.2 Linear Run

**Script**: `linearrun/linear_balance_single.m`

**Purpose**: Calculate quasilinear diffusion coefficients for full RMP amplitude.

**Type of Run**: `type_of_run = 'SingleStep'`

**Output**: Diffusion coefficients D_ql(r) for given RMP amplitude

### 3.3 Time Evolution Run

**Script**: `timeevol/timeevol_balance_run.m`

**Purpose**: Simulate profile evolution with time-dependent transport.

**Type of Run**: `type_of_run = 'TimeEvolution'`

**Output**: Full time history of profiles and transport coefficients

### 3.4 Parameter Scan

**Script**: `parameterscan/balancerun_parscan.m`

**Purpose**: Systematic scan over parameter space (e.g., RMP amplitude, rotation).

**Type of Run**: `type_of_run = 'ParameterScan'`

**Output**: Multiple runs with varying parameters

## 4. Python/MATLAB Interfaces

### 4.1 KiLCA Interface (Python)

**Location**: `/python/KiLCA_interface/`

**Main Class**: `KiLCA_interface`

**Workflow Example**:
```python
from KiLCA_interface import KiLCA_interface

# Initialize
kilca = KiLCA_interface(shot=39711, time=2000,
                        path='/path/to/run/',
                        rtype='flre', machine='AUG')

# Configure
kilca.set_modes(m=[5,6,7], n=[2,2,2])
kilca.set_ASDEX(nmodes=3)
kilca.set_profiles(ne_file, Te_file, Ti_file, vt_file)

# Run
kilca.write()  # Create input files
kilca.run()    # Execute KiLCA

# Post-process
from KiLCA_interface.KiLCA_postprocessor import KiLCA_postprocessor
post = KiLCA_postprocessor(kilca.path_of_run)
post.read_EB_dat()
post.plot_fields()
```

### 4.2 QL-Balance Interface (Python)

**Location**: `/python/balance_interface/`

**Main Class**: `QL_Balance_interface`

**Workflow Example**:
```python
from balance_interface import QL_Balance_interface

# Initialize
balance = QL_Balance_interface(run_path='/path/to/run/',
                                shot=39711, time=2000,
                                name='timeevol_run',
                                input_file='INPUT_39711_2000.hdf5')

# Configure
balance.set_type_of_run('TimeEvolution')
balance.set_modes(m=[5,6,7], n=[2,2,2])
balance.prepare_balance(Btor=-18000, a_minor=67.0)

# Run
balance.run_balance()

# Post-process
from postproc_class import balancepost
post = balancepost.BalancePost(balance.output_h5_file)
post.plot_profiles_evolution()
post.calculate_transport_coefficients()
```

### 4.3 KIM Interface (Python)

**Location**: `/python/KIMpy/`

**Main Class**: `KIMpy`

**Workflow Example**:
```python
from KIMpy import kimpy

# Initialize
kim = kimpy(config_file='KIM_config.nml')

# Set profiles
kim.set_profiles(ne_file, Te_file, Ti_file, q_file)

# Run
kim.run()

# Analyze dispersion relation
from KIMpy import kim_data
data = kim_data.KIM_data('kim_out.hdf5')
data.plot_dispersion_relation()
data.get_growth_rate()
```

### 4.4 MATLAB Balance Class

**Location**: `/matlab/balance/Balance.m`

**Usage**: See `script_prerun.m` example (lines 73-81)

## 5. Data Format Specifications

### 5.1 Profile File Format (.dat)

**Type**: ASCII text, two-column

**Format**:
```
# Header line (optional, starts with #)
r1  value1
r2  value2
...
rN  valueN
```

**Units**:
- Radial coordinate: cm
- Density: cm^-3
- Temperature: eV
- Velocity: cm/s
- Electric field: statvolt/cm

### 5.2 HDF5 File Standards

**Conventions**:
- Dataset names use lowercase with underscores
- Attributes include units, description, creation_date
- Metadata group includes git_version for reproducibility
- All physical quantities in CGS units (Gaussian)

**Required Attributes**:
```python
dataset.attrs['units'] = 'cm^-3'
dataset.attrs['description'] = 'Electron density profile'
dataset.attrs['creation_date'] = '2025-11-04T12:00:00'
```

### 5.3 Binary File Format (EB.dat, amn.dat)

**Type**: Fortran unformatted stream I/O

**Reading in Python**:
```python
import numpy as np

with open('EB.dat', 'rb') as f:
    nr = np.fromfile(f, dtype=np.int32, count=1)[0]
    ntheta = np.fromfile(f, dtype=np.int32, count=1)[0]
    nmodes = np.fromfile(f, dtype=np.int32, count=1)[0]

    for imode in range(nmodes):
        E_r = np.fromfile(f, dtype=np.complex128, count=nr*ntheta)
        E_r = E_r.reshape((nr, ntheta))
        # Continue reading other field components...
```

## 6. Summary of Key Files by Stage

### Preprocessing
- **Input**: GEQDSK, field.dat, experimental profiles (ne, Te, Ti, vt)
- **Output**: INPUT_[shot]_[time].hdf5, amn.dat, equil_r_q_psi.dat

### KiLCA
- **Input**: profiles/*.dat, amn.dat, equil_r_q_psi.dat, *.nml configs
- **Output**: EB.dat, formfactors.dat, resonance.dat

### KIM
- **Input**: profiles/*.dat, KIM_config.nml
- **Output**: kim_out.txt, kim_out.hdf5, kim_kernel.dat

### QL-Balance
- **Input**: INPUT_*.hdf5, flre/EB.dat, vacuum/EB.dat, balance_conf.nml
- **Output**: balance_out.hdf5, profiles_*.dat, diagnostics.dat

### Post-processing
- **Input**: All outputs from main codes
- **Output**: Plots, analysis results, physical quantities

## 7. Reproducibility and Metadata

All HDF5 outputs include:
- Git commit hash of code version
- Timestamp of run
- Input file checksums
- Command-line arguments used
- Compiler version and flags

This ensures complete reproducibility of all computational results.

## 8. Quick Reference: File Extensions

| Extension | Type | Used By |
|-----------|------|---------|
| `.hdf5`, `.h5` | HDF5 binary | All codes (standard interchange format) |
| `.dat` (text) | ASCII data | Profiles, equilibrium, diagnostics |
| `.dat` (binary) | Fortran unformatted | EB.dat, amn.dat (legacy format) |
| `.nml` | Fortran namelist | Configuration files (KIM, QL-Balance, KiLCA) |
| `.inp` | Text input | fouriermodes configuration |
| `.x` | Executable | All Fortran programs |
| `.m` | MATLAB script | Template workflows |
| `.py` | Python script | Interfaces and post-processing |
| `.txt`, `.log` | Text log | Execution logs and status |

## 9. Typical Workflow Command Sequence

```bash
# 1. Build all codes
make all

# 2. Run preprocessing (in run directory)
cd /path/to/run/
matlab -batch "run('/path/to/template_scripts/script_prerun.m')"

# 3. Run QL-Balance for time evolution
cd timeevol/
ln -s ../../INPUT_39711_2000.hdf5 .
mpirun -np 4 ../../build/QL-Balance/ql-balance.x

# 4. Post-process results
python
>>> from postproc_class import balancepost
>>> post = balancepost.BalancePost('out/balance_out.hdf5')
>>> post.plot_profiles_evolution()
```

## 10. Additional Resources

- **Mathematical background**: `/Documentation/` (LaTeX documents)
- **Python tutorials**: `/python/KiLCA_interface/tutorial/`
- **Example notebooks**: `/utility_scripts/python_utility/`
- **MATLAB examples**: `/template_scripts/`

---

**Last Updated**: 2025-11-04
**KAMEL Version**: See git commit hash in repository
