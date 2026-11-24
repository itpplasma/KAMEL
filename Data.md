# KAMEL Data Flow Documentation

This document describes how data flows through the KAMEL framework — what each code reads, what it produces, and how the components connect.

All codes use **CGS/Gaussian units**: lengths in cm, magnetic fields in Gauss, densities in cm⁻³, temperatures in eV, electric fields in statV/cm, velocities in cm/s (or cm/c where noted), currents in statamp.

---

## Complete Workflow

```
1. PREPROCESSING
   ├── Experimental data (GEQDSK, measured profiles)
   ├── fouriermodes.x
   │   └── amn.dat, equil_r_q_psi.dat
   └── Profile preparation (Python/MATLAB)
       └── INPUT.hdf5 (for QL-Balance)

2. PLASMA RESPONSE (choose one)
   ├── KiLCA (cylindrical FLR)
   │   ├── Vacuum run  → EB.dat per mode
   │   └── Plasma run  → EB.dat per mode + formfactors
   │
   └── KIM (integral model, standalone)
       └── fields, kernels, backgrounds (text or HDF5)

3. TRANSPORT (optional)
   └── QL-Balance (wave_code = 'KiLCA' or 'KIM')
       ├── KiLCA: reads KiLCA run directories + INPUT.hdf5
       ├── KIM:   calls KIM in-process with KIM_config.nml + INPUT.hdf5
       └── Writes time-evolved profiles to HDF5

4. POST-PROCESSING
   └── Python (KAMELpy) analysis and visualization
```

---

## 1. KiLCA

### Input

KiLCA uses plain-text configuration files (value before `#` comment) in the run directory:

| File | Purpose |
|------|---------|
| `background.in` | Machine geometry (R₀, a, B₀), profile paths, background options |
| `antenna.in` | RMP coil parameters (radius, current, frequency) |
| `modes.in` | List of (m, n) mode pairs, one per line: `(-6, 2)` |
| `output.in` | Flags controlling which outputs to compute and save |
| `zone_*_<type>.in` | Zone definitions (auto-discovered by glob pattern) |
| `eigmode.in` | Eigenmode solver settings (optional) |

**Profile files** are two-column text (radius, value) in a directory specified by `background.in`. Required profiles: `q`, `n`, `Ti`, `Te`, `Vth`, `Vz`, `Er` (plus `Bth`, `Bz`, `Jth`, `Jz` when reusing precomputed backgrounds).

### Output

KiLCA writes text and Fortran binary files — no HDF5.

```
run_directory/
├── background-data/          # equilibrium quantities
│   ├── b0th.dat, b0z.dat    # B-field components
│   ├── nui.dat, nue.dat     # collision frequencies
│   ├── ni_m.dat, ne_m.dat   # f₀ moments
│   └── ...                  # ~30 diagnostic files
│
├── linear-data/
│   └── m_{m}_n_{n}_flab_[re,im]/
│       ├── EB.dat           # EM fields (r + 6 complex components)
│       ├── mode_data.dat    # mode metadata (m, n, frequencies, r_res)
│       └── zone_*_*.dat     # per-zone diagnostics (power, flux, current)
│
├── dispersion-data/          # if dispersion enabled
│
├── formfactors.*.uff         # Fortran unformatted binary for QL-Balance
│   (modes_dim, grid, mode_index, complex form factor array)
│
└── equil_r_q_psi.dat        # equilibrium r, q, ψ mapping
```

**Key output consumed downstream:**
- `EB.dat` files → QL-Balance reads these per mode
- `formfactors.*.uff` → QL-Balance quasilinear transport coefficients
- `equil_r_q_psi.dat` → QL-Balance island width calculation

---

## 2. KIM

### Input

KIM uses a single Fortran namelist file (default `./KIM_config.nml`). Six namelist groups are required and read in this order; three further groups are optional or conditional. All options are described in the [KIM namelist reference](KIM/nmls/README.md).

| Namelist | Purpose |
|----------|---------|
| `KIM_CONFIG` | Run type, collision model, species count, plasma type |
| `WKB_DISPERSION` | WKB solver settings (mode, solver algorithm, tolerances) |
| `KIM_IO` | I/O paths, HDF5 flags, debug/status levels |
| `KIM_SETUP` | Machine parameters (B₀, R₀), mode numbers, boundary conditions |
| `KIM_GRID` | Radial grid, spline grid, theta integration method and tolerances |
| `KIM_PROFILES` | Profile coordinate type, input/output file paths, GEQDSK path |
| `KIM_PERIODIC` | Optional: forced-periodicity window and Fourier basis (`electrostatic_periodic`) |
| `KIM_FLR2` | Optional: term switches for the `flr2` run type |
| `KIM_species` | Ion charge and mass numbers; read when `read_species_from_namelist = .true.` |

**Run types** (`type_of_run` in `KIM_CONFIG`):

| Value | Description |
|-------|-------------|
| `'electrostatic'` | Poisson equation solve for Φ given external Br |
| `'electrostatic_periodic'` | Electrostatic solve in a forced-periodic radial window (Fourier basis) |
| `'electromagnetic'` | Coupled Poisson-Ampere 2N×2N solve for (Φ, A‖) with self-consistent Br |
| `'WKB_dispersion'` | Radial WKB dispersion relation kr(r) |
| `'flr2'` | Local second-order FLR response imported from KiLCA-FLR2, on KIM's background |
| `'flr2_benchmark'` | FLR2 approximation of KIM's non-local kernel, for validation |

**Profile input** is a two-stage pipeline:

1. **Preprocessing** (`profile_input_m.f90`): auto-detects coordinate type (`coord_type`), transforms √ψ_N → r_eff if needed using equilibrium/GEQDSK, calculates Er from force balance if `Er.dat` is missing, validates units and radial range.

2. **Loading** (`species_mod.f90`): reads the (possibly transformed) two-column `.dat` files from `profile_location`:

| File | Content |
|------|---------|
| `n.dat` | Electron density [cm⁻³] |
| `Te.dat` | Electron temperature [eV] |
| `Ti.dat` | Ion temperature [eV] |
| `q.dat` | Safety factor |
| `Er.dat` | Radial electric field [statV/cm] — optional, calculated if missing |
| `Vz.dat` | Toroidal rotation [cm/s] — optional |

### Output

KIM writes to `output_path/m{m}_n{n}/` in text format, HDF5, or both (controlled by `hdf5_output`).

```
output_path/m{m}_n{n}/
├── fields/
│   ├── Phi          # electrostatic potential
│   ├── Er, Etheta, Ez  # electric field components
│   ├── br           # radial magnetic perturbation
│   ├── jpar, rho    # parallel current, charge density
│   │   (electromagnetic adds:)
│   ├── Apar         # parallel vector potential
│   └── Br_selfconsistent  # self-consistent Br
│
├── kernel/          # sparse kernel matrix elements
├── backs/           # background equilibrium (B₀, pressure, etc.)
├── grid/            # radial grids (cell boundaries and centers)
└── dispersion/      # WKB roots and branches (WKB_dispersion mode)
```

**HDF5 output** (when `hdf5_output = .true.`): a single file (`h5_out_file`) mirrors the directory structure as HDF5 groups, with config/io/setup/grid namelists persisted in top-level groups.

### Data Flow

```
KIM_config.nml          Profile files (*.dat)
(namelist groups)       (possibly in √ψ_N coords)
     │                        │
     └────────┬───────────────┘
              │
     ┌────────▼──────────┐
     │  prepare_profiles │  coordinate transform, Er calc, validation
     └────────┬──────────┘
              │
     ┌────────▼────────┐
     │     KIM.x       │  init → build kernels → solve → postprocess
     └────────┬────────┘
              │
     ┌────────┴──────────────┐
     ▼                       ▼
  Text files              HDF5 file
  (fields/, kernel/,      (all data in one file)
   backs/, grid/)               │
                                ▼
                     Python post-processing
                     (kim.Result, KIMPoissonSolver)
```

---

## 3. QL-Balance

### Input

QL-Balance reads its configuration from `balance_conf.nml` (in the run directory) containing the `BALANCENML` namelist. Key settings include the wave code and its paths, machine parameters, grid control, time stepping, physics options, and HDF5 I/O paths. `TimeEvolutionNTV` runs also read an `&NEORT` group from the same file; `TimeEvolutionStellarator` runs also read `&BALANCE_STELL` from `stell_conf.nml` in the run directory (template: `QL-Balance/namelists/stell_conf.nml`).

**Run types** (`type_of_run` in `BALANCENML`):

| Value | Description |
|-------|-------------|
| `'SingleStep'` | One linear QL step — computes D_ql and \|Br\| at resonance |
| `'TimeEvolution'` | Full nonlinear time evolution with adaptive stepping |
| `'TimeEvolutionNTV'` | Time evolution with neoclassical toroidal viscosity torque from NEO-RT |
| `'ParameterScan'` | Sweep a parameter, doing a SingleStep at each point |
| `'TimeEvolutionStellarator'` | Time evolution for stellarator geometry (1/ν transport) |

**Wave code** (`wave_code` in `BALANCENML`) selects the source of the perturbation fields:

- `'KiLCA'` (template default): fields are read from KiLCA run directories (`flre_path`, `vac_path`).
- `'KIM'`: KIM is called in-process. `kim_config_path` points to a KIM namelist file and `kim_run_type` selects the KIM run type (default `'electromagnetic'`). With `kim_profiles_from_balance = .true.` (default), QL-Balance passes its own profiles to KIM; with `.false.`, KIM reads its own profile files and the modes are read from `modes.in`. With `kim_run_type = 'electrostatic_periodic'`, the modes are given by `kim_n_modes`, `kim_m_list` and `kim_n_list`; this mode requires `kim_profiles_from_balance = .true.` and supports the `SingleStep`, `TimeEvolution` and `ParameterScan` run types.

**Data sources:**

| Source | Content |
|--------|---------|
| `path2inp` (HDF5) | Initial profiles under `/preprocprof/`: r, q, n, Te, Ti, Vth, Vz, Er |
| `flre_path/` (KiLCA dir) | Plasma response — `modes.in`, EB.dat per mode, formfactors (`wave_code = 'KiLCA'`) |
| `vac_path/` (KiLCA dir) | Vacuum response — EB.dat per mode (`wave_code = 'KiLCA'`) |
| `kim_config_path` (namelist) | KIM configuration (`wave_code = 'KIM'`) |
| `equil_path` (text) | Equilibrium mapping: r, q, ψ, ψ_tor columns |

When restarting (`readfromtimestep > 0`), evolved profiles are read from `path2time` at HDF5 group `f_{m}_{n}/fort.1000/{1000+step}/`.

### Output

QL-Balance writes HDF5 to `path2out`. The structure depends on run type.

```
output.hdf5
├── init_params/              # initial profiles (n, Te, Ti, Vz, q, r, Er, Vth)
│
└── f_{m}_{n}/                # per-mode group (or multi_mode)
    ├── r_res                 # resonant surface radius
    ├── stopping_criterion    # why the run stopped
    │
    ├── (SingleStep)
    │   ├── dqle22_res        # QL diffusion at resonance
    │   ├── Br_abs_res        # |Br| at resonance
    │   └── dqle22(r)         # full radial profile
    │
    ├── (TimeEvolution)
    │   ├── br_abs_time(t)    # |Br| at resonance vs time
    │   ├── dqle22_res_time(t)  # D_ql at resonance vs time
    │   ├── T_tot_phi_e(t)    # electron toroidal torque vs time
    │   ├── T_tot_phi_i(t)    # ion toroidal torque vs time
    │   ├── Ipar(t)           # parallel current at resonance (complex)
    │   │
    │   ├── KinProfiles/{step}/   # snapshots every save_prof_time_step
    │   │   ├── rc, n, Te, Ti, Vz, Er, Vth
    │   │   └── sqg_btheta_overc
    │   │
    │   └── LinearProfiles/{step}/  # wave field profiles per linear solve
    │       ├── Br_abs, Br_Re, Br_Im
    │       ├── dqle22
    │       └── Jpe_abs, Jpi_abs
    │
    └── (ParameterScan)
        └── {scan_point}/
            ├── dqle22_res
            └── br_abs_res
```

### Data Flow

```
INPUT.hdf5           KiLCA FLRE run       KiLCA vacuum run    balance_conf.nml
(/preprocprof/)      (EB.dat, formfacs)   (EB.dat)            (BALANCENML)
                     ── or, with wave_code = 'KIM': KIM_config.nml ──
      │                     │                   │                    │
      └─────────────────────┴───────────────────┴────────────────────┘
                                    │
                            ┌───────▼───────┐
                            │ ql-balance.x  │
                            └───────┬───────┘
                                    │
                             Output HDF5 file
                          (profiles, time series,
                           transport coefficients)
                                    │
                                    ▼
                         Python post-processing
                         (QL_Balance_interface)
```

---

## File Format Summary

| Format | Where used | Notes |
|--------|-----------|-------|
| `.in` text files | KiLCA config | Value before `#` comment, line-oriented |
| `.nml` namelists | KIM, QL-Balance config | Standard Fortran namelists |
| `.dat` two-column text | All profile I/O | `radius  value`, one pair per line |
| `.uff` Fortran binary | KiLCA formfactors | Unformatted sequential, read by QL-Balance |
| `.hdf5` | KIM output, QL-Balance I/O | Structured hierarchical data |
