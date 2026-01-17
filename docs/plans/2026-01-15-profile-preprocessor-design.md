# Profile Preprocessor Module Design

**Date:** 2026-01-15
**Module:** `profile_preprocessor_m`
**Location:** `common/equil/profile_preprocessor.f90`

## Overview

The profile preprocessor module maps plasma profiles from normalized poloidal flux coordinates to effective radius coordinates. It reads input profiles (density, temperatures, rotation) as functions of psi_N or sqrt(psi_N), uses equilibrium data to establish the psi_N → r_eff mapping, and outputs profiles in the format expected by KiLCA/KIM/QL-Balance.

## Requirements

1. Process EQDSK equilibrium to obtain psi → r_eff mapping
2. Normalize poloidal flux by its maximum value (from equil_r_q_psi.dat)
3. Read input profiles: n_of_psiN.dat, Te_of_psiN.dat, Ti_of_psiN.dat, Vz_of_psiN.dat
4. Support both sqrt(psi_N) and psi_N as input coordinate (default: sqrt(psi_N))
5. Map profiles to r_eff using cubic spline interpolation
6. Write output files: n.dat, Te.dat, Ti.dat, Vz.dat
7. Provide programmatic access via getter functions

## Module Structure

### Public Constants

```fortran
integer, parameter :: COORD_SQRT_PSIN = 1  ! sqrt(psi/psi_max) - default
integer, parameter :: COORD_PSIN = 2       ! psi/psi_max
```

### Derived Type

```fortran
type, public :: profile_preprocessor_t
    ! Configuration
    character(len=512) :: equil_file = ''       ! Path to equil_r_q_psi.dat
    character(len=512) :: output_dir = '.'      ! Where to write output
    character(len=512) :: input_dir = '.'       ! Where to read input profiles
    integer :: coord_type = COORD_SQRT_PSIN     ! Input coordinate convention

    ! Input filenames (configurable)
    character(len=128) :: n_input_file = 'n_of_psiN.dat'
    character(len=128) :: Te_input_file = 'Te_of_psiN.dat'
    character(len=128) :: Ti_input_file = 'Ti_of_psiN.dat'
    character(len=128) :: Vz_input_file = 'Vz_of_psiN.dat'

    ! Equilibrium mapping arrays
    integer :: nrad = 0
    double precision, allocatable :: r_eff(:)    ! Effective radius grid
    double precision, allocatable :: psi_n(:)    ! Normalized poloidal flux
    double precision :: psi_max = 0.d0           ! Normalization value

    ! Output profiles on r_eff grid
    double precision, allocatable :: n(:), Te(:), Ti(:), Vz(:)

    logical :: initialized = .false.
contains
    procedure :: init, process, write_output, cleanup
    procedure :: get_n, get_Te, get_Ti, get_Vz
end type
```

### Type-Bound Procedures

| Procedure | Purpose |
|-----------|---------|
| `init(namelist_file, [equil_in])` | Read config, load/compute equilibrium mapping |
| `process()` | Read input profiles, interpolate to r_eff grid |
| `write_output()` | Write n.dat, Te.dat, Ti.dat, Vz.dat |
| `cleanup()` | Deallocate arrays |
| `get_n(r)`, `get_Te(r)`, `get_Ti(r)`, `get_Vz(r)` | Interpolate profile at arbitrary r_eff |

## Namelist Configuration

**Namelist group:** `&profile_preprocessor`

```fortran
namelist /profile_preprocessor/ &
    equil_file,        &  ! Path to equil_r_q_psi.dat (empty = recompute)
    output_dir,        &  ! Output directory for n.dat, Te.dat, etc.
    input_dir,         &  ! Input directory for profiles
    coord_type,        &  ! 1 = sqrt(psi_N), 2 = psi_N
    n_input_file,      &  ! Input filename for density
    Te_input_file,     &  ! Input filename for Te
    Ti_input_file,     &  ! Input filename for Ti
    Vz_input_file,     &  ! Input filename for Vz
    nsqpsi,            &  ! Grid points for equilibrium recomputation
    nstep,             &  ! Integration steps (if recomputing)
    nsurfmax,          &  ! Surface search points (if recomputing)
    niter                 ! Newton iterations (if recomputing)
```

**Example namelist file:**

```
&profile_preprocessor
    equil_file = '/path/to/equil_r_q_psi.dat'
    input_dir = '/path/to/input_profiles/'
    output_dir = '/path/to/output/'
    coord_type = 1
    n_input_file = 'ne_profile.dat'
/
```

## Data Flow

### Initialization (`init`)

1. Read namelist file
2. Obtain equilibrium mapping (one of three paths):
   - If `equil_in` argument provided → use it directly
   - Else if `equil_file` is set and exists → read equil_r_q_psi.dat
   - Else → call `init_equilibrium_field()`, create `equil_profiles_t`, compute profiles, write to `output_dir/equil_r_q_psi.dat`
3. Extract r_eff and psi arrays, compute `psi_n = psi / psi_max`
4. Allocate output profile arrays

### Processing (`process`)

For each profile (n, Te, Ti, Vz):

1. Construct full input path: `input_dir / filename`
2. If file doesn't exist → skip with warning, leave array unallocated
3. Read file, skipping header lines (starting with `#` or non-numeric)
4. Parse two columns: (coord_in, value)
5. Convert coord_in to psi_n based on coord_type:
   - `COORD_SQRT_PSIN`: psi_n_in = coord_in²
   - `COORD_PSIN`: psi_n_in = coord_in
6. Check bounds: if input extends beyond equilibrium → warn, truncate
7. Build cubic spline using libneo's spline routines
8. Interpolate onto equilibrium psi_n grid
9. Store in output array

### Output (`write_output`)

For each allocated profile array:

1. Open `output_dir / {n,Te,Ti,Vz}.dat`
2. Write two columns: `r_eff(i), profile(i)`
3. Close file

## Input File Format

Flexible two-column ASCII format:
- Lines starting with `#` or non-numeric characters are skipped (header)
- Data lines: `<coordinate> <value>` (whitespace separated)
- Coordinate must be monotonically increasing

Example:
```
# Electron density profile
# sqrt(psi_N)  n_e [cm^-3]
0.0  1.5e13
0.1  1.48e13
0.2  1.42e13
...
```

## Output File Format

Two-column ASCII, no header:
```
<r_eff>  <value>
```

Example `n.dat`:
```
0.0123  1.5e13
0.0456  1.48e13
...
```

## Error Handling

### File/IO Errors

| Situation | Action |
|-----------|--------|
| Namelist file not found | Stop with error |
| Namelist read error | Stop with error |
| equil_file specified but not found | Fall back to recomputation, warn |
| Input profile file not found | Skip with warning |
| Output directory doesn't exist | Create it; stop if creation fails |
| Write error | Stop with error |

### Data Validation

| Situation | Action |
|-----------|--------|
| Input profile has < 2 points | Stop with error |
| Input coordinate not monotonic | Stop with error |
| Input extends beyond equilibrium | Warn and truncate |
| Input doesn't cover equilibrium start | Warn; stop if gap > 10% |
| Negative density values | Warn (don't stop) |
| Negative temperature values | Warn (don't stop) |
| NaN or Inf in input | Stop with error |

### Runtime Errors

| Situation | Action |
|-----------|--------|
| `get_*` before `process` | Stop with error |
| `get_*` with r outside range | Warn, clamp to boundary |
| `write_output` with no profiles | Warn, write nothing |

Error message format:
```
[profile_preprocessor_m:procedure] ERROR/WARNING: message
```

## Dependencies

- `equil_profiles_m` - equilibrium computation fallback
- libneo spline routines - cubic interpolation
- Standard Fortran file I/O

## Build Integration

Add to `common/equil/CMakeLists.txt`:

```cmake
set(EQUIL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/rk4_integrator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mag_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/field_line_rhs.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/equil_profiles.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/profile_preprocessor.f90
)
```

## Usage Example

```fortran
use profile_preprocessor_m

type(profile_preprocessor_t) :: pp

! Initialize and process
call pp%init('profile_preprocessor.nml')
call pp%process()
call pp%write_output()

! Or use directly without files:
Te_at_r = pp%get_Te(0.25d0)

call pp%cleanup()
```
