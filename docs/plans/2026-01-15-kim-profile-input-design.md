# KIM Profile Input Module Design

**Date:** 2026-01-15
**Module:** `profile_input_m`
**Location:** `KIM/src/background_equilibrium/profile_input_m.f90`

## Overview

This module integrates the profile preprocessor into KIM to support input profiles in normalized poloidal flux coordinates (sqrt(psi_N)). It automatically detects the input coordinate system, performs coordinate transformation when needed, and calculates the radial electric field from force balance if not provided.

## Requirements

1. Support input profiles in both sqrt(psi_N) and r_eff coordinates
2. Automatic detection of coordinate type with manual override
3. Calculate Er from force balance when Er.dat is missing
4. Use existing btor/r0 from KIM_setup namelist
5. Validate equilibrium parameters against btor_rbig.dat

## Data Flow

```
Input profiles (n_of_psiN.dat, Te_of_psiN.dat, etc.)
              │
              ▼
   ┌─────────────────────────────┐
   │  Coordinate Detection       │
   │  (auto / sqrt_psiN / r_eff) │
   └─────────────────────────────┘
              │
     ┌────────┴────────┐
     │                 │
     ▼                 ▼
 sqrt(psi_N)        r_eff
     │              (max > 2)
     ▼                 │
┌─────────────┐        │
│ profile_    │        │
│ preprocessor│        │
│ (psi→r_eff) │        │
└─────────────┘        │
     │                 │
     └────────┬────────┘
              ▼
   ┌─────────────────────────────┐
   │  Check for Er.dat           │
   │  If missing → force balance │
   │  Write Er_no_Vpol.dat       │
   └─────────────────────────────┘
              │
              ▼
   ┌─────────────────────────────┐
   │  Output: n.dat, Te.dat,     │
   │  Ti.dat, Er.dat, q.dat      │
   │  (in r_eff coordinates)     │
   └─────────────────────────────┘
```

## Namelist Configuration

New namelist group `KIM_profiles` in `KIM_config.nml`:

```fortran
&KIM_profiles
    ! Coordinate type for input profiles
    ! 'auto'      - detect automatically (default)
    ! 'sqrt_psiN' - profiles are in sqrt(psi/psi_max) coordinates
    ! 'r_eff'     - profiles are already in effective radius [cm]
    coord_type = 'auto'

    ! Input directory for raw profiles (n_of_psiN.dat, etc.)
    input_profile_dir = './'

    ! Equilibrium file path (empty = compute from geqdsk)
    equil_file = ''
/
```

**Coordinate detection:**
- Threshold: `max(r_input) > 2.0` → r_eff (typically ~50 cm)
- Threshold: `max(r_input) ≤ 2.0` → sqrt(psi_N) (normalized, max ~1)

**btor/r0 handling:**
- Read from existing `KIM_setup` namelist (`btor`, `r0`)
- Fall back to `btor_rbig.dat` if not set in namelist
- Compare and warn on mismatch (>1%) if both sources exist

## Er Force Balance Calculation

When Er.dat is not found, calculate from force balance with k=0 (no poloidal rotation):

```
E_0r = (T_i / (e_i * n_i)) * ∂n_i/∂r + (1/e_i) * ∂T_i/∂r + (r * B_0 * V_tor) / (c * q * R_0)
```

**Terms:**
1. Density gradient: `(T_i / (e_i * n_i)) * dn_i/dr`
2. Temperature gradient: `(1/e_i) * dT_i/dr`
3. Toroidal rotation: `(r * B0 * Vz) / (c * q * R0)`

**Units:** CGS (statV/cm for Er)
- T_i in eV → convert to erg (multiply by e_charge)

**Output:**
- Write to `Er_no_Vpol.dat`
- Print warning: `"WARNING: Er calculated from force balance without V_pol (k=0)"`

## Module Interface

### Public Procedures

| Procedure | Purpose |
|-----------|---------|
| `prepare_profiles()` | Main entry point, called before `read_profiles()` |

### Internal Procedures

| Procedure | Purpose |
|-----------|---------|
| `read_kim_profiles_namelist()` | Read KIM_profiles namelist |
| `detect_coordinate_type()` | Auto-detect sqrt_psiN vs r_eff |
| `run_preprocessing()` | Call profile_preprocessor for coordinate transform |
| `calculate_Er_from_force_balance()` | Compute Er when missing |
| `validate_btor_rbig()` | Compare namelist vs file values |

### Main Entry Point

```fortran
subroutine prepare_profiles()
    ! Called early in KIM initialization, before read_profiles()

    ! 1. Read KIM_profiles namelist
    ! 2. Detect coordinate type (auto/sqrt_psiN/r_eff)
    ! 3. If sqrt_psiN:
    !    a. Initialize profile_preprocessor with equil_file
    !    b. If equil_file missing → compute equilibrium
    !    c. Process n, Te, Ti, Vz profiles → output to profile_location
    !    d. Check for Er.dat; if missing → calculate from force balance
    ! 4. If r_eff:
    !    a. Verify all required files exist (n.dat, Te.dat, Ti.dat, q.dat)
    !    b. Check for Er.dat; if missing → calculate from force balance
    ! 5. Compare btor/r0 from namelist with btor_rbig.dat, warn on mismatch
end subroutine
```

## Integration with KIM

```
KIM_init.f90
    │
    ├── call prepare_profiles()    ← NEW (before read_profiles)
    │
    ├── call read_profiles()       ← existing (reads n.dat, Te.dat, etc.)
    │
    └── ... rest of initialization
```

## Error Handling

### Coordinate Detection

| Situation | Action |
|-----------|--------|
| `coord_type = 'auto'` and no input profiles found | Stop with error |
| `coord_type = 'sqrt_psiN'` but max(r) > 2 | Warn: "Input appears to be r_eff, but sqrt_psiN specified" |
| `coord_type = 'r_eff'` but max(r) ≤ 2 | Warn: "Input appears to be sqrt_psiN, but r_eff specified" |

### Missing Profiles

| File | Required? | Action if missing |
|------|-----------|-------------------|
| n (or n_of_psiN) | Yes | Stop with error |
| Te (or Te_of_psiN) | Yes | Stop with error |
| Ti (or Ti_of_psiN) | Yes | Stop with error |
| q.dat | Yes | Stop with error (or read from equil_r_q_psi.dat) |
| Vz (or Vz_of_psiN) | No | Warn, set Vz=0 (Er loses rotation term) |
| Er.dat | No | Calculate from force balance, write Er_no_Vpol.dat |

### Equilibrium Validation

| Situation | Action |
|-----------|--------|
| equil_file specified but not found | Warn, fall back to computation |
| btor_rbig.dat not found and btor/r0 not in namelist | Stop with error |
| btor mismatch > 1% between namelist and file | Warn: "btor mismatch: namelist=X, file=Y" |
| r0 mismatch > 1% between namelist and file | Warn: "r0 mismatch: namelist=X, file=Y" |

## Files to Create/Modify

| File | Action |
|------|--------|
| `KIM/src/background_equilibrium/profile_input_m.f90` | Create - new module |
| `KIM/src/general/KIM_init.f90` | Modify - call `prepare_profiles()` |
| `KIM/src/setup/config_mod.f90` | Modify - add KIM_profiles namelist variables |
| `KIM/src/general/read_config.f90` | Modify - read KIM_profiles namelist |
| `KIM/nmls/KIM_config.nml` | Modify - add KIM_profiles namelist group |
| `KIM/CMakeLists.txt` | Modify - add new source file, link kamel_equil |
| `KIM/tests/test_profile_input.f90` | Create - unit tests |

## Testing Strategy

### Unit Tests

1. **Coordinate detection test**
   - Create synthetic profiles with max(r) = 0.95 → expect sqrt_psiN
   - Create synthetic profiles with max(r) = 45.0 → expect r_eff

2. **Er force balance test**
   - Use known analytical profiles (constant n, linear T, zero Vz)
   - Verify Er = (1/e) * dT/dr (pressure gradient only)
   - Compare against expected analytical result

3. **btor/r0 mismatch warning test**
   - Set namelist btor=20000, create btor_rbig.dat with 20500
   - Verify warning is emitted

### Integration Test

- Run KIM with sqrt_psiN input profiles
- Verify preprocessing creates correct r_eff output files
- Verify Er_no_Vpol.dat is created when Er.dat missing

## Dependencies

- `profile_preprocessor_m` from `common/equil/`
- `equil_profiles_m` from `common/equil/` (for equilibrium computation)
- Existing KIM modules: `config_m`, `setup_m`, `species_m`

## Build Integration

Add to `KIM/CMakeLists.txt`:

```cmake
# Link kamel_equil for profile preprocessing
target_link_libraries(KIM_lib PUBLIC kamel_equil)
```
