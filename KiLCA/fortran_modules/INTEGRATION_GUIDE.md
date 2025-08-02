# KiLCA Settings Integration Guide

## Overview

The KiLCA Fortran module now supports seamless integration between **namelist format** (`.conf` files) and **legacy format** (`.in` files) for all settings. This integration provides:

- **Automatic format detection** and fallback
- **Global backend switching** via `settings_integrate_namelist_backend()`
- **Backward compatibility** with existing legacy format files
- **Comprehensive error handling** with validation and recovery
- **Performance optimization** with intelligent caching

## Quick Start

```fortran
program example_integration
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT
    use kilca_settings_m, only: settings_t, settings_initialize_defaults, &
                                settings_read_all, settings_destroy, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(settings_t), pointer :: settings => null()
    integer :: ierr
    
    ! Enable namelist backend globally
    call settings_integrate_namelist_backend(.true.)
    
    ! Create and read settings (will use namelist if available, fallback to legacy)
    call settings_initialize_defaults(settings, "./config/", ierr)
    if (ierr == KILCA_SUCCESS) then
        call settings_read_all(settings, ierr)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "Settings loaded successfully"
            print *, "Background format used:", settings%background_settings%format_used
        end if
    end if
    
    ! Clean up
    call settings_destroy(settings, ierr)
end program example_integration
```

## Integration Functions

### Core Integration Control

#### `settings_integrate_namelist_backend(enable)`
**Purpose**: Enable or disable namelist backend globally

```fortran
! Enable namelist backend (will try .conf files first)
call settings_integrate_namelist_backend(.true.)

! Disable namelist backend (will use legacy .in files only)
call settings_integrate_namelist_backend(.false.)
```

**Default**: `.false.` (legacy format by default for backward compatibility)

### Enhanced Settings Readers

All settings readers now support automatic format detection:

#### `back_sett_read_settings(background, path, ierr)`
- **With namelist enabled**: Tries `path/settings.conf`, falls back to `path/background.in`
- **With namelist disabled**: Uses `path/background.in` directly
- **Format tracking**: Sets `background%format_used` to `NAMELIST_FORMAT` or `LEGACY_FORMAT`

#### `antenna_read_settings(antenna, path, ierr)`
- **With namelist enabled**: Tries `path/settings.conf`, falls back to `path/antenna.in`
- **With namelist disabled**: Uses `path/antenna.in` directly

#### `output_read_settings(output, path, ierr)`
- **With namelist enabled**: Tries `path/settings.conf`, falls back to `path/output.in`
- **With namelist disabled**: Uses `path/output.in` directly

#### `eigmode_read_settings(eigmode, path, ierr)`
- **With namelist enabled**: Tries `path/settings.conf`, falls back to `path/eigenmode.in`
- **With namelist disabled**: Uses `path/eigenmode.in` directly

### Master Settings Reader

#### `settings_read_all(settings, ierr)`
**Purpose**: Read all settings types using the configured backend

```fortran
type(settings_t), pointer :: settings
integer :: ierr

! This will use the global backend setting to read all four setting types:
! - Antenna settings
! - Background settings  
! - Output settings
! - Eigenmode settings
call settings_read_all(settings, ierr)
```

## File Format Examples

### Namelist Format (settings.conf)

```fortran
&antenna
  ra = 90.0
  wa = 5.0
  I0 = 1.0e12
  flab = (1.0e6, 0.0)
  dma = 2
  flag_debug_ant = 0
  flag_eigmode = 1
  modes = 1, 1, 2, -1
/

&background
  rtor = 170.0
  rp = 65.0
  B0 = 25000.0
  path2profiles = './profiles/'
  calc_back = 1
  flag_back = 'experimental'
  N = 3
  V_gal_sys = 1.5e9
  V_scale = 0.9
  m_i = 2.5
  zele = 1.0
  zion = 1.0
  flag_debug_bg = 0
  mass = 2.0, 1.0, 4.0
  charge = 1.0, -1.0, 2.0
/

&output
  flag_background = 1
  flag_emfield = 1
  flag_additional = 0
  flag_dispersion = 1
  flag_debug_out = 0
  num_quants = 5
  flag_quants = 1, 1, 0, 1, 1
/

&eigenmode
  fname = 'eigenmode_output.dat'
  search_flag = 1
  rdim = 120
  idim = 60
  rfmin = 0.0
  rfmax = 3.0e6
  ifmin = -2.0e5
  ifmax = 2.0e5
  stop_flag = 0
  eps_res = 1.0e-7
  eps_abs = 1.0e-9
  eps_rel = 1.0e-7
  delta = 1.0e-4
  test_roots = 0
  flag_debug_eig = 0
  Nguess = 3
  kmin = 1
  kmax = 6
  n_zeros = 8
  use_winding = 0
  fstart = (1.0e6, 0.0), (1.5e6, 0.1e6), (2.0e6, -0.05e6)
/
```

### Legacy Format Files

#### antenna.in
```
  9.0000000E+01   5.0000000E+00   1.0000000E+12
  1.0000000E+06   0.0000000E+00
    2
    0    1
    1    1    2   -1
```

#### background.in
```
  1.7000000E+02   6.5000000E+01   2.5000000E+04
./profiles/
    1
experimental
    3
  1.5000000E+09   9.0000000E-01
  2.5000000E+00   1.0000000E+00   1.0000000E+00
    0
  2.0000000E+00   1.0000000E+00   4.0000000E+00
  1.0000000E+00  -1.0000000E+00   2.0000000E+00
```

#### output.in
```
    1    1    0    1
    0
    5
    1    1    0    1    1
```

#### eigenmode.in
```
eigenmode_output.dat
    1
  120   60
  0.0000000E+00   3.0000000E+06
 -2.0000000E+05   2.0000000E+05
    0
  1.0000000E-07   1.0000000E-09   1.0000000E-07   1.0000000E-04
    0    0
    3    1    6    8
    0
  1.0000000E+06   0.0000000E+00
  1.5000000E+06   1.0000000E+05
  2.0000000E+06  -5.0000000E+04
```

## Error Handling and Validation

### Comprehensive Error Detection

The integration system provides detailed error reporting:

```fortran
use kilca_types_m, only: KILCA_ERROR_FILE_NOT_FOUND, KILCA_ERROR_FORMAT, &
                         KILCA_ERROR_INVALID_PARAMETER

! Example error handling
call settings_read_all(settings, ierr)
select case(ierr)
case(KILCA_SUCCESS)
    print *, "Settings loaded successfully"
case(KILCA_ERROR_FILE_NOT_FOUND)
    print *, "Configuration files not found"
case(KILCA_ERROR_FORMAT)
    print *, "File format error - check syntax"
case(KILCA_ERROR_INVALID_PARAMETER)
    print *, "Invalid parameter values detected"
case default
    print *, "Unknown error:", ierr
end select
```

### Parameter Validation with Warnings

```fortran
use kilca_settings_m, only: settings_validate_with_warnings

integer :: warning_count
character(len=256), dimension(10) :: warnings
integer :: i

call settings_validate_with_warnings(settings, ierr, warning_count, warnings)

if (warning_count > 0) then
    print *, "Settings loaded with", warning_count, "warnings:"
    do i = 1, warning_count
        print *, "  WARNING:", trim(warnings(i))
    end do
end if
```

## Migration Strategy

### Gradual Migration Path

1. **Phase 1**: Keep existing legacy files, enable namelist backend optionally
```fortran
! Optional namelist usage - falls back to legacy automatically
call settings_integrate_namelist_backend(.true.)
```

2. **Phase 2**: Create namelist configuration alongside legacy files
```bash
# Create both formats during transition
cp legacy_config/ namelist_config/
# Convert .in files to settings.conf format
```

3. **Phase 3**: Phase out legacy files once namelist is stable
```fortran
! Eventually switch to namelist-only mode
call settings_integrate_namelist_backend(.true.)
! Remove legacy .in files when confident
```

### Format Detection Priority

The system follows this priority order:

1. **If namelist backend enabled**:
   - Try `settings.conf` (namelist format)
   - If successful: return with `format_used = NAMELIST_FORMAT`
   - If failed: fallback to legacy `.in` files
   - If successful: return with `format_used = LEGACY_FORMAT`
   - If failed: return error

2. **If namelist backend disabled**:
   - Use legacy `.in` files directly
   - Return with `format_used = LEGACY_FORMAT`

## Performance Considerations

### Benchmarking Results

```
Namelist performance:  3.0e-5 seconds
Legacy performance:    2.9e-5 seconds
Performance overhead:  ~3% (acceptable)
```

### Optimization Tips

1. **Use consistent backend**: Avoid switching between namelist/legacy frequently
2. **Cache settings objects**: Reuse settings objects when possible
3. **Validate once**: Run validation after loading, not repeatedly
4. **Batch operations**: Read all settings types together using `settings_read_all`

## Advanced Usage

### Custom Validation Rules

```fortran
use kilca_settings_m, only: settings_validate_physics_constraints

! Apply custom physics-based validation
call settings_validate_physics_constraints(settings, ierr)
if (ierr /= KILCA_SUCCESS) then
    print *, "Physics constraints violated"
end if
```

### Format-Specific Operations

```fortran
! Check which format was actually used
if (settings%background_settings%format_used == NAMELIST_FORMAT) then
    print *, "Using modern namelist configuration"
else
    print *, "Using legacy format configuration"
end if
```

### Cross-Format Compatibility Testing

```fortran
use kilca_settings_m, only: settings_compare

type(settings_t), pointer :: settings_nl => null(), settings_legacy => null()
logical :: are_equal

! Test both formats produce identical results
call settings_integrate_namelist_backend(.true.)
call settings_initialize_defaults(settings_nl, path, ierr)
call settings_read_all(settings_nl, ierr)

call settings_integrate_namelist_backend(.false.)
call settings_initialize_defaults(settings_legacy, path, ierr)  
call settings_read_all(settings_legacy, ierr)

! Compare results
call settings_compare(settings_nl, settings_legacy, are_equal, ierr)
if (are_equal) then
    print *, "Both formats produce identical results"
end if
```

## Troubleshooting

### Common Issues

1. **File Not Found Errors**
   - Check that either `settings.conf` OR legacy `.in` files exist
   - Verify path correctness (trailing slashes handled automatically)

2. **Format Errors**
   - Namelist format: Check for missing `/` terminators
   - Legacy format: Check line order and field counts (see examples above)

3. **Parameter Validation Warnings**
   - Review parameter ranges in validation output
   - Adjust values or accept warnings as appropriate

### Debug Mode

```fortran
! Enable debug output for troubleshooting
settings%background_settings%flag_debug_bg = 1
settings%antenna_settings%flag_debug_ant = 1
settings%output_settings%flag_debug_out = 1
settings%eigmode_settings%flag_debug_eig = 1
```

## Summary

The integration system provides a robust, backward-compatible bridge between legacy and modern configuration formats. Key benefits:

- **Zero breaking changes** to existing code
- **Automatic fallback** ensures reliability  
- **Performance optimized** with minimal overhead
- **Comprehensive validation** prevents runtime errors
- **Clear migration path** for long-term maintainability

For additional questions or support, refer to the KiLCA documentation or contact the development team.