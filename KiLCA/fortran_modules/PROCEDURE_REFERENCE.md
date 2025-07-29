# KiLCA Settings Module - Procedure Reference

This document provides detailed reference documentation for all procedures in the `kilca_settings_m` module.

## Table of Contents

1. [Settings Lifecycle Management](#settings-lifecycle-management)
2. [Antenna Procedures](#antenna-procedures)
3. [Background Procedures](#background-procedures)
4. [Output Procedures](#output-procedures)
5. [Eigenmode Procedures](#eigenmode-procedures)
6. [Validation Procedures](#validation-procedures)
7. [Initialization Procedures](#initialization-procedures)
8. [Input/Output Procedures](#inputoutput-procedures)
9. [Error Handling](#error-handling)
10. [Error Codes](#error-codes)

## Settings Lifecycle Management

### settings_create(sd, path, ierr)
**Purpose**: Create and initialize a new settings structure.

**Parameters**:
- `sd` (type(settings_t), pointer, intent(out)): Pointer to settings structure
- `path` (character(*), intent(in)): Project path for settings
- `ierr` (integer, intent(out)): Error status

**Error Codes**:
- `KILCA_SUCCESS`: Success
- `KILCA_ERROR_MEMORY`: Memory allocation failed
- `KILCA_ERROR_INVALID_INPUT`: Empty path provided

**Usage**:
```fortran
type(settings_t), pointer :: sd
integer :: ierr
call settings_create(sd, "/path/to/project", ierr)
```

### settings_destroy(sd, ierr)
**Purpose**: Destroy settings structure and free memory.

**Parameters**:
- `sd` (type(settings_t), pointer, intent(inout)): Settings structure to destroy
- `ierr` (integer, intent(out)): Error status

## Antenna Procedures

### antenna_initialize_defaults(ant, ierr)
**Purpose**: Initialize antenna settings with default values.

**Parameters**:
- `ant` (type(antenna_t), intent(out)): Antenna structure
- `ierr` (integer, intent(out)): Error status

**Default Values**:
- `ra`: 0.0 cm (antenna radius)
- `wa`: 0.0 cm (current density layer width)
- `I0`: 0.0 statamps (current in antenna coils)
- `flab`: (0.0, 0.0) Hz (laboratory frequency)
- `dma`: 0 (dimension of modes array)
- `flag_debug`: 0 (debug flag off)
- `flag_eigmode`: 0 (eigenmode search off)

### antenna_initialize_custom(ant, ra, wa, I0, ierr)
**Purpose**: Initialize antenna with custom parameters.

**Parameters**:
- `ant` (type(antenna_t), intent(out)): Antenna structure
- `ra` (real(dp), optional): Antenna radius
- `wa` (real(dp), optional): Current layer width
- `I0` (real(dp), optional): Antenna current
- `ierr` (integer, intent(out)): Error status

### antenna_validate(ant, is_valid, error_msg, ierr)
**Purpose**: Validate antenna settings for physical consistency.

**Parameters**:
- `ant` (type(antenna_t), intent(in)): Antenna structure to validate
- `is_valid` (logical, intent(out)): Validation result
- `error_msg` (character(*), intent(out)): Error message if invalid
- `ierr` (integer, intent(out)): Error status

**Validation Rules**:
- `ra` must be non-negative
- `wa` must be non-negative
- `I0` must be non-negative
- If `dma` > 0, `modes` array must be allocated with size 2*dma

### antenna_print_settings(ant, ierr)
**Purpose**: Print antenna settings to standard output.

### antenna_print_settings_to_unit(ant, unit, ierr)
**Purpose**: Print antenna settings to specified unit.

### antenna_deep_copy(source, target, ierr)
**Purpose**: Create a deep copy of antenna settings.

### antenna_compare(ant1, ant2, is_equal, ierr)
**Purpose**: Compare two antenna structures for equality.

## Background Procedures

### back_sett_initialize_defaults(bs, ierr)
**Purpose**: Initialize background settings with meaningful defaults.

**Default Values**:
- `rtor`: 625.0 cm (ASDEX Upgrade major radius)
- `rp`: 200.0 cm (typical plasma minor radius)
- `B0`: 20000.0 G (typical toroidal field)
- `calc_back`: 1 (calculate background profiles)
- `N`: 5 (spline degree, must be odd)
- `V_gal_sys`: 0.0 cm/s (velocity of moving frame)
- `V_scale`: 1.0 (scale factor for Vz velocity profile)
- `m_i`: 1.0 (ion mass in proton masses)
- `zele`: 1.0 (electron collision coefficient)
- `zion`: 1.0 (ion collision coefficient)
- `flag_debug`: 0 (debug off)
- `huge_factor`: 1.0e30 (numerical cutoff)
- `flag_back`: "normal" (background calculation mode)

### back_sett_validate(bs, is_valid, error_msg, ierr)
**Purpose**: Validate background settings.

**Validation Rules**:
- `rtor` must be positive
- `rp` must be positive and less than `rtor`
- `B0` must be positive
- `N` must be odd and positive
- Physical parameters must be positive

## Output Procedures

### output_sett_initialize_defaults(os, ierr)
**Purpose**: Initialize output settings with defaults.

**Default Values**:
- `flag_background`: 1 (compute background data)
- `flag_emfield`: 1 (compute electromagnetic field data)
- `flag_additional`: 0 (don't compute additional quantities)
- `flag_dispersion`: 0 (don't compute dispersion)
- `num_quants`: 0 (no additional quantities)
- `flag_debug`: 0 (debug off)

### output_sett_validate(os, is_valid, error_msg, ierr)
**Purpose**: Validate output settings.

## Eigenmode Procedures

### eigmode_sett_initialize_defaults(es, ierr)
**Purpose**: Initialize eigenmode settings with defaults.

**Default Values**:
- `search_flag`: 0 (no eigenmode search)
- `rdim`: 100 (real frequency mesh dimension)
- `rfmin`: 0.0 Hz (minimum real frequency)
- `rfmax`: 1.0e9 Hz (maximum real frequency)
- `idim`: 100 (imaginary frequency mesh dimension)
- `ifmin`: -1.0e6 Hz (minimum imaginary frequency)
- `ifmax`: 1.0e6 Hz (maximum imaginary frequency)
- `stop_flag`: 0 (stopping criteria)
- `eps_res`: 1.0e-6 (residual error parameter)
- `eps_abs`: 1.0e-8 (absolute error parameter)
- `eps_rel`: 1.0e-6 (relative error parameter)
- `delta`: 1.0e-6 (delta for derivative)
- `test_roots`: 0 (don't test roots)
- `flag_debug`: 0 (debug off)
- `Nguess`: 0 (no initial guesses)
- `kmin`: 1 (minimum k value)
- `kmax`: 10 (maximum k value)
- `n_zeros`: 10 (number of zeros to find)
- `use_winding`: 0 (don't use winding number)
- `fname`: "eigenmode_output.dat" (output filename)

### eigmode_sett_validate(es, is_valid, error_msg, ierr)
**Purpose**: Validate eigenmode settings.

## Validation Procedures

### settings_validate(sd, is_valid, error_msg, ierr)
**Purpose**: Validate all settings in the structure.

### settings_validate_complete(sd, is_valid, error_msg, ierr)
**Purpose**: Complete validation including pointer consistency.

### settings_validate_consistency(sd, is_valid, error_msg, ierr)
**Purpose**: Validate consistency between different settings components.

### settings_validate_with_context(sd, context, error_msg, ierr)
**Purpose**: Validate settings with context information for better error messages.

## Initialization Procedures

### settings_initialize_defaults(sd, path, ierr)
**Purpose**: Initialize complete settings structure with defaults.

**Parameters**:
- `sd` (type(settings_t), pointer, intent(out)): Settings structure
- `path` (character(*), intent(in)): Project path
- `ierr` (integer, intent(out)): Error status

This procedure:
1. Creates the settings structure
2. Initializes all sub-components with defaults
3. Sets up internal pointers
4. Validates the result

## Input/Output Procedures

### settings_print_all(sd, ierr)
**Purpose**: Print all settings to standard output.

### antenna_print_settings(ant, ierr)
**Purpose**: Print antenna settings with formatted output.

**Output Format**:
```
antenna radius: <value> cm
antenna current layer width: <value> cm
antenna coils current: <value> statamps
antenna lab frequency: (<real>, <imag>) 1/s
dimension of modes array: <value>
flag for debugging mode: <value>
flag for eigmode search: <value>
array of mode numbers (m, n): <pairs>
```

### back_sett_print_settings(bs, ierr)
**Purpose**: Print background settings.

### output_sett_print_settings(os, ierr)
**Purpose**: Print output settings.

### eigmode_sett_print_settings(es, ierr)
**Purpose**: Print eigenmode settings.

All print procedures support both standard output and file output variants.

## Error Handling

### settings_format_error_message(error_code, context, parameter, error_msg)
**Purpose**: Format detailed error messages with context.

### settings_get_error_name(error_code, error_name, ierr)
**Purpose**: Get human-readable error name from error code.

### settings_attempt_recovery(ant, recovered, ierr)
**Purpose**: Attempt to recover from invalid antenna state.

**Recovery Actions**:
- Fix negative `ra` values
- Allocate `modes` array if `dma` > 0
- Fix array size inconsistencies
- Deallocate arrays when `dma` = 0

### settings_get_detailed_validation_errors(bs, detailed_error, ierr)
**Purpose**: Get detailed validation error information.

### settings_log_error(error_code, context, parameter, log_file, ierr)
**Purpose**: Log errors to file.

### settings_check_error_log_exists(log_file, logged, ierr)
**Purpose**: Check if error log file exists.

### settings_clear_error_log(log_file, ierr)
**Purpose**: Clear error log file.

## Error Codes

### Success Codes
- `KILCA_SUCCESS` (0): Operation completed successfully

### Error Codes
- `KILCA_ERROR` (-1): Generic error
- `KILCA_ERROR_MEMORY` (-2): Memory allocation failed
- `KILCA_ERROR_FILE` (-3): File I/O error
- `KILCA_ERROR_INVALID_INPUT` (-4): Invalid input parameter
- `KILCA_ERROR_CONVERGENCE` (-5): Numerical convergence failure
- `KILCA_ERROR_BOUNDS` (-6): Array bounds error
- `KILCA_ERROR_NOT_IMPLEMENTED` (-7): Feature not implemented
- `KILCA_ERROR_FORMAT` (-8): File format or parsing error

### Error Checking Pattern
Always check error codes after procedure calls:

```fortran
integer :: ierr
call some_procedure(args, ierr)
if (ierr /= KILCA_SUCCESS) then
    ! Handle error
    return
end if
```

### Error Recovery
Use error recovery procedures when possible:

```fortran
logical :: recovered
call settings_attempt_recovery(ant, recovered, ierr)
if (recovered) then
    ! Continue with corrected data
end if
```

## Threading and Memory Safety

All procedures in this module are thread-safe for read operations. Write operations require external synchronization. Memory management follows RAII principles - always call `settings_destroy` to free allocated memory.

## See Also

- [SETTINGS_USAGE.md](SETTINGS_USAGE.md) - Usage examples and tutorials
- [examples/](examples/) - Complete working examples