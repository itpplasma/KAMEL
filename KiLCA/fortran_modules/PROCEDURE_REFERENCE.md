# KiLCA Settings Module Procedure Reference

## Table of Contents

1. [Settings Lifecycle](#settings-lifecycle)
2. [Antenna Procedures](#antenna-procedures)
3. [Background Procedures](#background-procedures)
4. [Output Procedures](#output-procedures)
5. [Eigenmode Procedures](#eigenmode-procedures)
6. [Validation Procedures](#validation-procedures)
7. [I/O Procedures](#io-procedures)
8. [Error Handling Procedures](#error-handling-procedures)
9. [Utility Procedures](#utility-procedures)

---

## Settings Lifecycle

### `settings_create(sd, path, ierr)`
**Purpose**: Create and allocate settings structure  
**Parameters**:
- `sd` (out): `type(settings_t), pointer` - Settings structure pointer
- `path` (in): `character(len=*)` - Project path (must not be empty)
- `ierr` (out): `integer` - Error status

**Returns**: `KILCA_SUCCESS` or `KILCA_ERROR_INVALID_INPUT`/`KILCA_ERROR_MEMORY`

**Example**:
```fortran
type(settings_t), pointer :: sd
call settings_create(sd, "/path/to/project/", ierr)
```

### `settings_destroy(sd, ierr)`
**Purpose**: Destroy settings structure and free all memory  
**Parameters**:
- `sd` (inout): `type(settings_t), pointer` - Settings structure pointer
- `ierr` (out): `integer` - Error status

**Returns**: `KILCA_SUCCESS` or `KILCA_ERROR_INVALID_INPUT`

### `settings_initialize_defaults(sd, path, ierr)`
**Purpose**: Create settings and initialize all components with defaults  
**Parameters**:
- `sd` (out): `type(settings_t), pointer` - Settings structure pointer
- `path` (in): `character(len=*)` - Project path
- `ierr` (out): `integer` - Error status

---

## Antenna Procedures

### `antenna_initialize_defaults(ant, ierr)`
**Purpose**: Initialize antenna with default values  
**Parameters**:
- `ant` (inout): `type(antenna_t)` - Antenna structure
- `ierr` (out): `integer` - Error status

### `antenna_initialize_custom(ant, ra, wa, I0, ierr)`
**Purpose**: Initialize antenna with custom parameters  
**Parameters**:
- `ant` (inout): `type(antenna_t)` - Antenna structure
- `ra` (in, optional): `real(dp)` - Antenna radius (cm)
- `wa` (in, optional): `real(dp)` - Current layer width
- `I0` (in, optional): `real(dp)` - Current (statamp)
- `ierr` (out): `integer` - Error status

### `antenna_set_parameters(ant, ra, wa, I0, ierr)`
**Purpose**: Set antenna parameters  
**Parameters**:
- `ant` (inout): `type(antenna_t)` - Antenna structure
- `ra` (in, optional): `real(dp)` - Antenna radius (cm)
- `wa` (in, optional): `real(dp)` - Current layer width  
- `I0` (in, optional): `real(dp)` - Current (statamp)
- `ierr` (out): `integer` - Error status

### `antenna_validate(ant, is_valid, error_msg, ierr)`
**Purpose**: Validate antenna settings  
**Parameters**:
- `ant` (in): `type(antenna_t)` - Antenna structure
- `is_valid` (out): `logical` - Validation result
- `error_msg` (out): `character(len=*)` - Error message if invalid
- `ierr` (out): `integer` - Error status

**Validation Rules**:
- `ra > 0` (antenna radius must be positive)
- `wa >= 0` (current layer width must be non-negative)
- `I0 >= 0` (current must be non-negative)
- `dma >= 0` (mode dimension must be non-negative)
- If `dma > 0`, `modes` array must be allocated with size `dma`
- Mode values must be non-negative integers

### `antenna_deep_copy(src, dst, ierr)`
**Purpose**: Create deep copy of antenna settings  
**Parameters**:
- `src` (in): `type(antenna_t)` - Source antenna
- `dst` (inout): `type(antenna_t)` - Destination antenna
- `ierr` (out): `integer` - Error status

### `antenna_compare(ant1, ant2, is_equal, ierr)`
**Purpose**: Compare two antenna structures for equality  
**Parameters**:
- `ant1`, `ant2` (in): `type(antenna_t)` - Antenna structures to compare
- `is_equal` (out): `logical` - Comparison result
- `ierr` (out): `integer` - Error status

### `antenna_print_settings(ant, ierr)`
**Purpose**: Print antenna settings to stdout  
**Parameters**:
- `ant` (in): `type(antenna_t)` - Antenna structure
- `ierr` (out): `integer` - Error status

### `antenna_print_settings_to_unit(ant, unit, ierr)`
**Purpose**: Print antenna settings to specific file unit  
**Parameters**:
- `ant` (in): `type(antenna_t)` - Antenna structure
- `unit` (in): `integer` - File unit number
- `ierr` (out): `integer` - Error status

### `antenna_read_settings(ant, path, ierr)`
**Purpose**: Read antenna settings from file (stub implementation)  
**Parameters**:
- `ant` (inout): `type(antenna_t)` - Antenna structure
- `path` (in): `character(len=*)` - Project path
- `ierr` (out): `integer` - Error status

---

## Background Procedures

### `back_sett_initialize_defaults(bs, ierr)`
**Purpose**: Initialize background settings with ASDEX Upgrade defaults  
**Parameters**:
- `bs` (inout): `type(back_sett_t)` - Background settings structure
- `ierr` (out): `integer` - Error status

**Default Values**:
- `rtor = 625.0` cm (ASDEX Upgrade major radius)
- `rp = 200.0` cm (typical plasma radius)
- `B0 = 20000.0` G (2 Tesla toroidal field)
- `calc_back = 1`
- `N = 5` (spline order)
- `V_gal_sys = 0.0`
- `V_scale = 1.0`
- `m_i = 1.0` (proton mass)
- `zele = zion = 1.0`
- `flag_back = "normal"`
- `huge_factor = 1.0e30`

### `back_sett_validate(bs, is_valid, error_msg, ierr)`
**Purpose**: Validate background settings  
**Parameters**:
- `bs` (in): `type(back_sett_t)` - Background settings structure
- `is_valid` (out): `logical` - Validation result
- `error_msg` (out): `character(len=*)` - Error message if invalid
- `ierr` (out): `integer` - Error status

**Validation Rules**:
- `rtor > 0` (major radius must be positive)
- `rp > 0` (plasma radius must be positive)
- `rp < rtor` (plasma radius must be less than major radius)
- `B0 > 0` (magnetic field must be positive)
- `N` must be odd and >= 3 (spline order constraint)
- `m_i > 0` (ion mass must be positive)
- `zele, zion > 0` (collision coefficients must be positive)
- `huge_factor > 0` (huge factor must be positive)

### `back_sett_set_calc_flag(bs, calc_flag, ierr)`
**Purpose**: Set background calculation flag  
**Parameters**:
- `bs` (inout): `type(back_sett_t)` - Background settings structure
- `calc_flag` (in): `integer` - Calculation flag value
- `ierr` (out): `integer` - Error status

### `back_sett_deep_copy(src, dst, ierr)`
**Purpose**: Create deep copy of background settings  
**Parameters**:
- `src` (in): `type(back_sett_t)` - Source background settings
- `dst` (inout): `type(back_sett_t)` - Destination background settings
- `ierr` (out): `integer` - Error status

### `back_sett_compare(bs1, bs2, is_equal, ierr)`
**Purpose**: Compare two background settings for equality  
**Parameters**:
- `bs1`, `bs2` (in): `type(back_sett_t)` - Background settings to compare
- `is_equal` (out): `logical` - Comparison result
- `ierr` (out): `integer` - Error status

---

## Output Procedures

### `output_sett_initialize_defaults(os, ierr)`
**Purpose**: Initialize output settings with defaults  
**Parameters**:
- `os` (inout): `type(output_sett_t)` - Output settings structure
- `ierr` (out): `integer` - Error status

**Default Values**:
- `flag_background = 1` (compute background)
- `flag_emfield = 1` (compute EM fields)
- `flag_additional = 0` (skip additional quantities)
- `flag_dispersion = 0` (skip dispersion)
- `num_quants = 0`
- `flag_debug = 0`

### `output_sett_validate(os, is_valid, error_msg, ierr)`
**Purpose**: Validate output settings  
**Parameters**:
- `os` (in): `type(output_sett_t)` - Output settings structure
- `is_valid` (out): `logical` - Validation result
- `error_msg` (out): `character(len=*)` - Error message if invalid
- `ierr` (out): `integer` - Error status

**Validation Rules**:
- All flags must be 0, 1, or 2
- `num_quants >= 0`
- If `num_quants > 0`, `flag_quants` must be allocated with correct size
- All elements in `flag_quants` must be 0 or 1

### `output_sett_set_flags(os, flag_background, flag_emfield, ierr)`
**Purpose**: Set output flags  
**Parameters**:
- `os` (inout): `type(output_sett_t)` - Output settings structure
- `flag_background` (in, optional): `integer` - Background flag
- `flag_emfield` (in, optional): `integer` - EM field flag
- `ierr` (out): `integer` - Error status

---

## Eigenmode Procedures

### `eigmode_sett_initialize_defaults(es, ierr)`
**Purpose**: Initialize eigenmode settings with defaults  
**Parameters**:
- `es` (inout): `type(eigmode_sett_t)` - Eigenmode settings structure
- `ierr` (out): `integer` - Error status

**Default Values**:
- `search_flag = 0`
- `rdim = idim = 100` (frequency grid dimensions)
- `rfmin = 0.0, rfmax = 1.0e9` (real frequency range)
- `ifmin = -1.0e6, ifmax = 1.0e6` (imaginary frequency range)
- `eps_res = 1.0e-6` (residual tolerance)
- `eps_abs = 1.0e-8` (absolute tolerance)
- `eps_rel = 1.0e-6` (relative tolerance)
- `delta = 1.0e-6`
- `kmin = 1, kmax = 10`
- `n_zeros = 10`
- `fname = "eigenmode_output.dat"`

### `eigmode_sett_validate(es, is_valid, error_msg, ierr)`
**Purpose**: Validate eigenmode settings  
**Parameters**:
- `es` (in): `type(eigmode_sett_t)` - Eigenmode settings structure
- `is_valid` (out): `logical` - Validation result
- `error_msg` (out): `character(len=*)` - Error message if invalid
- `ierr` (out): `integer` - Error status

**Validation Rules**:
- `rdim, idim > 0` (grid dimensions must be positive)
- `rfmax > rfmin` (real frequency range must be valid)
- `ifmax > ifmin` (imaginary frequency range must be valid)
- All tolerance values must be positive
- `kmin <= kmax` and both must be positive
- `n_zeros > 0`
- If `Nguess > 0`, `fstart` must be allocated with size `Nguess`

### `eigmode_sett_set_search_flag(es, search_flag, ierr)`
**Purpose**: Set eigenmode search flag  
**Parameters**:
- `es` (inout): `type(eigmode_sett_t)` - Eigenmode settings structure
- `search_flag` (in): `integer` - Search flag value
- `ierr` (out): `integer` - Error status

---

## Validation Procedures

### `settings_validate(sd, ierr)`
**Purpose**: Validate all settings components  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
- `ierr` (out): `integer` - Error status

### `settings_validate_complete(sd, ierr)`
**Purpose**: Complete validation including consistency checks  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
- `ierr` (out): `integer` - Error status

### `settings_validate_consistency(sd, ierr)`
**Purpose**: Check inter-component consistency  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
- `ierr` (out): `integer` - Error status

**Consistency Checks**:
- Antenna frequency compatible with eigenmode search range
- Background spline order compatible with eigenmode grid resolution
- Output flags consistent with computation requirements

### `settings_validate_with_context(sd, context, error_msg, ierr)`
**Purpose**: Validate settings with contextual error reporting  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
- `context` (in): `character(len=*)` - Context description
- `error_msg` (out): `character(len=*)` - Contextual error message
- `ierr` (out): `integer` - Error status

---

## I/O Procedures

### `settings_read_all(sd, ierr)`
**Purpose**: Read all settings from configuration files  
**Parameters**:
- `sd` (inout): `type(settings_t)` - Settings structure
- `ierr` (out): `integer` - Error status

### `settings_print_all(sd, ierr)`
**Purpose**: Print all settings to stdout  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
- `ierr` (out): `integer` - Error status

---

## Error Handling Procedures

### `settings_format_error_message(error_code, context, parameter, error_msg)`
**Purpose**: Format detailed error message with context  
**Parameters**:
- `error_code` (in): `integer` - Error code
- `context` (in): `character(len=*)` - Error context
- `parameter` (in): `character(len=*)` - Parameter name
- `error_msg` (out): `character(len=*)` - Formatted error message

### `settings_get_error_name(error_code, error_name, ierr)`
**Purpose**: Get human-readable error code name  
**Parameters**:
- `error_code` (in): `integer` - Error code
- `error_name` (out): `character(len=*)` - Error name string
- `ierr` (out): `integer` - Error status

### `settings_attempt_recovery(ant, recovered, ierr)`
**Purpose**: Attempt to recover from invalid antenna settings  
**Parameters**:
- `ant` (inout): `type(antenna_t)` - Antenna structure to recover
- `recovered` (out): `logical` - Recovery success flag
- `ierr` (out): `integer` - Error status

**Recovery Actions**:
- Set negative values to zero
- Allocate missing arrays if needed
- Apply reasonable default values

### `settings_get_detailed_validation_errors(bs, detailed_error, ierr)`
**Purpose**: Get detailed validation error report  
**Parameters**:
- `bs` (in): `type(back_sett_t)` - Background settings to validate
- `detailed_error` (out): `character(len=*)` - Detailed error report
- `ierr` (out): `integer` - Error status

### `settings_log_error(error_code, context, parameter, log_file, ierr)`
**Purpose**: Log error to file with timestamp  
**Parameters**:
- `error_code` (in): `integer` - Error code
- `context` (in): `character(len=*)` - Error context
- `parameter` (in): `character(len=*)` - Parameter name
- `log_file` (in): `character(len=*)` - Log file path
- `ierr` (out): `integer` - Error status

### `settings_check_error_log_exists(log_file, logged, ierr)`
**Purpose**: Check if error log file exists  
**Parameters**:
- `log_file` (in): `character(len=*)` - Log file path
- `logged` (out): `logical` - File existence flag
- `ierr` (out): `integer` - Error status

### `settings_clear_error_log(log_file, ierr)`
**Purpose**: Clear/delete error log file  
**Parameters**:
- `log_file` (in): `character(len=*)` - Log file path
- `ierr` (out): `integer` - Error status

---

## Utility Procedures

### `settings_deep_copy(src, dst, ierr)`
**Purpose**: Create deep copy of complete settings  
**Parameters**:
- `src` (in): `type(settings_t)` - Source settings
- `dst` (inout): `type(settings_t)` - Destination settings
- `ierr` (out): `integer` - Error status

### `settings_compare(sd1, sd2, is_equal, ierr)`
**Purpose**: Compare two complete settings for equality  
**Parameters**:
- `sd1`, `sd2` (in): `type(settings_t)` - Settings to compare
- `is_equal` (out): `logical` - Comparison result
- `ierr` (out): `integer` - Error status

### `settings_get_antenna(sd)`
**Purpose**: Get pointer to antenna settings  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
**Returns**: `type(antenna_t), pointer` - Pointer to antenna settings

### `settings_get_background(sd)`
**Purpose**: Get pointer to background settings  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure  
**Returns**: `type(back_sett_t), pointer` - Pointer to background settings

### `settings_get_output(sd)`
**Purpose**: Get pointer to output settings  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
**Returns**: `type(output_sett_t), pointer` - Pointer to output settings

### `settings_get_eigmode(sd)`
**Purpose**: Get pointer to eigenmode settings  
**Parameters**:
- `sd` (in): `type(settings_t)` - Settings structure
**Returns**: `type(eigmode_sett_t), pointer` - Pointer to eigenmode settings

---

## Error Codes

All procedures use standard KiLCA error codes:

- `KILCA_SUCCESS` (0) - Operation successful
- `KILCA_ERROR_INVALID_INPUT` - Invalid input parameters
- `KILCA_ERROR_MEMORY` - Memory allocation failure
- `KILCA_ERROR_FILE` - File operation error
- `KILCA_ERROR_INVALID_DATA` - Invalid data encountered
- `KILCA_ERROR_UNKNOWN` - Unknown error

---

## Usage Notes

1. **Always check return codes**: Every procedure returns an error status that should be checked
2. **Memory management**: Always call `settings_destroy` when done with settings
3. **Validation**: Validate settings before use in calculations
4. **Deep copy**: Use provided copy procedures for safe copying of settings
5. **Error handling**: Use error handling procedures for robust error reporting
6. **Path validation**: Ensure paths are non-empty when creating settings
7. **Array management**: Check allocation status before accessing dynamic arrays