# KiLCA Settings Module Usage Guide

## Overview

The `kilca_settings_m` module provides comprehensive settings management for KiLCA Fortran implementation. It handles four main categories of settings:

1. **Antenna Settings** (`antenna_t`) - Antenna parameters and modes
2. **Background Settings** (`back_sett_t`) - Plasma background and machine parameters  
3. **Output Settings** (`output_sett_t`) - Output control flags and quantities
4. **Eigenmode Settings** (`eigmode_sett_t`) - Eigenmode search parameters

## Quick Start

### Basic Usage Pattern

```fortran
program example_settings_usage
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: sd
    integer :: ierr
    
    ! Create settings with project path
    call settings_create(sd, "/path/to/project/", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "Error creating settings:", ierr
        stop 1
    end if
    
    ! Initialize with defaults
    call settings_initialize_defaults(sd, "/path/to/project/", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "Error initializing defaults:", ierr
        stop 1
    end if
    
    ! Modify settings as needed
    sd%antenna_settings%ra = 50.0_dp  ! Set antenna radius
    sd%background_settings%B0 = 25000.0_dp  ! Set magnetic field
    
    ! Validate all settings
    call settings_validate_complete(sd, ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "Settings validation failed:", ierr
        stop 1
    end if
    
    ! Print current settings
    call settings_print_all(sd, ierr)
    
    ! Clean up
    call settings_destroy(sd, ierr)
    
end program example_settings_usage
```

## Detailed Usage Instructions

### 1. Settings Creation and Destruction

#### Creating Settings
```fortran
type(settings_t), pointer :: sd
integer :: ierr

! Create with project path
call settings_create(sd, "/home/user/plasma_run/", ierr)
```

#### Destroying Settings  
```fortran
! Always destroy settings to free memory
call settings_destroy(sd, ierr)
```

### 2. Initialization Strategies

#### Default Initialization
```fortran
! Initialize all subsettings with scientific defaults
call settings_initialize_defaults(sd, "/path/to/project/", ierr)
```

#### Custom Initialization
```fortran
! Initialize individual components with custom parameters
call antenna_initialize_custom(sd%antenna_settings, &
                               ra=45.0_dp, wa=2.0_dp, I0=500.0_dp, ierr=ierr)
```

#### Reading from Files
```fortran
! Read all settings from configuration files
call settings_read_all(sd, ierr)
```

### 3. Antenna Settings (`antenna_t`)

#### Setting Antenna Parameters
```fortran
! Direct assignment
sd%antenna_settings%ra = 50.0_dp        ! Antenna radius (cm)
sd%antenna_settings%wa = 1.5_dp         ! Current layer width
sd%antenna_settings%I0 = 1000.0_dp      ! Current (statamp)
sd%antenna_settings%flab = (1.0e5_dp, 0.0_dp)  ! Frequency (Hz)

! Using procedure
call antenna_set_parameters(sd%antenna_settings, &
    ra=50.0_dp, wa=1.5_dp, I0=1000.0_dp, ierr=ierr)
```

#### Setting Antenna Modes
```fortran
! Allocate and set mode array
sd%antenna_settings%dma = 4
allocate(sd%antenna_settings%modes(sd%antenna_settings%dma))
sd%antenna_settings%modes = [1, 1, 2, 1]  ! [m1, n1, m2, n2]
```

#### Antenna Validation and Output
```fortran
logical :: is_valid
character(len=1024) :: error_msg

! Validate antenna settings
call antenna_validate(sd%antenna_settings, is_valid, error_msg, ierr)
if (.not. is_valid) then
    print *, "Antenna validation error:", trim(error_msg)
end if

! Print antenna settings
call antenna_print_settings(sd%antenna_settings, ierr)
```

### 4. Background Settings (`back_sett_t`)

#### Setting Machine Parameters
```fortran
! ASDEX Upgrade typical parameters
sd%background_settings%rtor = 625.0_dp      ! Major radius (cm)
sd%background_settings%rp = 200.0_dp        ! Plasma radius (cm)  
sd%background_settings%B0 = 25000.0_dp      ! Magnetic field (G)
```

#### Setting Plasma Parameters
```fortran
sd%background_settings%m_i = 2.0_dp         ! Ion mass (proton masses)
sd%background_settings%zele = 1.0_dp        ! Electron collision coeff
sd%background_settings%zion = 1.0_dp        ! Ion collision coeff
sd%background_settings%V_scale = 1.0_dp     ! Velocity scale factor
```

#### Setting Background Calculation Mode
```fortran
! Set background calculation flag
call back_sett_set_calc_flag(sd%background_settings, 1, ierr)

! Set background type
if (allocated(sd%background_settings%flag_back)) &
    deallocate(sd%background_settings%flag_back)
sd%background_settings%flag_back = "normal"  ! or "homogeneous"
```

### 5. Output Settings (`output_sett_t`)

#### Setting Output Flags
```fortran
sd%output_settings%flag_background = 1      ! Compute background
sd%output_settings%flag_emfield = 1         ! Compute EM fields
sd%output_settings%flag_additional = 0      ! Skip additional quantities
sd%output_settings%flag_dispersion = 0      ! Skip dispersion

! Using procedure
call output_sett_set_flags(sd%output_settings, &
    flag_background=1, flag_emfield=1, ierr=ierr)
```

#### Setting Quantity Flags
```fortran
! Set number of quantities and allocate array
sd%output_settings%num_quants = 3
if (allocated(sd%output_settings%flag_quants)) &
    deallocate(sd%output_settings%flag_quants)
allocate(sd%output_settings%flag_quants(3))
sd%output_settings%flag_quants = [1, 1, 0]  ! Compute first two quantities
```

### 6. Eigenmode Settings (`eigmode_sett_t`)

#### Setting Search Parameters
```fortran
sd%eigmode_settings%search_flag = 1         ! Enable search
sd%eigmode_settings%rdim = 200              ! Real dimension
sd%eigmode_settings%idim = 200              ! Imaginary dimension
sd%eigmode_settings%rfmin = 0.0_dp          ! Real frequency min
sd%eigmode_settings%rfmax = 1.0e6_dp        ! Real frequency max
sd%eigmode_settings%ifmin = -1.0e5_dp       ! Imaginary frequency min
sd%eigmode_settings%ifmax = 1.0e5_dp        ! Imaginary frequency max
```

#### Setting Solver Parameters  
```fortran
sd%eigmode_settings%eps_res = 1.0e-8_dp     ! Residual tolerance
sd%eigmode_settings%eps_abs = 1.0e-10_dp    ! Absolute tolerance
sd%eigmode_settings%eps_rel = 1.0e-8_dp     ! Relative tolerance
sd%eigmode_settings%delta = 1.0e-8_dp       ! Delta parameter
```

#### Setting Output File
```fortran
if (allocated(sd%eigmode_settings%fname)) &
    deallocate(sd%eigmode_settings%fname)
sd%eigmode_settings%fname = "eigenmode_results.dat"
```

## Advanced Usage

### 7. Settings Validation

#### Individual Component Validation
```fortran
logical :: is_valid
character(len=1024) :: error_msg
integer :: ierr

! Validate individual components
call antenna_validate(sd%antenna_settings, is_valid, error_msg, ierr)
call back_sett_validate(sd%background_settings, is_valid, error_msg, ierr)
call output_sett_validate(sd%output_settings, is_valid, error_msg, ierr)
call eigmode_sett_validate(sd%eigmode_settings, is_valid, error_msg, ierr)
```

#### Complete Settings Validation
```fortran
! Validate all settings and check consistency
call settings_validate_complete(sd, ierr)
if (ierr /= KILCA_SUCCESS) then
    print *, "Settings validation failed"
end if

! Validate with context information
character(len=2048) :: context_error
call settings_validate_with_context(sd, "main simulation", context_error, ierr)
```

#### Detailed Error Reporting
```fortran
character(len=4096) :: detailed_errors

! Get detailed validation errors for background settings
call settings_get_detailed_validation_errors(sd%background_settings, &
                                             detailed_errors, ierr)
if (ierr /= KILCA_SUCCESS .and. len_trim(detailed_errors) > 0) then
    print *, "Detailed errors:"
    print *, trim(detailed_errors)
end if
```

### 8. Settings Copying and Comparison

#### Deep Copy
```fortran
type(settings_t), pointer :: sd_copy
integer :: ierr

! Create copy structure
call settings_create(sd_copy, sd%path2project, ierr)

! Perform deep copy
call settings_deep_copy(sd, sd_copy, ierr)
```

#### Comparison
```fortran
logical :: is_equal

! Compare two settings structures
call settings_compare(sd, sd_copy, is_equal, ierr)
if (is_equal) then
    print *, "Settings are identical"
else
    print *, "Settings differ"
end if
```

### 9. Error Handling

#### Error Message Formatting
```fortran
character(len=1024) :: error_msg

! Format detailed error message
call settings_format_error_message(KILCA_ERROR_INVALID_INPUT, &
                                  "antenna validation", &
                                  "ra parameter", error_msg)
print *, trim(error_msg)
```

#### Error Recovery
```fortran
logical :: recovered

! Attempt to recover from invalid settings
call settings_attempt_recovery(sd%antenna_settings, recovered, ierr)
if (recovered) then
    print *, "Successfully recovered from invalid settings"
else
    print *, "Could not recover, manual intervention required"
end if
```

#### Error Logging
```fortran
! Log error to file
call settings_log_error(KILCA_ERROR_INVALID_INPUT, &
                       "simulation startup", &
                       "background B0 field", &
                       "error.log", ierr)
```

### 10. File I/O

#### Reading Settings from Files
```fortran
! Read individual components
call antenna_read_settings(sd%antenna_settings, sd%path2project, ierr)
call back_sett_read_settings(sd%background_settings, sd%path2project, ierr)

! Read all settings
call settings_read_all(sd, ierr)
```

#### Printing Settings
```fortran
! Print to stdout
call settings_print_all(sd, ierr)

! Print individual components
call antenna_print_settings(sd%antenna_settings, ierr)
call back_sett_print_settings(sd%background_settings, ierr)

! Print to file
open(unit=10, file="settings_output.txt")
call antenna_print_settings_to_unit(sd%antenna_settings, 10, ierr)
call back_sett_print_settings_to_unit(sd%background_settings, 10, ierr)
close(10)
```

## Error Codes

The module uses the following error codes from `kilca_types_m`:

- `KILCA_SUCCESS` - Operation completed successfully
- `KILCA_ERROR_INVALID_INPUT` - Invalid input parameters
- `KILCA_ERROR_MEMORY` - Memory allocation failure
- `KILCA_ERROR_FILE` - File operation error
- `KILCA_ERROR_INVALID_DATA` - Invalid data encountered
- `KILCA_ERROR_UNKNOWN` - Unknown error

## Best Practices

### 1. Always Check Return Codes
```fortran
call some_settings_procedure(sd, ierr)
if (ierr /= KILCA_SUCCESS) then
    ! Handle error appropriately
    print *, "Error in settings procedure:", ierr
    ! Consider cleanup and exit
end if
```

### 2. Validate Before Use
```fortran
! Always validate settings before using them in calculations
call settings_validate_complete(sd, ierr)
if (ierr /= KILCA_SUCCESS) then
    print *, "Settings are invalid, cannot proceed"
    call settings_destroy(sd, ierr)
    stop 1
end if
```

### 3. Proper Memory Management
```fortran
! Always deallocate arrays before reassigning
if (allocated(sd%antenna_settings%modes)) then
    deallocate(sd%antenna_settings%modes)
end if

! Always destroy settings when done
call settings_destroy(sd, ierr)
```

### 4. Use Initialization Procedures
```fortran
! Use provided initialization procedures rather than manual assignment
call settings_initialize_defaults(sd, path, ierr)
! Rather than manually setting each field
```

### 5. Handle String Allocations Carefully
```fortran
! Properly handle allocatable strings
if (allocated(sd%background_settings%flag_back)) then
    deallocate(sd%background_settings%flag_back)
end if
sd%background_settings%flag_back = "normal"
```

## Common Pitfalls

1. **Forgetting to validate settings** - Always validate before use
2. **Memory leaks** - Always call `settings_destroy` when done
3. **Array bounds** - Check `dma` before accessing `modes` array
4. **String allocation** - Deallocate before reassigning allocatable strings
5. **Error handling** - Always check return codes from procedures

## Examples

See the test files for comprehensive usage examples:
- `tests/test_settings_initialization.f90` - Initialization examples
- `tests/test_settings_validation.f90` - Validation examples  
- `tests/test_settings_copy_compare.f90` - Copy and comparison examples
- `tests/test_settings_error_handling.f90` - Error handling examples

## Support

For questions or issues with the settings module, refer to:
1. This documentation
2. Test files for working examples
3. Module source code comments
4. KiLCA project documentation