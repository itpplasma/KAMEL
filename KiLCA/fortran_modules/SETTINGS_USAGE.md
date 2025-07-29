# KiLCA Settings Module - Usage Guide

This document provides practical examples and usage patterns for the `kilca_settings_m` module.

## Table of Contents

1. [Basic Usage](#basic-usage)
2. [Advanced Usage Patterns](#advanced-usage-patterns)
3. [Error Handling](#error-handling)
4. [Validation and Debugging](#validation-and-debugging)
5. [Performance Considerations](#performance-considerations)
6. [Common Patterns](#common-patterns)
7. [Troubleshooting](#troubleshooting)

## Basic Usage

### Creating and Initializing Settings

The most common usage pattern is to create settings with defaults and then customize as needed:

```fortran
program basic_settings_example
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: settings
    integer :: ierr
    
    ! Create settings with defaults
    call settings_initialize_defaults(settings, "/path/to/project", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "Failed to initialize settings"
        stop 1
    end if
    
    ! Customize antenna settings
    settings%antenna_settings%ra = 10.5_dp      ! Antenna radius (cm)
    settings%antenna_settings%wa = 1.2_dp       ! Current layer width (cm)
    settings%antenna_settings%I0 = 1000.0_dp    ! Antenna current (statamps)
    settings%antenna_settings%flab = cmplx(50.0e6_dp, 2000.0_dp)  ! Frequency (Hz)
    
    ! Set up modes array
    settings%antenna_settings%dma = 2
    allocate(settings%antenna_settings%modes(4))  ! 2*dma
    settings%antenna_settings%modes = [1, 1, 2, 2]  ! (m,n) pairs
    
    ! Customize background settings
    settings%background_settings%rtor = 625.0_dp    ! Major radius (cm)
    settings%background_settings%rp = 200.0_dp      ! Minor radius (cm)
    settings%background_settings%B0 = 25000.0_dp    ! Toroidal field (G)
    
    ! Print all settings
    call settings_print_all(settings, ierr)
    
    ! Clean up
    call settings_destroy(settings, ierr)
end program basic_settings_example
```

### Working with Individual Components

You can also work with individual settings components:

```fortran
program component_example
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(antenna_t) :: antenna
    type(back_sett_t) :: background
    integer :: ierr
    logical :: is_valid
    character(len=1024) :: error_msg
    
    ! Initialize antenna with defaults
    call antenna_initialize_defaults(antenna, ierr)
    if (ierr /= KILCA_SUCCESS) stop 1
    
    ! Customize specific parameters
    call antenna_initialize_custom(antenna, ra=15.0_dp, wa=2.0_dp, I0=1500.0_dp, ierr=ierr)
    if (ierr /= KILCA_SUCCESS) stop 1
    
    ! Validate settings
    call antenna_validate(antenna, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "Antenna validation failed: ", trim(error_msg)
        stop 1
    end if
    
    ! Print antenna settings
    call antenna_print_settings(antenna, ierr)
    
    ! Initialize background settings
    call back_sett_initialize_defaults(background, ierr)
    call back_sett_print_settings(background, ierr)
    
    ! Clean up arrays if allocated
    if (allocated(antenna%modes)) deallocate(antenna%modes)
    if (allocated(background%flag_back)) deallocate(background%flag_back)
end program component_example
```

## Advanced Usage Patterns

### Complete Plasma Configuration

Here's a more advanced example setting up a complete ASDEX Upgrade-like plasma configuration:

```fortran
program advanced_plasma_setup
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: settings
    integer :: ierr, i
    logical :: is_valid
    character(len=1024) :: error_msg
    
    ! Initialize with project path
    call settings_create(settings, "/home/user/KAMEL_runs/shot_12345", ierr)
    if (ierr /= KILCA_SUCCESS) stop 1
    
    ! === ANTENNA CONFIGURATION ===
    ! ICRF antenna at plasma edge
    settings%antenna_settings%ra = 18.5_dp       ! Close to plasma boundary
    settings%antenna_settings%wa = 0.8_dp        ! Narrow current layer
    settings%antenna_settings%I0 = 800.0_dp      ! Moderate current
    settings%antenna_settings%flab = cmplx(30.0e6_dp, 0.0_dp)  ! 30 MHz
    settings%antenna_settings%flag_debug = 1     ! Enable debugging
    
    ! Set up toroidal mode numbers for ICRF heating
    settings%antenna_settings%dma = 3
    allocate(settings%antenna_settings%modes(6))
    settings%antenna_settings%modes = [0, 15, 0, 16, 0, 17]  ! n=15,16,17
    
    ! === BACKGROUND PLASMA ===
    ! ASDEX Upgrade parameters
    settings%background_settings%rtor = 165.0_dp    ! ASDEX major radius
    settings%background_settings%rp = 50.0_dp       ! ASDEX minor radius  
    settings%background_settings%B0 = 25000.0_dp    ! 2.5 T toroidal field
    settings%background_settings%calc_back = 1      ! Calculate profiles
    settings%background_settings%N = 7              ! Higher order splines
    
    ! Plasma composition (D plasma)
    allocate(settings%background_settings%mass(2))
    allocate(settings%background_settings%charge(2))
    settings%background_settings%mass = [1.0_dp, 2.0_dp]    ! e-, D+
    settings%background_settings%charge = [-1.0_dp, 1.0_dp]  ! charges
    
    ! Set background profile path
    settings%background_settings%path2profiles = &
        trim(settings%path2project) // "/profiles/"
    
    ! === OUTPUT CONFIGURATION ===
    settings%output_settings%flag_background = 1
    settings%output_settings%flag_emfield = 1
    settings%output_settings%flag_additional = 1
    settings%output_settings%num_quants = 5
    
    allocate(settings%output_settings%flag_quants(5))
    settings%output_settings%flag_quants = [1, 1, 0, 1, 1]  ! Select quantities
    
    ! === EIGENMODE SETTINGS ===
    settings%eigmode_settings%search_flag = 1       ! Enable eigenmode search
    settings%eigmode_settings%rdim = 200           ! High resolution
    settings%eigmode_settings%idim = 200
    settings%eigmode_settings%rfmin = 20.0e6_dp    ! Around antenna frequency
    settings%eigmode_settings%rfmax = 40.0e6_dp
    settings%eigmode_settings%ifmin = -1.0e3_dp    ! Small damping range
    settings%eigmode_settings%ifmax = 1.0e3_dp
    settings%eigmode_settings%n_zeros = 20         ! Find many modes
    
    ! Set output filename
    settings%eigmode_settings%fname = &
        trim(settings%path2project) // "/eigenmode_results.dat"
    
    ! Provide initial guesses near antenna frequency
    settings%eigmode_settings%Nguess = 3
    allocate(settings%eigmode_settings%fstart(3))
    settings%eigmode_settings%fstart = [ &
        cmplx(29.5e6_dp, 100.0_dp), &
        cmplx(30.0e6_dp, 50.0_dp), &
        cmplx(30.5e6_dp, 150.0_dp) &
    ]
    
    ! === VALIDATION ===
    call settings_validate_complete(settings, is_valid, error_msg, ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "Validation failed with error code:", ierr
        stop 1
    end if
    
    if (.not. is_valid) then
        print *, "Settings validation failed:"
        print *, trim(error_msg)
        stop 1
    end if
    
    ! Check consistency between components  
    call settings_validate_consistency(settings, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "Consistency check failed:"
        print *, trim(error_msg)
        stop 1
    end if
    
    print *, "=== PLASMA CONFIGURATION SUMMARY ==="
    call settings_print_all(settings, ierr)
    
    ! Clean up
    call settings_destroy(settings, ierr)
    
end program advanced_plasma_setup
```

## Error Handling

### Comprehensive Error Handling Pattern

```fortran
program error_handling_example
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: settings
    integer :: ierr
    character(len=1024) :: error_msg
    character(len=256) :: error_name
    logical :: recovered
    
    ! Attempt to create settings
    call settings_create(settings, "", ierr)  ! Empty path - should fail
    
    if (ierr /= KILCA_SUCCESS) then
        ! Get error information
        call settings_get_error_name(ierr, error_name, ierr)
        call settings_format_error_message(ierr, "settings creation", "empty path", error_msg)
        
        print *, "Error occurred:"
        print *, "  Code: ", ierr
        print *, "  Name: ", trim(error_name)  
        print *, "  Message: ", trim(error_msg)
        
        ! Log error
        call settings_log_error(ierr, "settings creation", "empty path", "error.log", ierr)
    end if
    
    ! Create with valid path
    call settings_create(settings, "/tmp/test", ierr)
    if (ierr /= KILCA_SUCCESS) stop 1
    
    ! Intentionally create invalid antenna state
    settings%antenna_settings%ra = -10.0_dp      ! Invalid negative radius
    settings%antenna_settings%dma = 2            ! But no modes array
    
    ! Try to recover
    call settings_attempt_recovery(settings%antenna_settings, recovered, ierr)
    if (recovered) then
        print *, "Successfully recovered from invalid antenna state"
        print *, "  Fixed ra =", settings%antenna_settings%ra
        if (allocated(settings%antenna_settings%modes)) then
            print *, "  Allocated modes array with size", size(settings%antenna_settings%modes)
        end if
    else
        print *, "Could not recover from invalid state"
    end if
    
    call settings_destroy(settings, ierr)
    
end program error_handling_example
```

### Error Logging and Recovery

```fortran
subroutine robust_settings_setup(settings, project_path, success)
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer, intent(out) :: settings
    character(len=*), intent(in) :: project_path
    logical, intent(out) :: success
    
    integer :: ierr, attempt
    logical :: is_valid, recovered
    character(len=1024) :: error_msg, detailed_error
    character(len=256) :: log_file
    
    success = .false.
    log_file = trim(project_path) // "/setup_errors.log"
    
    ! Clear any previous error log
    call settings_clear_error_log(log_file, ierr)
    
    do attempt = 1, 3  ! Try up to 3 times
        ! Attempt to create and initialize
        call settings_initialize_defaults(settings, project_path, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            call settings_log_error(ierr, "initialization", "attempt", log_file, ierr)
            cycle  ! Try again
        end if
        
        ! Validate the settings
        call settings_validate_complete(settings, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call settings_log_error(ierr, "validation", "complete check", log_file, ierr)
            call settings_destroy(settings, ierr)
            cycle
        end if
        
        if (.not. is_valid) then
            ! Try to get detailed error information
            call settings_get_detailed_validation_errors( &
                settings%background_settings, detailed_error, ierr)
            
            ! Log the validation failure
            call settings_log_error(KILCA_ERROR_INVALID_INPUT, &
                "validation", trim(error_msg), log_file, ierr)
                
            ! Try recovery
            call settings_attempt_recovery(settings%antenna_settings, recovered, ierr)
            if (recovered) then
                ! Re-validate after recovery
                call settings_validate_complete(settings, is_valid, error_msg, ierr)
                if (is_valid) then
                    success = .true.
                    exit
                end if
            end if
            
            call settings_destroy(settings, ierr)
            cycle
        end if
        
        ! Success!
        success = .true.
        exit
    end do
    
    if (.not. success) then
        print *, "Failed to set up settings after", attempt, "attempts"
        print *, "Check error log:", trim(log_file)
    end if
    
end subroutine robust_settings_setup
```

## Validation and Debugging

### Comprehensive Validation

```fortran
program validation_example
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: settings
    integer :: ierr
    logical :: is_valid
    character(len=2048) :: error_msg, detailed_error
    
    call settings_initialize_defaults(settings, "/tmp/test", ierr)
    if (ierr /= KILCA_SUCCESS) stop 1
    
    ! Set up some potentially problematic values
    settings%antenna_settings%ra = 25.0_dp       ! Larger than plasma
    settings%background_settings%rp = 20.0_dp    ! Small plasma
    settings%background_settings%rtor = 30.0_dp  ! Small torus
    
    ! Individual component validation
    print *, "=== INDIVIDUAL COMPONENT VALIDATION ==="
    
    call antenna_validate(settings%antenna_settings, is_valid, error_msg, ierr)
    print *, "Antenna valid:", is_valid
    if (.not. is_valid) print *, "  Error:", trim(error_msg)
    
    call back_sett_validate(settings%background_settings, is_valid, error_msg, ierr)
    print *, "Background valid:", is_valid  
    if (.not. is_valid) print *, "  Error:", trim(error_msg)
    
    ! Complete validation
    print *, ""
    print *, "=== COMPLETE VALIDATION ==="
    call settings_validate_complete(settings, is_valid, error_msg, ierr)
    print *, "Complete settings valid:", is_valid
    if (.not. is_valid) print *, "  Error:", trim(error_msg)
    
    ! Consistency validation
    print *, ""
    print *, "=== CONSISTENCY VALIDATION ==="
    call settings_validate_consistency(settings, is_valid, error_msg, ierr)
    print *, "Settings consistent:", is_valid
    if (.not. is_valid) print *, "  Error:", trim(error_msg)
    
    ! Detailed error analysis
    print *, ""
    print *, "=== DETAILED ERROR ANALYSIS ==="
    call settings_get_detailed_validation_errors( &
        settings%background_settings, detailed_error, ierr)
    print *, "Detailed errors:", trim(detailed_error)
    
    call settings_destroy(settings, ierr)
    
end program validation_example
```

## Performance Considerations

### Efficient Settings Management

```fortran
program performance_example
    use kilca_types_m  
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: template_settings, run_settings
    integer :: ierr, i
    
    ! Create a template settings once
    call settings_initialize_defaults(template_settings, "/project/base", ierr)
    if (ierr /= KILCA_SUCCESS) stop 1
    
    ! Configure common parameters
    template_settings%background_settings%rtor = 165.0_dp
    template_settings%background_settings%rp = 50.0_dp
    template_settings%background_settings%B0 = 25000.0_dp
    
    ! For multiple runs, copy from template
    do i = 1, 100
        ! Deep copy from template (faster than re-initializing)
        call settings_deep_copy(template_settings, run_settings, ierr)
        if (ierr /= KILCA_SUCCESS) stop 1
        
        ! Modify run-specific parameters
        run_settings%antenna_settings%I0 = 500.0_dp + i * 10.0_dp
        run_settings%antenna_settings%flab = cmplx(30.0e6_dp + i * 1.0e5_dp, 0.0_dp)
        
        ! Use settings for calculation...
        ! (simulation code would go here)
        
        ! Clean up this run
        call settings_destroy(run_settings, ierr)
    end do
    
    ! Clean up template
    call settings_destroy(template_settings, ierr)
    
end program performance_example
```

## Common Patterns

### Pattern 1: Parameter Sweep

```fortran
subroutine antenna_frequency_sweep(base_settings, frequencies, n_freq)
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer, intent(in) :: base_settings
    real(dp), intent(in) :: frequencies(:)
    integer, intent(in) :: n_freq
    
    type(settings_t), pointer :: sweep_settings
    integer :: i, ierr
    character(len=256) :: filename
    
    do i = 1, n_freq
        ! Copy base settings
        call settings_deep_copy(base_settings, sweep_settings, ierr)
        if (ierr /= KILCA_SUCCESS) cycle
        
        ! Modify frequency
        sweep_settings%antenna_settings%flab = cmplx(frequencies(i), 0.0_dp)
        
        ! Update output filename
        write(filename, '("freq_sweep_",I0,".dat")') i
        sweep_settings%eigmode_settings%fname = filename
        
        ! Run calculation (placeholder)
        call run_calculation(sweep_settings)
        
        call settings_destroy(sweep_settings, ierr)
    end do
    
end subroutine antenna_frequency_sweep
```

### Pattern 2: Configuration Variants

```fortran
subroutine create_configuration_variants(base_path)
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    character(len=*), intent(in) :: base_path
    
    type(settings_t), pointer :: settings
    integer :: ierr
    
    ! High-performance configuration
    call settings_initialize_defaults(settings, trim(base_path)//"/high_perf", ierr)
    settings%background_settings%N = 7                    ! High-order splines
    settings%eigmode_settings%rdim = 500                  ! High resolution
    settings%eigmode_settings%idim = 500
    settings%output_settings%flag_additional = 1         ! All diagnostics
    call save_configuration(settings, "high_performance.json")
    call settings_destroy(settings, ierr)
    
    ! Fast configuration
    call settings_initialize_defaults(settings, trim(base_path)//"/fast", ierr)
    settings%background_settings%N = 3                    ! Low-order splines
    settings%eigmode_settings%rdim = 50                   ! Low resolution
    settings%eigmode_settings%idim = 50
    settings%output_settings%flag_additional = 0         ! Minimal diagnostics
    call save_configuration(settings, "fast.json")
    call settings_destroy(settings, ierr)
    
end subroutine create_configuration_variants
```

## Troubleshooting

### Common Issues and Solutions

1. **Memory Allocation Failures**
   ```fortran
   ! Check available memory before large allocations
   if (settings%eigmode_settings%rdim * settings%eigmode_settings%idim > 1000000) then
       print *, "Warning: Large frequency mesh may cause memory issues"
   end if
   ```

2. **Validation Failures**
   ```fortran
   ! Always validate after manual modifications
   call settings_validate_consistency(settings, is_valid, error_msg, ierr)
   if (.not. is_valid) then
       print *, "Fix required:", trim(error_msg)
       ! Apply fixes...
   end if
   ```

3. **File Path Issues**
   ```fortran
   ! Ensure directories exist
   logical :: dir_exists
   inquire(file=trim(settings%path2project), exist=dir_exists)
   if (.not. dir_exists) then
       call system("mkdir -p " // trim(settings%path2project))
   end if
   ```

4. **Array Size Mismatches**
   ```fortran
   ! Always check array sizes before allocation
   if (settings%antenna_settings%dma > 0) then
       if (allocated(settings%antenna_settings%modes)) then
           if (size(settings%antenna_settings%modes) /= 2 * settings%antenna_settings%dma) then
               deallocate(settings%antenna_settings%modes)
               allocate(settings%antenna_settings%modes(2 * settings%antenna_settings%dma))
           end if
       else
           allocate(settings%antenna_settings%modes(2 * settings%antenna_settings%dma))
       end if
   end if
   ```

### Debugging Tips

1. **Enable Debug Output**
   ```fortran
   settings%antenna_settings%flag_debug = 1
   settings%background_settings%flag_debug = 1
   settings%output_settings%flag_debug = 1
   settings%eigmode_settings%flag_debug = 1
   ```

2. **Use Validation Context**
   ```fortran
   call settings_validate_with_context(settings, "main initialization", error_msg, ierr)
   ```

3. **Check Error Logs**
   ```fortran
   logical :: log_exists
   call settings_check_error_log_exists("debug.log", log_exists, ierr)
   if (log_exists) then
       print *, "Check debug.log for detailed error information"
   end if
   ```

## See Also

- [PROCEDURE_REFERENCE.md](PROCEDURE_REFERENCE.md) - Complete procedure documentation
- [examples/](examples/) - Working code examples
- KiLCA main documentation for physics background