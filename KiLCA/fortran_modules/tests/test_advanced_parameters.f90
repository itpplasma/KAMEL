program test_advanced_parameters
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: settings_t, settings_read_namelist
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "==========================================="
    print *, "Testing Advanced Complex Numbers and Arrays [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_advanced_parameter_features()
    
    if (test_status == 0) then
        print *, ""
        print *, "Advanced parameter tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Advanced parameter tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test advanced complex numbers and allocatable array features
    subroutine test_advanced_parameter_features()
        print *, "Testing advanced complex and array parameter features..."
        
        ! Create advanced test configuration file
        call create_advanced_config()
        
        ! This should FAIL - advanced features not fully supported yet
        print *, "Attempting to read advanced parameter configuration..."
        call settings_read_namelist("advanced_config.conf", settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read advanced configuration, error:", ierr
            test_status = test_status + 1
            return
        else
            print *, "SUCCESS: Advanced configuration read correctly"
        end if
        
        ! Test advanced complex number scenarios
        call test_complex_number_precision()
        call test_dynamic_array_sizing()
        call test_complex_array_edge_cases()
        call test_string_parameter_edge_cases()
        
    end subroutine test_advanced_parameter_features
    
    !> Test high-precision complex number handling
    subroutine test_complex_number_precision()
        print *, "Testing high-precision complex number handling..."
        
        ! Test very small and very large complex numbers
        if (abs(real(settings%antenna_settings%flab) - 1.23456789e-12_dp) > 1.0e-20_dp) then
            print *, "FAIL: High-precision real part not handled correctly:", real(settings%antenna_settings%flab)
            test_status = test_status + 1
        else
            print *, "PASS: High-precision real part"
        end if
        
        if (abs(aimag(settings%antenna_settings%flab) - (-9.87654321e10_dp)) > 1.0e2_dp) then
            print *, "FAIL: Large imaginary part not handled correctly:", aimag(settings%antenna_settings%flab)
            test_status = test_status + 1
        else
            print *, "PASS: Large imaginary part"
        end if
        
    end subroutine test_complex_number_precision
    
    !> Test dynamic array sizing based on parameter values
    subroutine test_dynamic_array_sizing()
        print *, "Testing dynamic array sizing..."
        
        ! Test that array sizes match parameter specifications
        if (.not. allocated(settings%eigmode_settings%fstart)) then
            print *, "FAIL: fstart array not allocated"
            test_status = test_status + 1
            return
        end if
        
        ! Size should match Nguess parameter (15 in config)
        if (size(settings%eigmode_settings%fstart) /= settings%eigmode_settings%Nguess) then
            print *, "FAIL: fstart array size doesn't match Nguess:", &
                     size(settings%eigmode_settings%fstart), " vs ", settings%eigmode_settings%Nguess
            test_status = test_status + 1
        else
            print *, "PASS: Dynamic array sizing matches parameter"
        end if
        
        ! Test antenna modes array sizing
        if (.not. allocated(settings%antenna_settings%modes)) then
            print *, "FAIL: antenna modes array not allocated"
            test_status = test_status + 1
            return
        end if
        
        if (size(settings%antenna_settings%modes) /= settings%antenna_settings%dma * 2) then
            print *, "FAIL: modes array size doesn't match dma*2:", &
                     size(settings%antenna_settings%modes), " vs ", settings%antenna_settings%dma * 2
            test_status = test_status + 1
        else
            print *, "PASS: Modes array sizing matches dma parameter"
        end if
        
    end subroutine test_dynamic_array_sizing
    
    !> Test complex array edge cases
    subroutine test_complex_array_edge_cases()
        print *, "Testing complex array edge cases..."
        
        ! Test zero values in complex array
        if (abs(settings%eigmode_settings%fstart(1)) > 1.0e-15_dp) then
            print *, "FAIL: Zero complex value not handled correctly:", settings%eigmode_settings%fstart(1)
            test_status = test_status + 1
        else
            print *, "PASS: Zero complex value"
        end if
        
        ! Test pure imaginary number
        if (abs(real(settings%eigmode_settings%fstart(2))) > 1.0e-15_dp .or. &
            abs(aimag(settings%eigmode_settings%fstart(2)) - 1.0e9_dp) > 1.0e-6_dp) then
            print *, "FAIL: Pure imaginary value not handled correctly:", settings%eigmode_settings%fstart(2)
            test_status = test_status + 1
        else
            print *, "PASS: Pure imaginary value"
        end if
        
        ! Test negative complex number
        if (abs(real(settings%eigmode_settings%fstart(3)) - (-5.5e6_dp)) > 1.0e-6_dp .or. &
            abs(aimag(settings%eigmode_settings%fstart(3)) - (-2.2e5_dp)) > 1.0e-6_dp) then
            print *, "FAIL: Negative complex value not handled correctly:", settings%eigmode_settings%fstart(3)
            test_status = test_status + 1
        else
            print *, "PASS: Negative complex value"
        end if
        
    end subroutine test_complex_array_edge_cases
    
    !> Test string parameter edge cases
    subroutine test_string_parameter_edge_cases()
        print *, "Testing string parameter edge cases..."
        
        ! Test long path with special characters
        if (.not. allocated(settings%background_settings%path2profiles)) then
            print *, "FAIL: path2profiles not allocated"
            test_status = test_status + 1
            return
        end if
        
        if (settings%background_settings%path2profiles /= "/very/long/path/with-special_chars/and.dots/profiles/") then
            print *, "FAIL: Long path with special characters not handled correctly:"
            print *, "Expected: '/very/long/path/with-special_chars/and.dots/profiles/'"
            print *, "Got:      '", settings%background_settings%path2profiles, "'"
            test_status = test_status + 1
        else
            print *, "PASS: Long path with special characters"
        end if
        
        ! Test string with spaces and quotes
        if (.not. allocated(settings%eigmode_settings%fname)) then
            print *, "FAIL: fname not allocated"
            test_status = test_status + 1
            return
        end if
        
        if (settings%eigmode_settings%fname /= "output file with spaces.dat") then
            print *, "FAIL: String with spaces not handled correctly:"
            print *, "Expected: 'output file with spaces.dat'"
            print *, "Got:      '", settings%eigmode_settings%fname, "'"
            test_status = test_status + 1
        else
            print *, "PASS: String with spaces"
        end if
        
    end subroutine test_string_parameter_edge_cases
    
    !> Create advanced test configuration with edge cases
    subroutine create_advanced_config()
        integer :: unit
        integer :: iostat
        
        open(newunit=unit, file="advanced_config.conf", status="replace", action="write", iostat=iostat)
        if (iostat /= 0) then
            print *, "ERROR: Could not create advanced config file, iostat =", iostat
            return
        end if
        
        write(unit, '(a)') "! Advanced Parameter Configuration - Edge cases and precision tests"
        write(unit, '(a)') ""
        
        ! Antenna section with high precision complex number
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 0.0"
        write(unit, '(a)') "  I0 = 1.0e13"
        write(unit, '(a)') "  flab = (1.23456789e-12, -9.87654321e10)"
        write(unit, '(a)') "  dma = 5"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1, 3, 2, -2, 3, 4, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        
        ! Background section with long path and special characters
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.69"
        write(unit, '(a)') "  rp = 70.0"
        write(unit, '(a)') "  B0 = 23176.46"
        write(unit, '(a)') "  path2profiles = '/very/long/path/with-special_chars/and.dots/profiles/'"
        write(unit, '(a)') "  calc_back = 1"
        write(unit, '(a)') "  flag_back = 'homogeneous'"
        write(unit, '(a)') "  N = 9"
        write(unit, '(a)') "  V_gal_sys = 1.0e9"
        write(unit, '(a)') "  V_scale = 1.0"
        write(unit, '(a)') "  m_i = 2.0"
        write(unit, '(a)') "  zele = 1.0"
        write(unit, '(a)') "  zion = 1.0"
        write(unit, '(a)') "  flag_debug_bg = 0"
        write(unit, '(a)') "  mass = 2.0, 1.0, 4.0"
        write(unit, '(a)') "  charge = 1.0, -1.0, 2.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        
        ! Output section
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 2"
        write(unit, '(a)') "  flag_emfield = 2"
        write(unit, '(a)') "  flag_additional = 2"
        write(unit, '(a)') "  flag_dispersion = 0"
        write(unit, '(a)') "  flag_debug_out = 0"
        write(unit, '(a)') "  num_quants = 8"
        write(unit, '(a)') "  flag_quants = 1, 1, 1, 1, 1, 1, 1, 0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        
        ! Eigenmode section with large complex array and string with spaces
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  fname = 'output file with spaces.dat'"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  rdim = 100"
        write(unit, '(a)') "  rfmin = 0.0"
        write(unit, '(a)') "  rfmax = 1.0e6"
        write(unit, '(a)') "  idim = 100"
        write(unit, '(a)') "  ifmin = -1.0e5"
        write(unit, '(a)') "  ifmax = 1.0e5"
        write(unit, '(a)') "  stop_flag = 0"
        write(unit, '(a)') "  eps_res = 1.0e-6"
        write(unit, '(a)') "  eps_abs = 1.0e-8"
        write(unit, '(a)') "  eps_rel = 1.0e-6"
        write(unit, '(a)') "  delta = 1.0e-3"
        write(unit, '(a)') "  test_roots = 0"
        write(unit, '(a)') "  flag_debug_eig = 0"
        write(unit, '(a)') "  Nguess = 15"
        write(unit, '(a)') "  kmin = 1"
        write(unit, '(a)') "  kmax = 10"
        write(unit, '(a)') "  n_zeros = 10"
        write(unit, '(a)') "  use_winding = 0"
        ! Large complex array with edge cases: zero, pure imaginary, negative values
        write(unit, '(a)') "  fstart = (0.0, 0.0), (0.0, 1.0e9), (-5.5e6, -2.2e5), (1.0e6, 0.0), (2.0e6, 1.0e5), (3.0e6, 2.0e5), (4.0e6, 3.0e5), (5.0e6, 4.0e5), (6.0e6, 5.0e5), (7.0e6, -1.0e5), (8.0e6, -2.0e5), (9.0e6, -3.0e5), (1.0e7, -4.0e5), (1.1e7, -5.0e5), (1.2e7, 0.0)"
        write(unit, '(a)') "/"
        
        close(unit)
        print *, "Advanced config file created successfully"
    end subroutine create_advanced_config

end program test_advanced_parameters