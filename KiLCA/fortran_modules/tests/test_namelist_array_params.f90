program test_namelist_array_params
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: settings_t, settings_read_namelist
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "==========================================="
    print *, "Testing Namelist Array Parameters [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_array_parameters_reading()
    
    if (test_status == 0) then
        print *, ""
        print *, "Array parameter namelist tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Array parameter namelist tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test array parameter reading functionality
    subroutine test_array_parameters_reading()
        integer :: i
        
        print *, "Testing array parameter reading functionality..."
        
        ! Create test configuration file with array parameters
        call create_test_array_config()
        
        ! This should FAIL - array parameters not supported yet
        print *, "Attempting to read namelist settings with arrays..."
        call settings_read_namelist("test_array_config.conf", settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read array parameter settings, error:", ierr
            test_status = test_status + 1
            return
        else
            print *, "SUCCESS: Array parameter settings read correctly"
        end if
        
        ! Verify output flag_quants array read correctly
        if (.not. allocated(settings%output_settings%flag_quants)) then
            print *, "FAIL: Output flag_quants array not allocated"
            test_status = test_status + 1
            return
        end if
        
        if (size(settings%output_settings%flag_quants) /= 8) then
            print *, "FAIL: Output flag_quants array wrong size:", size(settings%output_settings%flag_quants)
            test_status = test_status + 1
            return
        end if
        
        ! Test specific array values
        if (settings%output_settings%flag_quants(1) /= 1) then
            print *, "FAIL: flag_quants(1) wrong value:", settings%output_settings%flag_quants(1)
            test_status = test_status + 1
        else
            print *, "PASS: flag_quants(1) read correctly"
        end if
        
        if (settings%output_settings%flag_quants(8) /= 0) then
            print *, "FAIL: flag_quants(8) wrong value:", settings%output_settings%flag_quants(8)
            test_status = test_status + 1
        else
            print *, "PASS: flag_quants(8) read correctly"
        end if
        
        ! Verify background mass and charge arrays (if they exist in config)
        if (allocated(settings%background_settings%mass)) then
            if (size(settings%background_settings%mass) /= 3) then
                print *, "FAIL: Background mass array wrong size:", size(settings%background_settings%mass)
                test_status = test_status + 1
            else
                print *, "PASS: Background mass array size correct"
            end if
            
            if (abs(settings%background_settings%mass(1) - 2.0_dp) > 1.0e-12_dp) then
                print *, "FAIL: mass(1) wrong value:", settings%background_settings%mass(1)
                test_status = test_status + 1
            else
                print *, "PASS: mass(1) read correctly"
            end if
        else
            print *, "FAIL: Background mass array not allocated"
            test_status = test_status + 1
        end if
        
        if (allocated(settings%background_settings%charge)) then
            if (size(settings%background_settings%charge) /= 3) then
                print *, "FAIL: Background charge array wrong size:", size(settings%background_settings%charge)
                test_status = test_status + 1
            else
                print *, "PASS: Background charge array size correct"
            end if
        else
            print *, "FAIL: Background charge array not allocated"
            test_status = test_status + 1
        end if
        
    end subroutine test_array_parameters_reading
    
    !> Create test namelist configuration file with array parameters
    subroutine create_test_array_config()
        integer :: unit
        integer :: iostat
        
        open(newunit=unit, file="test_array_config.conf", status="replace", action="write", iostat=iostat)
        if (iostat /= 0) then
            print *, "ERROR: Could not create test array config file, iostat =", iostat
            return
        end if
        
        write(unit, '(a)') "! KiLCA Test Configuration - Array parameters"
        write(unit, '(a)') ""
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 0.0"
        write(unit, '(a)') "  I0 = 1.0e13"
        write(unit, '(a)') "  flab = (1.0, 0.0)"
        write(unit, '(a)') "  dma = 5"
        write(unit, '(a)') "  flag_debug_ant = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.69"
        write(unit, '(a)') "  rp = 70.0"
        write(unit, '(a)') "  B0 = 23176.46"
        write(unit, '(a)') "  V_gal_sys = 1.0e9"
        write(unit, '(a)') "  m_i = 2.0"
        write(unit, '(a)') "  calc_back = 1"
        write(unit, '(a)') "  flag_debug_bg = 0"
        write(unit, '(a)') "  mass = 2.0, 1.0, 4.0"
        write(unit, '(a)') "  charge = 1.0, -1.0, 2.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 2"
        write(unit, '(a)') "  flag_emfield = 2"
        write(unit, '(a)') "  flag_additional = 2"
        write(unit, '(a)') "  flag_dispersion = 0"
        write(unit, '(a)') "  num_quants = 8"
        write(unit, '(a)') "  flag_quants = 1, 1, 1, 1, 1, 1, 1, 0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  rdim = 100"
        write(unit, '(a)') "  idim = 100"
        write(unit, '(a)') "  rfmin = 0.0"
        write(unit, '(a)') "  rfmax = 1.0e6"
        write(unit, '(a)') "  ifmin = -1.0e5"
        write(unit, '(a)') "  ifmax = 1.0e5"
        write(unit, '(a)') "/"
        
        close(unit)
        print *, "Test array config file created successfully"
    end subroutine create_test_array_config

end program test_namelist_array_params