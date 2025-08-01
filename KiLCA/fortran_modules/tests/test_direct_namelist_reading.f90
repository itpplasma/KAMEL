program test_direct_namelist_reading
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: settings_t, settings_read_namelist
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Direct Namelist Reading [RED PHASE - SHOULD FAIL]"
    print *, "=========================================="
    
    call test_basic_namelist_reading()
    
    if (test_status == 0) then
        print *, ""
        print *, "Direct namelist reading tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Direct namelist reading tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test basic namelist reading functionality
    subroutine test_basic_namelist_reading()
        print *, "Testing basic namelist reading functionality..."
        
        ! Create test configuration file first
        call create_test_namelist_config()
        
        ! This should now PASS - procedure implemented
        print *, "Attempting to read namelist settings..."
        call settings_read_namelist("test_kilca.conf", settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read namelist settings, error:", ierr
            test_status = test_status + 1
            return
        else
            print *, "SUCCESS: Namelist settings read correctly"
        end if
        
        ! Verify all basic parameters read correctly
        if (abs(settings%antenna_settings%ra - 90.0_dp) > 1.0e-12_dp) then
            print *, "FAIL: Antenna ra parameter not read correctly:", settings%antenna_settings%ra
            test_status = test_status + 1
        else
            print *, "PASS: Antenna ra parameter read correctly"
        end if
        
        if (abs(settings%background_settings%rtor - 170.69_dp) > 1.0e-12_dp) then
            print *, "FAIL: Background rtor parameter not read correctly:", settings%background_settings%rtor
            test_status = test_status + 1
        else
            print *, "PASS: Background rtor parameter read correctly"
        end if
        
        if (settings%output_settings%flag_background /= 2) then
            print *, "FAIL: Output flag_background not read correctly:", settings%output_settings%flag_background
            test_status = test_status + 1
        else
            print *, "PASS: Output flag_background read correctly"
        end if
        
        if (settings%eigmode_settings%search_flag /= 1) then
            print *, "FAIL: Eigenmode search_flag not read correctly:", settings%eigmode_settings%search_flag
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode search_flag read correctly"
        end if
        
    end subroutine test_basic_namelist_reading
    
    !> Create test namelist configuration file
    subroutine create_test_namelist_config()
        integer :: unit
        integer :: iostat
        
        open(newunit=unit, file="test_kilca.conf", status="replace", action="write", iostat=iostat)
        if (iostat /= 0) then
            print *, "ERROR: Could not create test config file, iostat =", iostat
            return
        end if
        
        write(unit, '(a)') "! KiLCA Test Configuration - Direct namelist format"
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
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 2"
        write(unit, '(a)') "  flag_emfield = 2"
        write(unit, '(a)') "  flag_additional = 2"
        write(unit, '(a)') "  flag_dispersion = 0"
        write(unit, '(a)') "  num_quants = 8"
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
        print *, "Test config file created successfully"
    end subroutine create_test_namelist_config

end program test_direct_namelist_reading