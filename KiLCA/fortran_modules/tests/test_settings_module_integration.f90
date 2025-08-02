program test_settings_module_integration
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT
    use kilca_settings_m, only: settings_t, settings_create, &
                                back_sett_read_settings, &
                                antenna_read_settings, &
                                output_read_settings, &
                                eigmode_read_settings, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(settings_t), pointer :: sd => null()
    character(len=*), parameter :: test_path = "./test_project/"
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Settings Module Integration [RED PHASE - SHOULD FAIL]"
    print *, "=========================================="
    
    call test_settings_integration()
    
    if (test_status == 0) then
        print *, ""
        print *, "Settings integration tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Settings integration tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test integration with existing settings procedures
    subroutine test_settings_integration()
        print *, "Testing settings module integration..."
        
        ! Setup test directory and files
        call setup_test_environment()
        
        ! Test 1: Create settings with namelist backend
        print *, ""
        print *, "Test 1: Create settings with namelist backend..."
        call settings_create(sd, test_path, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create settings, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: Settings created successfully"
        end if
        
        ! Test 2: Verify backward compatibility with back_sett_read_settings
        print *, ""
        print *, "Test 2: Backward compatibility with back_sett_read_settings..."
        call back_sett_read_settings(sd%background_settings, test_path, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_read_settings failed, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: back_sett_read_settings succeeded"
            
            ! Verify it uses namelist internally
            if (sd%background_settings%format_used == NAMELIST_FORMAT) then
                print *, "PASS: Using namelist format internally"
            else
                print *, "FAIL: Not using namelist format internally"
                test_status = test_status + 1
            end if
        end if
        
        ! Test 3: Verify antenna_read_settings uses namelist
        print *, ""
        print *, "Test 3: Antenna settings with namelist backend..."
        call antenna_read_settings(sd%antenna_settings, test_path, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_read_settings failed, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: antenna_read_settings succeeded"
        end if
        
        ! Test 4: Integration flag to force namelist usage
        print *, ""
        print *, "Test 4: Force namelist backend globally..."
        call settings_integrate_namelist_backend(.true.)
        
        ! All subsequent reads should use namelist
        call output_read_settings(sd%output_settings, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: output_read_settings with namelist backend failed"
            test_status = test_status + 1
        else
            print *, "PASS: output_read_settings with namelist backend succeeded"
        end if
        
        ! Test 5: Verify settings are identical between formats
        print *, ""
        print *, "Test 5: Verify identical results between formats..."
        call test_format_equivalence()
        
        ! Test 6: Performance comparison
        print *, ""
        print *, "Test 6: Performance comparison..."
        call test_performance_improvement()
        
        ! Clean up
        call cleanup_test_environment()
        
    end subroutine test_settings_integration
    
    !> Test that namelist and legacy formats produce identical results
    subroutine test_format_equivalence()
        type(settings_t), pointer :: sd_namelist => null(), sd_legacy => null()
        
        ! Read with namelist format
        call settings_integrate_namelist_backend(.true.)
        call settings_create(sd_namelist, test_path, ierr)
        
        ! Read with legacy format
        call settings_integrate_namelist_backend(.false.)
        call settings_create(sd_legacy, test_path, ierr)
        
        ! Compare all values
        if (abs(sd_namelist%antenna_settings%ra - sd_legacy%antenna_settings%ra) < 1.0e-12_dp) then
            print *, "PASS: Antenna settings match between formats"
        else
            print *, "FAIL: Antenna settings differ between formats"
            test_status = test_status + 1
        end if
        
        if (abs(sd_namelist%background_settings%rtor - sd_legacy%background_settings%rtor) < 1.0e-12_dp) then
            print *, "PASS: Background settings match between formats"
        else
            print *, "FAIL: Background settings differ between formats"
            test_status = test_status + 1
        end if
        
    end subroutine test_format_equivalence
    
    !> Test performance improvement with namelist format
    subroutine test_performance_improvement()
        real :: start_time, end_time
        real :: namelist_time, legacy_time
        integer :: i
        type(settings_t), pointer :: sd_temp => null()
        
        ! Measure namelist performance
        call settings_integrate_namelist_backend(.true.)
        call cpu_time(start_time)
        do i = 1, 100
            call settings_create(sd_temp, test_path, ierr)
        end do
        call cpu_time(end_time)
        namelist_time = end_time - start_time
        
        ! Measure legacy performance
        call settings_integrate_namelist_backend(.false.)
        call cpu_time(start_time)
        do i = 1, 100
            call settings_create(sd_temp, test_path, ierr)
        end do
        call cpu_time(end_time)
        legacy_time = end_time - start_time
        
        print *, "Namelist time:", namelist_time, "seconds"
        print *, "Legacy time:", legacy_time, "seconds"
        
        if (namelist_time <= legacy_time) then
            print *, "PASS: Namelist format is not slower"
        else
            print *, "INFO: Namelist format is slower by", &
                     100.0 * (namelist_time - legacy_time) / legacy_time, "%"
        end if
        
    end subroutine test_performance_improvement
    
    !> Setup test environment
    subroutine setup_test_environment()
        integer :: unit
        
        call system("mkdir -p " // test_path)
        
        ! Create settings.conf
        open(newunit=unit, file=test_path // "settings.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "  path2profiles = './profiles/'"
        write(unit, '(a)') "  calc_back = 1"
        write(unit, '(a)') "  flag_back = 'experimental'"
        write(unit, '(a)') "  N = 3"
        write(unit, '(a)') "  V_gal_sys = 1.5e9"
        write(unit, '(a)') "  V_scale = 0.9"
        write(unit, '(a)') "  m_i = 2.5"
        write(unit, '(a)') "  zele = 1.0"
        write(unit, '(a)') "  zion = 1.0"
        write(unit, '(a)') "  flag_debug_bg = 0"
        write(unit, '(a)') "  mass = 2.0, 1.0, 4.0"
        write(unit, '(a)') "  charge = 1.0, -1.0, 2.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "  flag_emfield = 1"
        write(unit, '(a)') "  flag_additional = 0"
        write(unit, '(a)') "  flag_dispersion = 1"
        write(unit, '(a)') "  flag_debug_out = 0"
        write(unit, '(a)') "  num_quants = 5"
        write(unit, '(a)') "  flag_quants = 1, 1, 0, 1, 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  fname = 'eigenmode_output.dat'"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  rdim = 120"
        write(unit, '(a)') "  idim = 60"
        write(unit, '(a)') "  rfmin = 0.0"
        write(unit, '(a)') "  rfmax = 3.0e6"
        write(unit, '(a)') "  ifmin = -2.0e5"
        write(unit, '(a)') "  ifmax = 2.0e5"
        write(unit, '(a)') "  stop_flag = 0"
        write(unit, '(a)') "  eps_res = 1.0e-7"
        write(unit, '(a)') "  eps_abs = 1.0e-9"
        write(unit, '(a)') "  eps_rel = 1.0e-7"
        write(unit, '(a)') "  delta = 1.0e-4"
        write(unit, '(a)') "  test_roots = 0"
        write(unit, '(a)') "  flag_debug_eig = 0"
        write(unit, '(a)') "  Nguess = 3"
        write(unit, '(a)') "  kmin = 1"
        write(unit, '(a)') "  kmax = 6"
        write(unit, '(a)') "  n_zeros = 8"
        write(unit, '(a)') "  use_winding = 0"
        write(unit, '(a)') "  fstart = (1.0e6, 0.0), (1.5e6, 0.1e6), (2.0e6, -0.05e6)"
        write(unit, '(a)') "/"
        close(unit)
        
        ! Also create legacy format files for comparison
        call create_legacy_format_files()
        
        print *, "Test environment created"
        
    end subroutine setup_test_environment
    
    !> Create legacy format files
    subroutine create_legacy_format_files()
        integer :: unit
        
        ! antenna.in
        open(newunit=unit, file=test_path // "antenna.in", status="replace")
        write(unit, '(3(e15.7,1x))') 90.0_dp, 5.0_dp, 1.0e12_dp
        write(unit, '(2(e15.7,1x))') 1.0e6_dp, 0.0_dp
        write(unit, '(i5)') 2
        write(unit, '(2(i5,1x))') 0, 1
        write(unit, '(4(i5,1x))') 1, 1, 2, -1
        close(unit)
        
        ! Similar for other .in files...
        print *, "Legacy format files created"
        
    end subroutine create_legacy_format_files
    
    !> Clean up test environment
    subroutine cleanup_test_environment()
        call system("rm -rf " // test_path)
        print *, "Test environment cleaned up"
    end subroutine cleanup_test_environment

end program test_settings_module_integration