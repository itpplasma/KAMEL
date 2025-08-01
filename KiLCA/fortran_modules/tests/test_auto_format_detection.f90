program test_auto_format_detection
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT, LEGACY_FORMAT
    use kilca_settings_m, only: settings_t, read_settings_auto_detect
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    character(len=*), parameter :: test_path = "./test_auto_format/"
    
    print *, "==========================================="
    print *, "Testing Automatic Format Detection [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_automatic_format_detection()
    
    if (test_status == 0) then
        print *, ""
        print *, "Automatic format detection tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Automatic format detection tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test automatic format detection functionality
    subroutine test_automatic_format_detection()
        print *, "Testing automatic format detection functionality..."
        
        ! Setup test directory
        call setup_test_directory()
        
        ! Test 1: Namelist format present - should detect and use namelist
        print *, "Test 1: Namelist format detection..."
        call create_namelist_settings()
        
        ! This should FAIL - auto format detection not implemented yet
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read settings with auto detection, error:", ierr
            test_status = test_status + 1
            return
        else
            print *, "SUCCESS: Auto detection read correctly"
        end if
        
        ! Check format used
        call validate_format_detection(1)  ! 1 = NAMELIST_FORMAT
        
        ! Test 2: Legacy format fallback - no namelist present
        print *, "Test 2: Legacy format fallback..."
        call cleanup_namelist_files()
        call create_legacy_settings()
        
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not fallback to legacy format, error:", ierr
            test_status = test_status + 1
            return
        else
            print *, "SUCCESS: Legacy fallback worked"
        end if
        
        ! Check format used
        call validate_format_detection(2)  ! 2 = LEGACY_FORMAT
        
        ! Test 3: No settings files present - should fail gracefully
        print *, "Test 3: No settings files present..."
        call cleanup_all_settings_files()
        
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Should have failed when no settings files present"
            test_status = test_status + 1
        else
            print *, "PASS: Correctly failed when no settings files present, error:", ierr
        end if
        
        ! Cleanup
        call cleanup_test_directory()
        
    end subroutine test_automatic_format_detection
    
    !> Setup test directory structure
    subroutine setup_test_directory()
        call system("mkdir -p " // test_path)
        print *, "Test directory created: ", test_path
    end subroutine setup_test_directory
    
    !> Create namelist format settings file
    subroutine create_namelist_settings()
        integer :: unit
        
        open(newunit=unit, file=test_path // "settings.conf", status="replace", action="write")
        
        write(unit, '(a)') "! Test namelist configuration"
        write(unit, '(a)') ""
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
        write(unit, '(a)') "  path2profiles = '../profiles/'"
        write(unit, '(a)') "  calc_back = 1"  
        write(unit, '(a)') "  flag_back = 'normal'"
        write(unit, '(a)') "  N = 5"
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
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "  flag_emfield = 1"
        write(unit, '(a)') "  flag_additional = 1"
        write(unit, '(a)') "  flag_dispersion = 0"
        write(unit, '(a)') "  flag_debug_out = 0"
        write(unit, '(a)') "  num_quants = 4"
        write(unit, '(a)') "  flag_quants = 1, 1, 0, 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  fname = 'auto_test.dat'"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  rdim = 100"
        write(unit, '(a)') "  rfmin = 0.0"
        write(unit, '(a)') "  rfmax = 1.0e6"
        write(unit, '(a)') "  idim = 50"
        write(unit, '(a)') "  ifmin = -1.0e5"
        write(unit, '(a)') "  ifmax = 1.0e5"
        write(unit, '(a)') "  stop_flag = 0"
        write(unit, '(a)') "  eps_res = 1.0e-6"
        write(unit, '(a)') "  eps_abs = 1.0e-8"
        write(unit, '(a)') "  eps_rel = 1.0e-6"
        write(unit, '(a)') "  delta = 1.0e-3"
        write(unit, '(a)') "  test_roots = 0"
        write(unit, '(a)') "  flag_debug_eig = 0"
        write(unit, '(a)') "  Nguess = 2"
        write(unit, '(a)') "  kmin = 1"
        write(unit, '(a)') "  kmax = 5"
        write(unit, '(a)') "  n_zeros = 5"
        write(unit, '(a)') "  use_winding = 0"
        write(unit, '(a)') "  fstart = (1.0e6, 0.0), (2.0e6, 0.0)"
        write(unit, '(a)') "/"
        
        close(unit)
        print *, "Namelist settings file created"
    end subroutine create_namelist_settings
    
    !> Create legacy format settings files
    subroutine create_legacy_settings()
        integer :: unit
        
        ! Create legacy antenna.in file
        open(newunit=unit, file=test_path // "antenna.in", status="replace", action="write")
        write(unit, '(a)') "90.0 5.0 1.0e12"  ! ra wa I0
        write(unit, '(a)') "1.0e6 0.0"         ! flab (real, imag)
        write(unit, '(a)') "2"                 ! dma
        write(unit, '(a)') "0 1"               ! flag_debug flag_eigmode
        write(unit, '(a)') "1 1 2 -1"          ! modes
        close(unit)
        
        ! Create legacy background.in file  
        open(newunit=unit, file=test_path // "background.in", status="replace", action="write")
        write(unit, '(a)') "170.0 65.0 25000.0"      ! rtor rp B0
        write(unit, '(a)') "../profiles/"             ! path2profiles
        write(unit, '(a)') "1"                        ! calc_back
        write(unit, '(a)') "normal"                   ! flag_back
        write(unit, '(a)') "5"                        ! N
        write(unit, '(a)') "1.0e9 1.0 2.0 1.0 1.0"   ! V_gal_sys V_scale m_i zele zion
        write(unit, '(a)') "0"                        ! flag_debug
        write(unit, '(a)') "2.0 1.0 4.0"             ! mass array
        write(unit, '(a)') "1.0 -1.0 2.0"            ! charge array
        close(unit)
        
        print *, "Legacy settings files created"
    end subroutine create_legacy_settings
    
    !> Validate format detection result
    subroutine validate_format_detection(expected_format)
        integer, intent(in) :: expected_format
        
        ! Now implemented - check format_used field
        if (settings%format_used /= expected_format) then
            print *, "FAIL: Wrong format detected, expected:", expected_format, "got:", settings%format_used
            test_status = test_status + 1
        else
            print *, "PASS: Format detection correct"
        end if
        
    end subroutine validate_format_detection
    
    !> Remove namelist files
    subroutine cleanup_namelist_files()
        call system("rm -f " // test_path // "settings.conf")
    end subroutine cleanup_namelist_files
    
    !> Remove all settings files
    subroutine cleanup_all_settings_files()
        call system("rm -f " // test_path // "settings.conf")
        call system("rm -f " // test_path // "antenna.in")
        call system("rm -f " // test_path // "background.in")
        call system("rm -f " // test_path // "output.in")
        call system("rm -f " // test_path // "eigenmode.in")
    end subroutine cleanup_all_settings_files
    
    !> Remove test directory
    subroutine cleanup_test_directory()
        call system("rm -rf " // test_path)
        print *, "Test directory cleaned up"
    end subroutine cleanup_test_directory
    

end program test_auto_format_detection