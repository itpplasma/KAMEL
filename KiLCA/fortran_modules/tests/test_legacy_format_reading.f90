program test_legacy_format_reading
    use kilca_types_m, only: dp, KILCA_SUCCESS, LEGACY_FORMAT
    use kilca_settings_m, only: settings_t, read_settings_auto_detect
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    character(len=*), parameter :: test_path = "./test_legacy_format/"
    
    print *, "==========================================="
    print *, "Testing Legacy Format Reading [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_legacy_format_functionality()
    
    if (test_status == 0) then
        print *, ""
        print *, "Legacy format reading tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Legacy format reading tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test legacy format reading functionality
    subroutine test_legacy_format_functionality()
        print *, "Testing legacy format reading functionality..."
        
        ! Setup test directory
        call setup_test_directory()
        
        ! Test 1: Read complete legacy format
        print *, ""
        print *, "Test 1: Complete legacy format reading..."
        call create_complete_legacy_files()
        
        ! This should work since auto-detect is implemented
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read legacy format, error:", ierr
            test_status = test_status + 1
        else
            print *, "SUCCESS: Legacy format read without error"
            
            ! But actual values should be wrong since implementation is incomplete
            call validate_legacy_antenna_values()
            call validate_legacy_background_values()
            call validate_legacy_output_values()
            call validate_legacy_eigenmode_values()
        end if
        
        ! Test 2: Partial legacy files (missing files)
        print *, ""
        print *, "Test 2: Partial legacy files handling..."
        call cleanup_test_files()
        call create_partial_legacy_files()
        
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Should have failed with partial legacy files"
            test_status = test_status + 1
        else
            print *, "PASS: Correctly failed with partial legacy files, error:", ierr
        end if
        
        ! Test 3: Corrupt legacy file format
        print *, ""
        print *, "Test 3: Corrupt legacy file handling..."
        call cleanup_test_files()
        call create_corrupt_legacy_files()
        
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Should have failed with corrupt legacy files"
            test_status = test_status + 1
        else
            print *, "PASS: Correctly failed with corrupt legacy files, error:", ierr
        end if
        
        ! Test 4: Legacy format with extra parameters
        print *, ""
        print *, "Test 4: Legacy format with extra parameters..."
        call cleanup_test_files()
        call create_extended_legacy_files()
        
        call read_settings_auto_detect(test_path, settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Should handle extra parameters gracefully, error:", ierr
            test_status = test_status + 1
        else
            print *, "SUCCESS: Handled extra parameters gracefully"
        end if
        
        ! Cleanup
        call cleanup_test_directory()
        
    end subroutine test_legacy_format_functionality
    
    !> Validate antenna values from legacy format
    subroutine validate_legacy_antenna_values()
        real(dp) :: expected_ra = 85.5_dp
        real(dp) :: expected_wa = 3.0_dp
        real(dp) :: expected_I0 = 2.5e12_dp
        
        if (abs(settings%antenna_settings%ra - expected_ra) > 1.0e-10_dp) then
            print *, "FAIL: Antenna ra incorrect:", settings%antenna_settings%ra, "expected:", expected_ra
            test_status = test_status + 1
        else
            print *, "PASS: Antenna ra correct"
        end if
        
        if (abs(settings%antenna_settings%wa - expected_wa) > 1.0e-10_dp) then
            print *, "FAIL: Antenna wa incorrect:", settings%antenna_settings%wa, "expected:", expected_wa
            test_status = test_status + 1
        else
            print *, "PASS: Antenna wa correct"
        end if
        
        if (abs(settings%antenna_settings%I0 - expected_I0) > 1.0e-10_dp) then
            print *, "FAIL: Antenna I0 incorrect:", settings%antenna_settings%I0, "expected:", expected_I0
            test_status = test_status + 1
        else
            print *, "PASS: Antenna I0 correct"
        end if
        
        ! Complex number check
        if (abs(real(settings%antenna_settings%flab) - 2.0e6_dp) > 1.0e-10_dp .or. &
            abs(aimag(settings%antenna_settings%flab) - 0.5e6_dp) > 1.0e-10_dp) then
            print *, "FAIL: Antenna flab incorrect:", settings%antenna_settings%flab
            test_status = test_status + 1
        else
            print *, "PASS: Antenna flab correct"
        end if
        
        ! Array check
        if (.not. allocated(settings%antenna_settings%modes)) then
            print *, "FAIL: Antenna modes not allocated"
            test_status = test_status + 1
        else if (size(settings%antenna_settings%modes) /= 6) then
            print *, "FAIL: Antenna modes wrong size:", size(settings%antenna_settings%modes)
            test_status = test_status + 1
        else if (settings%antenna_settings%modes(1) /= 1 .or. settings%antenna_settings%modes(6) /= -2) then
            print *, "FAIL: Antenna modes values incorrect"
            test_status = test_status + 1
        else
            print *, "PASS: Antenna modes correct"
        end if
        
    end subroutine validate_legacy_antenna_values
    
    !> Validate background values from legacy format
    subroutine validate_legacy_background_values()
        if (abs(settings%background_settings%rtor - 165.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background rtor incorrect:", settings%background_settings%rtor
            test_status = test_status + 1
        else
            print *, "PASS: Background rtor correct"
        end if
        
        if (abs(settings%background_settings%rp - 55.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background rp incorrect:", settings%background_settings%rp
            test_status = test_status + 1
        else
            print *, "PASS: Background rp correct"
        end if
        
        if (abs(settings%background_settings%B0 - 28000.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background B0 incorrect:", settings%background_settings%B0
            test_status = test_status + 1
        else
            print *, "PASS: Background B0 correct"
        end if
        
        ! String check
        if (.not. allocated(settings%background_settings%path2profiles)) then
            print *, "FAIL: Background path2profiles not allocated"
            test_status = test_status + 1
        else if (settings%background_settings%path2profiles /= "./profiles/") then
            print *, "FAIL: Background path2profiles incorrect:", settings%background_settings%path2profiles
            test_status = test_status + 1
        else
            print *, "PASS: Background path2profiles correct"
        end if
        
        ! Arrays check
        if (.not. allocated(settings%background_settings%mass)) then
            print *, "FAIL: Background mass array not allocated"
            test_status = test_status + 1
        else if (abs(settings%background_settings%mass(1) - 2.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background mass values incorrect"
            test_status = test_status + 1
        else
            print *, "PASS: Background mass array correct"
        end if
        
    end subroutine validate_legacy_background_values
    
    !> Validate output values from legacy format
    subroutine validate_legacy_output_values()
        if (settings%output_settings%flag_background /= 1) then
            print *, "FAIL: Output flag_background incorrect:", settings%output_settings%flag_background
            test_status = test_status + 1
        else
            print *, "PASS: Output flag_background correct"
        end if
        
        if (settings%output_settings%num_quants /= 5) then
            print *, "FAIL: Output num_quants incorrect:", settings%output_settings%num_quants
            test_status = test_status + 1
        else
            print *, "PASS: Output num_quants correct"
        end if
        
        ! Array check
        if (.not. allocated(settings%output_settings%flag_quants)) then
            print *, "FAIL: Output flag_quants not allocated"
            test_status = test_status + 1
        else if (size(settings%output_settings%flag_quants) /= 5) then
            print *, "FAIL: Output flag_quants wrong size"
            test_status = test_status + 1
        else
            print *, "PASS: Output flag_quants correct"
        end if
        
    end subroutine validate_legacy_output_values
    
    !> Validate eigenmode values from legacy format
    subroutine validate_legacy_eigenmode_values()
        if (.not. allocated(settings%eigmode_settings%fname)) then
            print *, "FAIL: Eigenmode fname not allocated"
            test_status = test_status + 1
        else if (settings%eigmode_settings%fname /= "legacy_output.dat") then
            print *, "FAIL: Eigenmode fname incorrect:", settings%eigmode_settings%fname
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode fname correct"
        end if
        
        if (settings%eigmode_settings%search_flag /= 1) then
            print *, "FAIL: Eigenmode search_flag incorrect:", settings%eigmode_settings%search_flag
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode search_flag correct"
        end if
        
        if (abs(settings%eigmode_settings%eps_res - 1.0e-7_dp) > 1.0e-15_dp) then
            print *, "FAIL: Eigenmode eps_res incorrect:", settings%eigmode_settings%eps_res
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode eps_res correct"
        end if
        
        ! Complex array check
        if (.not. allocated(settings%eigmode_settings%fstart)) then
            print *, "FAIL: Eigenmode fstart not allocated"
            test_status = test_status + 1
        else if (size(settings%eigmode_settings%fstart) /= 3) then
            print *, "FAIL: Eigenmode fstart wrong size:", size(settings%eigmode_settings%fstart)
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode fstart allocated correctly"
        end if
        
    end subroutine validate_legacy_eigenmode_values
    
    !> Setup test directory
    subroutine setup_test_directory()
        call system("mkdir -p " // test_path)
        print *, "Test directory created: ", test_path
    end subroutine setup_test_directory
    
    !> Create complete legacy format files
    subroutine create_complete_legacy_files()
        integer :: unit
        
        ! Create antenna.in
        open(newunit=unit, file=test_path // "antenna.in", status="replace", action="write")
        write(unit, '(3(e15.7,1x))') 85.5_dp, 3.0_dp, 2.5e12_dp      ! ra wa I0
        write(unit, '(2(e15.7,1x))') 2.0e6_dp, 0.5e6_dp               ! flab real imag
        write(unit, '(i5)') 3                                          ! dma
        write(unit, '(2(i5,1x))') 0, 1                                 ! flag_debug flag_eigmode
        write(unit, '(6(i5,1x))') 1, 1, 2, -1, 3, -2                  ! modes
        close(unit)
        
        ! Create background.in
        open(newunit=unit, file=test_path // "background.in", status="replace", action="write")
        write(unit, '(3(e15.7,1x))') 165.0_dp, 55.0_dp, 28000.0_dp    ! rtor rp B0
        write(unit, '(a)') "./profiles/"                               ! path2profiles
        write(unit, '(i5)') 1                                          ! calc_back
        write(unit, '(a)') "experimental"                              ! flag_back
        write(unit, '(i5)') 6                                          ! N
        write(unit, '(5(e15.7,1x))') 1.5e9_dp, 0.9_dp, 2.5_dp, 1.0_dp, 1.0_dp  ! V_gal_sys V_scale m_i zele zion
        write(unit, '(i5)') 0                                          ! flag_debug
        write(unit, '(3(e15.7,1x))') 2.0_dp, 1.0_dp, 4.0_dp          ! mass
        write(unit, '(3(e15.7,1x))') 1.0_dp, -1.0_dp, 2.0_dp         ! charge
        close(unit)
        
        ! Create output.in
        open(newunit=unit, file=test_path // "output.in", status="replace", action="write")
        write(unit, '(4(i5,1x))') 1, 2, 1, 0                          ! flags
        write(unit, '(i5)') 0                                          ! flag_debug
        write(unit, '(i5)') 5                                          ! num_quants
        write(unit, '(5(i5,1x))') 1, 1, 0, 1, 1                       ! flag_quants
        close(unit)
        
        ! Create eigenmode.in
        open(newunit=unit, file=test_path // "eigenmode.in", status="replace", action="write")
        write(unit, '(a)') "legacy_output.dat"                         ! fname
        write(unit, '(i5)') 1                                          ! search_flag
        write(unit, '(2(i5,1x))') 120, 60                             ! rdim idim
        write(unit, '(2(e15.7,1x))') 0.0_dp, 3.0e6_dp                 ! rfmin rfmax
        write(unit, '(2(e15.7,1x))') -2.0e5_dp, 2.0e5_dp             ! ifmin ifmax
        write(unit, '(i5)') 0                                          ! stop_flag
        write(unit, '(4(e15.7,1x))') 1.0e-7_dp, 1.0e-9_dp, 1.0e-7_dp, 1.0e-4_dp  ! eps_res eps_abs eps_rel delta
        write(unit, '(2(i5,1x))') 0, 0                                 ! test_roots flag_debug
        write(unit, '(4(i5,1x))') 3, 1, 6, 8                          ! Nguess kmin kmax n_zeros
        write(unit, '(i5)') 0                                          ! use_winding
        ! fstart array (3 complex numbers)
        write(unit, '(2(e15.7,1x))') 1.0e6_dp, 0.0_dp
        write(unit, '(2(e15.7,1x))') 1.5e6_dp, 0.1e6_dp
        write(unit, '(2(e15.7,1x))') 2.0e6_dp, -0.05e6_dp
        close(unit)
        
        print *, "Complete legacy files created"
    end subroutine create_complete_legacy_files
    
    !> Create partial legacy files (missing some required files)
    subroutine create_partial_legacy_files()
        integer :: unit
        
        ! Only create antenna.in
        open(newunit=unit, file=test_path // "antenna.in", status="replace", action="write")
        write(unit, '(3(e15.7,1x))') 85.5_dp, 3.0_dp, 2.5e12_dp
        write(unit, '(2(e15.7,1x))') 2.0e6_dp, 0.5e6_dp
        write(unit, '(i5)') 3
        write(unit, '(2(i5,1x))') 0, 1
        write(unit, '(6(i5,1x))') 1, 1, 2, -1, 3, -2
        close(unit)
        
        ! Missing background.in!
        
        print *, "Partial legacy files created (missing background.in)"
    end subroutine create_partial_legacy_files
    
    !> Create corrupt legacy files
    subroutine create_corrupt_legacy_files()
        integer :: unit
        
        ! Create antenna.in with bad format
        open(newunit=unit, file=test_path // "antenna.in", status="replace", action="write")
        write(unit, '(a)') "This is not valid numeric data"
        write(unit, '(a)') "85.5 three point zero"  ! Text instead of numbers
        close(unit)
        
        ! Create background.in with missing data
        open(newunit=unit, file=test_path // "background.in", status="replace", action="write")
        write(unit, '(2(e15.7,1x))') 165.0_dp, 55.0_dp  ! Missing B0
        write(unit, '(a)') "./profiles/"
        ! Missing rest of data
        close(unit)
        
        print *, "Corrupt legacy files created"
    end subroutine create_corrupt_legacy_files
    
    !> Create extended legacy files with extra parameters
    subroutine create_extended_legacy_files()
        integer :: unit
        
        ! Create antenna.in with extra lines
        open(newunit=unit, file=test_path // "antenna.in", status="replace", action="write")
        write(unit, '(3(e15.7,1x))') 85.5_dp, 3.0_dp, 2.5e12_dp
        write(unit, '(2(e15.7,1x))') 2.0e6_dp, 0.5e6_dp
        write(unit, '(i5)') 3
        write(unit, '(2(i5,1x))') 0, 1
        write(unit, '(6(i5,1x))') 1, 1, 2, -1, 3, -2
        write(unit, '(a)') "# Extra comment line"
        write(unit, '(e15.7)') 999.0_dp  ! Extra parameter
        close(unit)
        
        ! Create background.in with extra data
        open(newunit=unit, file=test_path // "background.in", status="replace", action="write")
        write(unit, '(3(e15.7,1x))') 165.0_dp, 55.0_dp, 28000.0_dp
        write(unit, '(a)') "./profiles/"
        write(unit, '(i5)') 1
        write(unit, '(a)') "experimental"
        write(unit, '(i5)') 6
        write(unit, '(5(e15.7,1x))') 1.5e9_dp, 0.9_dp, 2.5_dp, 1.0_dp, 1.0_dp
        write(unit, '(i5)') 0
        write(unit, '(3(e15.7,1x))') 2.0_dp, 1.0_dp, 4.0_dp
        write(unit, '(3(e15.7,1x))') 1.0_dp, -1.0_dp, 2.0_dp
        write(unit, '(a)') "Extra data that should be ignored"
        close(unit)
        
        ! Create basic output.in
        open(newunit=unit, file=test_path // "output.in", status="replace", action="write")
        write(unit, '(4(i5,1x))') 1, 2, 1, 0
        write(unit, '(i5)') 0
        write(unit, '(i5)') 5
        write(unit, '(5(i5,1x))') 1, 1, 0, 1, 1
        close(unit)
        
        ! Create basic eigenmode.in
        open(newunit=unit, file=test_path // "eigenmode.in", status="replace", action="write")
        write(unit, '(a)') "legacy_output.dat"
        write(unit, '(i5)') 1
        write(unit, '(2(i5,1x))') 120, 60
        write(unit, '(2(e15.7,1x))') 0.0_dp, 3.0e6_dp
        write(unit, '(2(e15.7,1x))') -2.0e5_dp, 2.0e5_dp
        write(unit, '(i5)') 0
        write(unit, '(4(e15.7,1x))') 1.0e-7_dp, 1.0e-9_dp, 1.0e-7_dp, 1.0e-4_dp
        write(unit, '(2(i5,1x))') 0, 0
        write(unit, '(4(i5,1x))') 3, 1, 6, 8
        write(unit, '(i5)') 0
        write(unit, '(2(e15.7,1x))') 1.0e6_dp, 0.0_dp
        write(unit, '(2(e15.7,1x))') 1.5e6_dp, 0.1e6_dp
        write(unit, '(2(e15.7,1x))') 2.0e6_dp, -0.05e6_dp
        close(unit)
        
        print *, "Extended legacy files created"
    end subroutine create_extended_legacy_files
    
    !> Clean up test files
    subroutine cleanup_test_files()
        call system("rm -f " // test_path // "*.in")
        call system("rm -f " // test_path // "*.conf")
    end subroutine cleanup_test_files
    
    !> Clean up test directory
    subroutine cleanup_test_directory()
        call system("rm -rf " // test_path)
        print *, "Test directory cleaned up"
    end subroutine cleanup_test_directory

end program test_legacy_format_reading