program test_settings_file_parsing
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/tmp/kilca_test/"
    
    ! Create test directory and files
    call setup_test_files()
    
    ! Test 1: Parse antenna settings
    call test_parse_antenna_settings()
    
    ! Test 2: Parse background settings
    call test_parse_background_settings()
    
    ! Test 3: Parse output settings
    call test_parse_output_settings()
    
    ! Test 4: Parse eigenmode settings
    call test_parse_eigenmode_settings()
    
    ! Test 5: Parse complete settings
    call test_parse_complete_settings()
    
    ! Test 6: Handle missing files
    call test_missing_files()
    
    ! Test 7: Handle malformed files
    call test_malformed_files()
    
    ! Clean up test files
    call cleanup_test_files()
    
    if (test_status == 0) then
        print *, "All file parsing tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine setup_test_files()
        integer :: unit, iostat
        
        ! Create test directory
        call system("mkdir -p " // trim(test_path))
        
        ! Create antenna.in file
        open(newunit=unit, file=trim(test_path) // "antenna.in", status='replace', iostat=iostat)
        if (iostat == 0) then
            write(unit, '(a)') "# Antenna settings:"
            write(unit, '(a)') "10.5  # ra: Small radius of antenna location (cm)"
            write(unit, '(a)') "1.2   # wa: Current density layer width (cm)"  
            write(unit, '(a)') "1000.0 # I0: Current in antenna coils (statamp)"
            write(unit, '(a)') "(5.0e7, 0.0) # flab: Frequency (Hz) in lab frame"
            write(unit, '(a)') "2     # dma: Dimension of modes array"
            write(unit, '(a)') "1     # flag_debug: Debug flag"
            write(unit, '(a)') "0     # flag_eigmode: Eigenmode search flag"
            write(unit, '(a)') "# End of antenna settings"
            close(unit)
        end if
        
        ! Create modes.in file
        open(newunit=unit, file=trim(test_path) // "modes.in", status='replace', iostat=iostat)
        if (iostat == 0) then
            write(unit, '(a)') "(1, 1)"
            write(unit, '(a)') "(2, 2)"
            close(unit)
        end if
        
        ! Create background.in file
        open(newunit=unit, file=trim(test_path) // "background.in", status='replace', iostat=iostat)
        if (iostat == 0) then
            write(unit, '(a)') "# Machine settings:"
            write(unit, '(a)') "625.0    # rtor: Big torus radius (cm)"
            write(unit, '(a)') "99.0     # rp: Plasma radius (cm)"
            write(unit, '(a)') "20000.0  # B0: Toroidal magnetic field (G)"
            write(unit, '(a)') "# Background field and plasma settings:"
            write(unit, '(a)') "./profiles/  # path2profiles: Path to profiles"
            write(unit, '(a)') "1        # calc_back: Calculate background flag"
            write(unit, '(a)') "normal   # flag_back: Background type"
            write(unit, '(a)') "5        # N: Splines degree"
            write(unit, '(a)') "0.0      # V_gal_sys: Velocity of moving frame"
            write(unit, '(a)') "1.0      # V_scale: Scale of Vz velocity profile"
            write(unit, '(a)') "1.0      # m_i: Ions mass in proton mass units"
            write(unit, '(a)') "1.0      # zele: Collision coefficient electrons"
            write(unit, '(a)') "1.0      # zion: Collision coefficient ions"
            write(unit, '(a)') "# Checking settings:"
            write(unit, '(a)') "0        # flag_debug: Debug flag"
            write(unit, '(a)') "# End of background settings"
            close(unit)
        end if
        
        ! Create output.in file
        open(newunit=unit, file=trim(test_path) // "output.in", status='replace', iostat=iostat)
        if (iostat == 0) then
            write(unit, '(a)') "# Output settings:"
            write(unit, '(a)') "1  # flag_background: Background output flag"
            write(unit, '(a)') "1  # flag_emfield: EM field output flag"
            write(unit, '(a)') "0  # flag_additional: Additional quantities flag"
            write(unit, '(a)') "0  # flag_dispersion: Dispersion output flag"
            write(unit, '(a)') "3  # num_quants: Number of quantities"
            write(unit, '(a)') "1  # flag_quant[0]: First quantity flag"
            write(unit, '(a)') "0  # flag_quant[1]: Second quantity flag"
            write(unit, '(a)') "1  # flag_quant[2]: Third quantity flag"
            write(unit, '(a)') "0  # flag_debug: Debug flag"
            close(unit)
        end if
        
        ! Create eigmode.in file
        open(newunit=unit, file=trim(test_path) // "eigmode.in", status='replace', iostat=iostat)
        if (iostat == 0) then
            write(unit, '(a)') "# Eigenmode settings:"
            write(unit, '(a)') "output.dat   # fname: Output file name"
            write(unit, '(a)') "1            # search_flag: Search option"
            write(unit, '(a)') "100          # rdim: Real dimension"
            write(unit, '(a)') "0.0          # rfmin: Real minimum"
            write(unit, '(a)') "1.0e9        # rfmax: Real maximum"
            write(unit, '(a)') "100          # idim: Imaginary dimension"
            write(unit, '(a)') "-1.0e6       # ifmin: Imaginary minimum"
            write(unit, '(a)') "1.0e6        # ifmax: Imaginary maximum"
            write(unit, '(a)') "0            # stop_flag: Stop flag"
            write(unit, '(a)') "1.0e-6       # eps_res: Residual accuracy"
            write(unit, '(a)') "1.0e-8       # eps_abs: Absolute accuracy"
            write(unit, '(a)') "1.0e-6       # eps_rel: Relative accuracy"
            write(unit, '(a)') "1.0e-6       # delta: Delta for derivatives"
            write(unit, '(a)') "0            # test_roots: Test roots flag"
            write(unit, '(a)') "0            # flag_debug: Debug flag"
            write(unit, '(a)') "2            # Nguess: Number of guess values"
            write(unit, '(a)') "1            # kmin: Minimum k value"
            write(unit, '(a)') "10           # kmax: Maximum k value"
            write(unit, '(a)') "(1.0e6, 0.0) # fstart[0]: First guess frequency"
            write(unit, '(a)') "(2.0e6, 0.0) # fstart[1]: Second guess frequency"
            write(unit, '(a)') "10           # n_zeros: Number of zeros to find"
            write(unit, '(a)') "0            # use_winding: Use winding number"
            close(unit)
        end if
        
    end subroutine setup_test_files
    
    subroutine test_parse_antenna_settings()
        type(antenna_t) :: ant
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing antenna settings parsing..."
        
        call antenna_read_settings_full(ant, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not parse antenna settings, error code:", ierr
            test_status = test_status + 1
            return
        end if
        
        ! Check parsed values
        if (abs(ant%ra - 10.5_dp) > 1.0e-10_dp) then
            print *, "FAIL: Antenna ra incorrect:", ant%ra
            test_status = test_status + 1
        end if
        
        if (abs(ant%wa - 1.2_dp) > 1.0e-10_dp) then
            print *, "FAIL: Antenna wa incorrect:", ant%wa
            test_status = test_status + 1
        end if
        
        if (abs(ant%I0 - 1000.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Antenna I0 incorrect:", ant%I0
            test_status = test_status + 1
        end if
        
        if (abs(real(ant%flab) - 5.0e7_dp) > 1.0e-3_dp) then
            print *, "FAIL: Antenna flab real part incorrect:", real(ant%flab)
            test_status = test_status + 1
        end if
        
        if (ant%dma /= 2) then
            print *, "FAIL: Antenna dma incorrect:", ant%dma
            test_status = test_status + 1
        end if
        
        if (ant%flag_debug /= 1) then
            print *, "FAIL: Antenna flag_debug incorrect:", ant%flag_debug
            test_status = test_status + 1
        end if
        
        if (.not. allocated(ant%modes)) then
            print *, "FAIL: Antenna modes not allocated"
            test_status = test_status + 1
        else if (size(ant%modes) /= 4) then
            print *, "FAIL: Antenna modes size incorrect:", size(ant%modes)
            test_status = test_status + 1
        else if (ant%modes(1) /= 1 .or. ant%modes(2) /= 1 .or. &
                 ant%modes(3) /= 2 .or. ant%modes(4) /= 2) then
            print *, "FAIL: Antenna modes values incorrect:", ant%modes
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(ant%modes)) deallocate(ant%modes)
        
        print *, "test_parse_antenna_settings completed"
    end subroutine test_parse_antenna_settings
    
    subroutine test_parse_background_settings()
        type(back_sett_t) :: bs
        integer :: ierr
        
        print *, "Testing background settings parsing..."
        
        call back_sett_read_settings_full(bs, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not parse background settings"
            test_status = test_status + 1
            return
        end if
        
        ! Check parsed values
        if (abs(bs%rtor - 625.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background rtor incorrect:", bs%rtor
            test_status = test_status + 1
        end if
        
        if (abs(bs%rp - 99.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background rp incorrect:", bs%rp
            test_status = test_status + 1
        end if
        
        if (abs(bs%B0 - 20000.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Background B0 incorrect:", bs%B0
            test_status = test_status + 1
        end if
        
        if (bs%calc_back /= 1) then
            print *, "FAIL: Background calc_back incorrect:", bs%calc_back
            test_status = test_status + 1
        end if
        
        if (allocated(bs%flag_back)) then
            if (trim(bs%flag_back) /= "normal") then
                print *, "FAIL: Background flag_back incorrect:", trim(bs%flag_back)
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Background flag_back not allocated"
            test_status = test_status + 1
        end if
        
        if (bs%N /= 5) then
            print *, "FAIL: Background N incorrect:", bs%N
            test_status = test_status + 1
        end if
        
        print *, "test_parse_background_settings completed"
    end subroutine test_parse_background_settings
    
    subroutine test_parse_output_settings()
        type(output_sett_t) :: os
        integer :: ierr
        
        print *, "Testing output settings parsing..."
        
        call output_sett_read_settings_full(os, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not parse output settings"
            test_status = test_status + 1
            return
        end if
        
        ! Check parsed values
        if (os%flag_background /= 1) then
            print *, "FAIL: Output flag_background incorrect:", os%flag_background
            test_status = test_status + 1
        end if
        
        if (os%flag_emfield /= 1) then
            print *, "FAIL: Output flag_emfield incorrect:", os%flag_emfield
            test_status = test_status + 1
        end if
        
        if (os%num_quants /= 3) then
            print *, "FAIL: Output num_quants incorrect:", os%num_quants
            test_status = test_status + 1
        end if
        
        if (allocated(os%flag_quants)) then
            if (size(os%flag_quants) /= 3) then
                print *, "FAIL: Output flag_quants size incorrect:", size(os%flag_quants)
                test_status = test_status + 1
            else if (os%flag_quants(1) /= 1 .or. os%flag_quants(2) /= 0 .or. &
                     os%flag_quants(3) /= 1) then
                print *, "FAIL: Output flag_quants values incorrect:", os%flag_quants
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Output flag_quants not allocated"
            test_status = test_status + 1
        end if
        
        print *, "test_parse_output_settings completed"
    end subroutine test_parse_output_settings
    
    subroutine test_parse_eigenmode_settings()
        type(eigmode_sett_t) :: es
        integer :: ierr
        
        print *, "Testing eigenmode settings parsing..."
        
        call eigmode_sett_read_settings_full(es, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not parse eigenmode settings"
            test_status = test_status + 1
            return
        end if
        
        ! Check parsed values
        if (es%search_flag /= 1) then
            print *, "FAIL: Eigenmode search_flag incorrect:", es%search_flag
            test_status = test_status + 1
        end if
        
        if (allocated(es%fname)) then
            if (trim(es%fname) /= "output.dat") then
                print *, "FAIL: Eigenmode fname incorrect:", trim(es%fname)
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Eigenmode fname not allocated"
            test_status = test_status + 1
        end if
        
        if (es%rdim /= 100) then
            print *, "FAIL: Eigenmode rdim incorrect:", es%rdim
            test_status = test_status + 1
        end if
        
        if (es%Nguess /= 2) then
            print *, "FAIL: Eigenmode Nguess incorrect:", es%Nguess
            test_status = test_status + 1
        end if
        
        if (allocated(es%fstart)) then
            if (size(es%fstart) /= 2) then
                print *, "FAIL: Eigenmode fstart size incorrect:", size(es%fstart)
                test_status = test_status + 1
            else if (abs(real(es%fstart(1)) - 1.0e6_dp) > 1.0e-3_dp .or. &
                     abs(real(es%fstart(2)) - 2.0e6_dp) > 1.0e-3_dp) then
                print *, "FAIL: Eigenmode fstart values incorrect:", es%fstart
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Eigenmode fstart not allocated"
            test_status = test_status + 1
        end if
        
        print *, "test_parse_eigenmode_settings completed"
    end subroutine test_parse_eigenmode_settings
    
    subroutine test_parse_complete_settings()
        type(settings_t), pointer :: sd
        integer :: ierr
        
        print *, "Testing complete settings parsing..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create settings"
            test_status = test_status + 1
            return
        end if
        
        call settings_read_all_full(sd, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not parse complete settings"
            test_status = test_status + 1
            call settings_destroy(sd, ierr)
            return
        end if
        
        ! Check that all subsettings were parsed
        if (abs(sd%antenna_settings%ra - 10.5_dp) > 1.0e-10_dp) then
            print *, "FAIL: Complete antenna settings incorrect"
            test_status = test_status + 1
        end if
        
        if (abs(sd%background_settings%rtor - 625.0_dp) > 1.0e-10_dp) then
            print *, "FAIL: Complete background settings incorrect"
            test_status = test_status + 1
        end if
        
        if (sd%output_settings%flag_background /= 1) then
            print *, "FAIL: Complete output settings incorrect"
            test_status = test_status + 1
        end if
        
        if (sd%eigmode_settings%search_flag /= 1) then
            print *, "FAIL: Complete eigenmode settings incorrect"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_parse_complete_settings completed"
    end subroutine test_parse_complete_settings
    
    subroutine test_missing_files()
        type(antenna_t) :: ant
        integer :: ierr
        
        print *, "Testing missing file handling..."
        
        call antenna_read_settings_full(ant, "/nonexistent/path/", ierr)
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Missing file not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected missing file"
        end if
        
        print *, "test_missing_files completed"
    end subroutine test_missing_files
    
    subroutine test_malformed_files()
        integer :: unit, iostat
        type(antenna_t) :: ant
        integer :: ierr
        character(len=256) :: bad_path = "/tmp/kilca_bad_test/"
        
        print *, "Testing malformed file handling..."
        
        ! Create bad test directory and malformed file
        call system("mkdir -p " // trim(bad_path))
        
        open(newunit=unit, file=trim(bad_path) // "antenna.in", status='replace', iostat=iostat)
        if (iostat == 0) then
            write(unit, '(a)') "# Bad antenna file"
            write(unit, '(a)') "not_a_number  # ra: This should fail"
            close(unit)
        end if
        
        call antenna_read_settings_full(ant, bad_path, ierr)
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Malformed file not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected malformed file"
        end if
        
        ! Clean up
        call system("rm -rf " // trim(bad_path))
        
        print *, "test_malformed_files completed"
    end subroutine test_malformed_files
    
    subroutine cleanup_test_files()
        call system("rm -rf " // trim(test_path))
    end subroutine cleanup_test_files

end program test_settings_file_parsing