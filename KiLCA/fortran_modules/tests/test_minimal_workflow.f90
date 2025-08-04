program test_minimal_workflow
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT, LEGACY_FORMAT
    use kilca_settings_m, only: settings_t, settings_initialize_defaults, &
                                settings_read_all, settings_destroy, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(settings_t), pointer :: sd_legacy => null()
    integer :: ierr
    character(len=*), parameter :: test_path = "./test_minimal_workflow/"
    
    print *, "Creating minimal test environment..."
    call setup_minimal_test()
    
    print *, "Step 1: Disable namelist backend"
    call settings_integrate_namelist_backend(.false.)
    
    print *, "Step 2: Create settings object"
    call settings_initialize_defaults(sd_legacy, test_path, ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "FAILED to create settings:", ierr
        stop 1
    end if
    print *, "Settings created successfully"
    
    print *, "Step 3: Read settings with legacy backend"
    call settings_read_all(sd_legacy, ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "FAILED to read settings:", ierr
        print *, "Expected this to work with legacy files"
    else
        print *, "SUCCESS: Read settings with legacy backend"
    end if
    
    print *, "Cleaning up..."
    call settings_destroy(sd_legacy, ierr)
    call system("rm -rf " // test_path)
    
contains

    subroutine setup_minimal_test()
        integer :: unit
        
        call system("mkdir -p " // test_path)
        
        ! Create ONLY legacy format files (no settings.conf)
        
        ! antenna.in
        open(newunit=unit, file=test_path // "antenna.in", status="replace")
        write(unit, '(3(e15.7,1x))') 90.0_dp, 5.0_dp, 1.0e12_dp
        write(unit, '(2(e15.7,1x))') 1.0e6_dp, 0.0_dp
        write(unit, '(i5)') 2
        write(unit, '(2(i5,1x))') 0, 1
        write(unit, '(4(i5,1x))') 1, 1, 2, -1
        close(unit)
        
        ! background.in (corrected order)
        open(newunit=unit, file=test_path // "background.in", status="replace")
        write(unit, '(3(e15.7,1x))') 170.0_dp, 65.0_dp, 25000.0_dp  ! rtor, rp, B0
        write(unit, '(a)') "./profiles/"                              ! path2profiles
        write(unit, '(i5)') 1                                        ! calc_back
        write(unit, '(a)') "experimental"                            ! flag_back
        write(unit, '(i5)') 3                                        ! N
        write(unit, '(2(e15.7,1x))') 1.5e9_dp, 0.9_dp               ! V_gal_sys, V_scale
        write(unit, '(3(e15.7,1x))') 2.5_dp, 1.0_dp, 1.0_dp         ! m_i, zele, zion
        write(unit, '(i5)') 0                                        ! flag_debug_bg
        write(unit, '(3(e15.7,1x))') 2.0_dp, 1.0_dp, 4.0_dp         ! mass
        write(unit, '(3(e15.7,1x))') 1.0_dp, -1.0_dp, 2.0_dp        ! charge
        close(unit)
        
        ! output.in (corrected format)
        open(newunit=unit, file=test_path // "output.in", status="replace")
        write(unit, '(4(i5,1x))') 1, 1, 0, 1        ! flag_background, flag_emfield, flag_additional, flag_dispersion
        write(unit, '(i5)') 0                       ! flag_debug
        write(unit, '(i5)') 5                       ! num_quants
        write(unit, '(5(i5,1x))') 1, 1, 0, 1, 1     ! flag_quants array
        close(unit)
        
        ! eigenmode.in (corrected format)
        open(newunit=unit, file=test_path // "eigenmode.in", status="replace")
        write(unit, '(a)') "eigenmode.dat"                           ! fname
        write(unit, '(i5)') 1                                        ! search_flag
        write(unit, '(2(i5,1x))') 120, 60                            ! rdim, idim
        write(unit, '(2(e15.7,1x))') 0.0_dp, 3.0e6_dp               ! rfmin, rfmax
        write(unit, '(2(e15.7,1x))') -2.0e5_dp, 2.0e5_dp            ! ifmin, ifmax
        write(unit, '(i5)') 0                                        ! stop_flag
        write(unit, '(4(e15.7,1x))') 1.0e-7_dp, 1.0e-9_dp, 1.0e-7_dp, 1.0e-4_dp  ! eps_res, eps_abs, eps_rel, delta
        write(unit, '(2(i5,1x))') 0, 0                               ! test_roots, flag_debug_eig
        write(unit, '(4(i5,1x))') 3, 1, 6, 8                        ! Nguess, kmin, kmax, n_zeros
        write(unit, '(i5)') 0                                        ! use_winding
        ! fstart array (Nguess=3 pairs, each pair on separate line)
        write(unit, '(2(e15.7,1x))') 1.0e6_dp, 0.0_dp      ! fstart[1] = (1.0e6, 0.0)
        write(unit, '(2(e15.7,1x))') 1.5e6_dp, 0.1e6_dp    ! fstart[2] = (1.5e6, 0.1e6)
        write(unit, '(2(e15.7,1x))') 2.0e6_dp, -0.05e6_dp  ! fstart[3] = (2.0e6, -0.05e6)
        close(unit)
        
        print *, "Created legacy format files only"
        
    end subroutine setup_minimal_test
    
end program test_minimal_workflow