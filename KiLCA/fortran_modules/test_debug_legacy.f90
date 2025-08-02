program test_debug_legacy
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT, LEGACY_FORMAT
    use kilca_settings_m, only: settings_t, settings_initialize_defaults, &
                                settings_read_all, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(settings_t), pointer :: sd => null()
    integer :: ierr
    
    print *, "Testing legacy workflow debug..."
    
    ! Create test directory with legacy files
    call system("mkdir -p ./test_debug_legacy/")
    call create_legacy_files()
    
    ! Create settings
    call settings_initialize_defaults(sd, "./test_debug_legacy/", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "FAILED: Could not create settings:", ierr
        stop 1
    end if
    
    print *, "Created settings successfully"
    
    ! Disable namelist backend (use legacy)
    call settings_integrate_namelist_backend(.false.)
    print *, "Disabled namelist backend (legacy mode)"
    
    ! Try to read all settings
    call settings_read_all(sd, ierr)
    
    if (ierr == KILCA_SUCCESS) then
        print *, "SUCCESS: Read all settings successfully"
        print *, "Background rtor:", sd%background_settings%rtor
    else
        print *, "FAILED: Error code:", ierr
    end if
    
    ! Clean up
    call system("rm -rf ./test_debug_legacy/")
    
contains

    subroutine create_legacy_files()
        integer :: unit
        
        ! antenna.in
        open(newunit=unit, file="./test_debug_legacy/antenna.in", status="replace")
        write(unit, '(3(e15.7,1x))') 90.0_dp, 5.0_dp, 1.0e12_dp
        write(unit, '(2(e15.7,1x))') 1.0e6_dp, 0.0_dp
        write(unit, '(i5)') 2
        write(unit, '(2(i5,1x))') 0, 1
        write(unit, '(4(i5,1x))') 1, 1, 2, -1
        close(unit)
        
        ! background.in
        open(newunit=unit, file="./test_debug_legacy/background.in", status="replace")
        write(unit, '(3(e15.7,1x))') 170.0_dp, 65.0_dp, 25000.0_dp
        write(unit, '(a)') "./profiles/"
        write(unit, '(a)') "experimental"
        write(unit, '(i5)') 1
        write(unit, '(i5)') 3
        write(unit, '(2(e15.7,1x))') 1.5e9_dp, 0.9_dp
        write(unit, '(3(e15.7,1x))') 2.5_dp, 1.0_dp, 1.0_dp
        write(unit, '(i5)') 0
        write(unit, '(3(e15.7,1x))') 2.0_dp, 1.0_dp, 4.0_dp
        write(unit, '(3(e15.7,1x))') 1.0_dp, -1.0_dp, 2.0_dp
        close(unit)
        
        ! output.in
        open(newunit=unit, file="./test_debug_legacy/output.in", status="replace")
        write(unit, '(5(i5,1x))') 1, 1, 0, 1, 0
        write(unit, '(i5)') 5
        write(unit, '(5(i5,1x))') 1, 1, 0, 1, 1
        close(unit)
        
        ! eigmode.in
        open(newunit=unit, file="./test_debug_legacy/eigmode.in", status="replace")
        write(unit, '(a)') "eigenmode.dat"
        write(unit, '(i5)') 1
        write(unit, '(2(i5,1x))') 120, 60
        write(unit, '(4(e15.7,1x))') 0.0_dp, 3.0e6_dp, -2.0e5_dp, 2.0e5_dp
        write(unit, '(i5)') 0
        write(unit, '(4(e15.7,1x))') 1.0e-7_dp, 1.0e-9_dp, 1.0e-7_dp, 1.0e-4_dp
        write(unit, '(2(i5,1x))') 0, 0
        write(unit, '(4(i5,1x))') 3, 1, 6, 8
        write(unit, '(i5)') 0
        write(unit, '(6(e15.7,1x))') 1.0e6_dp, 0.0_dp, 1.5e6_dp, 0.1e6_dp, 2.0e6_dp, -0.05e6_dp
        close(unit)
        
        print *, "Created legacy format files"
        
    end subroutine create_legacy_files
    
end program test_debug_legacy