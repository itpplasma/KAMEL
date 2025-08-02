program test_debug_individual
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT, LEGACY_FORMAT
    use kilca_settings_m, only: antenna_t, back_sett_t, output_sett_t, eigmode_sett_t, &
                                antenna_read_settings, back_sett_read_settings, &
                                output_read_settings, eigmode_read_settings, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(antenna_t) :: ant
    type(back_sett_t) :: bg
    type(output_sett_t) :: out
    type(eigmode_sett_t) :: eig
    integer :: ierr
    character(len=*), parameter :: test_path = "./test_debug_legacy/"
    
    print *, "Testing individual settings readers..."
    
    ! Disable namelist backend
    call settings_integrate_namelist_backend(.false.)
    print *, "Disabled namelist backend"
    
    ! Test antenna settings
    print *, "Testing antenna_read_settings..."
    call antenna_read_settings(ant, test_path, ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "PASS: antenna_read_settings"
    else
        print *, "FAIL: antenna_read_settings, error:", ierr
    end if
    
    ! Test background settings  
    print *, "Testing back_sett_read_settings..."
    call back_sett_read_settings(bg, test_path, ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "PASS: back_sett_read_settings"
    else
        print *, "FAIL: back_sett_read_settings, error:", ierr
    end if
    
    ! Test output settings
    print *, "Testing output_read_settings..."
    call output_read_settings(out, test_path, ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "PASS: output_read_settings"
    else
        print *, "FAIL: output_read_settings, error:", ierr
    end if
    
    ! Test eigenmode settings
    print *, "Testing eigmode_read_settings..."
    call eigmode_read_settings(eig, test_path, ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "PASS: eigmode_read_settings"
    else
        print *, "FAIL: eigmode_read_settings, error:", ierr
    end if
    
end program test_debug_individual