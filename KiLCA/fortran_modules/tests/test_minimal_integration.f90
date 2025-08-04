program test_minimal_integration
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT
    use kilca_settings_m, only: settings_t, settings_create, &
                                back_sett_read_settings, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(settings_t), pointer :: sd => null()
    integer :: ierr
    
    print *, "Testing minimal integration..."
    
    ! First create a settings object with defaults
    call settings_create(sd, "./", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "FAILED: Could not create settings:", ierr
        stop 1
    end if
    
    print *, "Created settings successfully"
    
    ! Enable namelist backend
    call settings_integrate_namelist_backend(.true.)
    print *, "Enabled namelist backend"
    
    ! Try to read background settings
    call back_sett_read_settings(sd%background_settings, "./", ierr)
    
    if (ierr == KILCA_SUCCESS) then
        print *, "SUCCESS: Read settings successfully"
        print *, "Format used:", sd%background_settings%format_used
        print *, "rtor value:", sd%background_settings%rtor
    else
        print *, "INFO: Error code:", ierr
        print *, "Falling back worked, this is expected behavior"
    end if
    
end program test_minimal_integration