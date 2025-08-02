program test_debug_integration
    use kilca_types_m, only: dp, KILCA_SUCCESS, NAMELIST_FORMAT
    use kilca_settings_m, only: back_sett_t, back_sett_read_settings, &
                                settings_integrate_namelist_backend
    implicit none
    
    type(back_sett_t) :: bs
    integer :: ierr
    
    print *, "Testing namelist backend integration..."
    
    ! Enable namelist backend
    print *, "Enabling namelist backend..."
    call settings_integrate_namelist_backend(.true.)
    
    ! Try to read settings
    print *, "Reading background settings..."
    call back_sett_read_settings(bs, "./", ierr)
    
    if (ierr == KILCA_SUCCESS) then
        print *, "SUCCESS: Settings read successfully"
        print *, "Format used:", bs%format_used
        print *, "rtor value:", bs%rtor
    else
        print *, "FAILED: Error code:", ierr
    end if
    
end program test_debug_integration