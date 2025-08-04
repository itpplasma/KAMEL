program test_namelist_direct
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: settings_t, settings_read_namelist
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    
    print *, "Testing direct namelist reading..."
    
    call settings_read_namelist("./settings.conf", settings, ierr)
    
    if (ierr == KILCA_SUCCESS) then
        print *, "SUCCESS: Settings read successfully"
        print *, "rtor value:", settings%background_settings%rtor
        print *, "ra value:", settings%antenna_settings%ra
    else
        print *, "FAILED: Error code:", ierr
    end if
    
end program test_namelist_direct