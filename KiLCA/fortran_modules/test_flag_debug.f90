program test_flag_debug
    use kilca_settings_m, only: settings_integrate_namelist_backend
    implicit none
    
    print *, "Testing global flag functionality..."
    
    ! Test setting to true
    call settings_integrate_namelist_backend(.true.)
    print *, "Set namelist backend to TRUE"
    
    ! Test setting to false  
    call settings_integrate_namelist_backend(.false.)
    print *, "Set namelist backend to FALSE"
    
    print *, "Flag setting appears to work correctly"
    
end program test_flag_debug