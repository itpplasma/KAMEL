program test_enhanced_validation_simple
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: back_sett_t, back_sett_validate, back_sett_initialize_defaults
    implicit none
    
    type(back_sett_t) :: bs
    logical :: is_valid
    character(len=1024) :: error_msg
    integer :: ierr
    
    print *, "Testing Enhanced Validation..."
    
    ! Test torus radius range validation
    call back_sett_initialize_defaults(bs, ierr)
    bs%rtor = 30.0_dp  ! Below minimum range (should fail)
    bs%rp = 20.0_dp
    bs%B0 = 25000.0_dp
    bs%m_i = 2.0_dp
    
    call back_sett_validate(bs, is_valid, error_msg, ierr)
    if (is_valid) then
        print *, "FAIL: Small torus radius validation not working"
    else
        print *, "PASS: Small torus radius correctly rejected:", trim(error_msg)
    end if
    
    ! Test ion mass range validation
    bs%rtor = 170.0_dp  ! Fix radius
    bs%m_i = 0.2_dp     ! Below minimum range (should fail)
    
    call back_sett_validate(bs, is_valid, error_msg, ierr)
    if (is_valid) then
        print *, "FAIL: Small ion mass validation not working"
    else
        print *, "PASS: Small ion mass correctly rejected:", trim(error_msg)
    end if
    
    print *, "Enhanced validation test completed"
end program test_enhanced_validation_simple