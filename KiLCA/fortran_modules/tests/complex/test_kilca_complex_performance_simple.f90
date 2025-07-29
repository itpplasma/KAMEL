program test_kilca_complex_performance_simple
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, result
    complex(real64), allocatable :: array1(:), array2(:), result_array(:)
    real(real64) :: real_result
    integer :: n = 100
    
    test_passed = .true.
    
    print *, "Testing complex number performance optimization interface compilation..."
    print *, ""
    
    ! Test 1: Check that performance optimization functions exist and can be called
    print *, "Testing performance optimization function interface existence..."
    
    ! Test functions
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    z2 = cmplx(1.0_real64, 2.0_real64, real64)
    
    ! Test fast magnitude squared
    real_result = cmplx_abs2_fast(z1)
    print *, "PASS: cmplx_abs2_fast interface available"
    
    ! Test fast real division
    result = cmplx_div_real_fast(z1, 2.0_real64)
    print *, "PASS: cmplx_div_real_fast interface available"
    
    ! Test fast multiplication by i
    result = cmplx_mult_i_fast(z1)
    print *, "PASS: cmplx_mult_i_fast interface available"
    
    ! Test fast conjugate multiplication
    result = cmplx_mult_conj_fast(z1, z2)
    print *, "PASS: cmplx_mult_conj_fast interface available"
    
    ! Test fast squaring
    result = cmplx_square_fast(z1)
    print *, "PASS: cmplx_square_fast interface available"
    
    ! Test fast cubing
    result = cmplx_cube_fast(z1)
    print *, "PASS: cmplx_cube_fast interface available"
    
    ! Test fast reciprocal
    result = cmplx_reciprocal_fast(z1)
    print *, "PASS: cmplx_reciprocal_fast interface available"
    
    ! Test cached trigonometric values
    result = cmplx_exp_i_pi_half()
    print *, "PASS: cmplx_exp_i_pi_half interface available"
    
    result = cmplx_exp_i_pi()
    print *, "PASS: cmplx_exp_i_pi interface available"
    
    ! Test in-place operations
    call cmplx_add_inplace(z1, z2)
    print *, "PASS: cmplx_add_inplace interface available"
    
    z1 = cmplx(3.0_real64, 4.0_real64, real64)  ! Reset
    call cmplx_mult_inplace(z1, z2)
    print *, "PASS: cmplx_mult_inplace interface available"
    
    ! Test array operations
    allocate(array1(n), array2(n), result_array(n))
    array1 = cmplx(1.0_real64, 2.0_real64, real64)
    array2 = cmplx(3.0_real64, 4.0_real64, real64)
    
    call cmplx_array_add_fast(array1, array2, result_array)
    print *, "PASS: cmplx_array_add_fast interface available"
    
    call cmplx_array_mult_fast(array1, array2, result_array)
    print *, "PASS: cmplx_array_mult_fast interface available"
    
    deallocate(array1, array2, result_array)
    
    print *, ""
    print *, "Interface functions available:"
    print *, "- cmplx_abs2_fast(z)          : Fast magnitude squared calculation"
    print *, "- cmplx_div_real_fast(z, r)   : Fast division by real number"
    print *, "- cmplx_mult_i_fast(z)        : Fast multiplication by imaginary unit"
    print *, "- cmplx_mult_conj_fast(z1,z2) : Fast conjugate multiplication"
    print *, "- cmplx_square_fast(z)        : Fast squaring operation"
    print *, "- cmplx_cube_fast(z)          : Fast cubing operation"
    print *, "- cmplx_reciprocal_fast(z)    : Fast reciprocal calculation"
    print *, "- cmplx_exp_i_pi_half()       : Pre-computed exp(i*π/2)"
    print *, "- cmplx_exp_i_pi()            : Pre-computed exp(i*π)"
    print *, "- cmplx_add_inplace(z1, z2)   : In-place addition"
    print *, "- cmplx_mult_inplace(z1, z2)  : In-place multiplication"
    print *, "- cmplx_array_add_fast()      : Fast array addition"
    print *, "- cmplx_array_mult_fast()     : Fast array multiplication"
    
    print *, ""
    if (test_passed) then
        print *, "All complex number performance optimization interface tests PASSED!"
    else
        print *, "Some complex number performance optimization interface tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_performance_simple