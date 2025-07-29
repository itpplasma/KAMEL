program test_kilca_complex_hypergeometric_simple
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2
    
    test_passed = .true.
    
    print *, "Testing complex hypergeometric function interface compilation..."
    print *, ""
    
    ! Test 1: Check that functions exist and can be called
    print *, "Testing function interface existence..."
    
    ! Test function calls
    z1 = cmplx(0.5_real64, 0.0_real64, real64)
    
    ! Test 1F1 function interface
    z2 = cmplx_hyp1f1(1.0_real64, 1.0_real64, z1)
    print *, "PASS: cmplx_hyp1f1 interface available"
    
    ! Test 1F1 derivative interface
    z2 = cmplx_hyp1f1_derivative(1.0_real64, 2.0_real64, z1, 1)
    print *, "PASS: cmplx_hyp1f1_derivative interface available"
    
    ! Test 2F1 function interface
    z2 = cmplx_hyp2f1(1.0_real64, 1.0_real64, 2.0_real64, z1)
    print *, "PASS: cmplx_hyp2f1 interface available"
    
    ! Test 0F1 function interface
    z2 = cmplx_hyp0f1(2.0_real64, z1)
    print *, "PASS: cmplx_hyp0f1 interface available"
    
    ! Test U function interface
    z2 = cmplx_hypU(1.0_real64, 2.0_real64, z1)
    print *, "PASS: cmplx_hypU interface available"
    
    print *, ""
    print *, "Interface functions available:"
    print *, "- cmplx_hyp1f1(a, b, z)           : Confluent hypergeometric 1F1(a,b,z)"
    print *, "- cmplx_hyp1f1_derivative(a,b,z,n): n-th derivative of 1F1(a,b,z)"
    print *, "- cmplx_hyp2f1(a, b, c, z)        : Gauss hypergeometric 2F1(a,b,c,z)"
    print *, "- cmplx_hyp0f1(b, z)              : Confluent hypergeometric 0F1(;b;z)"
    print *, "- cmplx_hypU(a, b, z)             : Tricomi confluent hypergeometric U(a,b,z)"
    
    print *, ""
    if (test_passed) then
        print *, "All hypergeometric function interface tests PASSED!"
    else
        print *, "Some hypergeometric function interface tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_hypergeometric_simple