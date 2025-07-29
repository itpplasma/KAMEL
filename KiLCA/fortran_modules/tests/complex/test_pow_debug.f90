program test_pow_debug
    use iso_fortran_env, only: real64
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    complex(real64) :: z, result
    real(real64) :: magnitude
    
    ! Test i^pi
    z = cmplx_I
    result = cmplx_pow_real(z, pi)
    magnitude = abs(result)
    
    print *, "i = ", z
    print *, "i^pi = ", result
    print *, "|i^pi| = ", magnitude
    print *, "Expected: |i^pi| = exp(-pi^2/2) = ", exp(-pi**2/2.0_real64)
    
    ! Let's verify: i^pi = exp(pi * log(i)) = exp(pi * i*pi/2) = exp(i*pi^2/2)
    ! So |i^pi| = |exp(i*pi^2/2)| = 1, not < 1
    
    print *, ""
    print *, "Actually, |exp(i*theta)| = 1 for any real theta"
    print *, "So |i^pi| should equal 1, not be < 1"
    
end program test_pow_debug