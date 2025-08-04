program simple_bessel_test
    use iso_fortran_env, only: real64
    use bessel_gsl_m
    implicit none
    
    real(real64) :: result
    integer :: n = 1
    real(real64) :: x = 2.0_real64
    
    write(*, '(A)') "Testing simple bessel_j_n call..."
    
    result = bessel_j_n(n, x)
    
    write(*, '(A, F12.8)') "bessel_j_n(1, 2.0) = ", result
    
end program simple_bessel_test