program test_gauss_int

    use KIM_kinds, only: dp
    use gauss_quad, only: gauss_config_t, init_gauss_int!, gauss_integrate

    implicit none

    type(gauss_config_t) :: gauss_conf

    print *, "Testing Gauss-Legendre Integration"

    gauss_conf%n = 10

    call init_gauss_int(gauss_conf)
    !call gauss_integrate(integ, 0.0d0, 3.14159265358979d0, result, gauss_conf)

    !if (abs(result-2.0d0) > 1.0d-6) then
        !print *, 'Error: Integral result is incorrect.'
        !stop
    !end if
    !print *, 'Success: Integral result is correct.'

    !contains


end program
