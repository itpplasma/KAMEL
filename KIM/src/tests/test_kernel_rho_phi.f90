program test_kernel_rho_phi

    use kernels, only: kernel_rho_phi_of_kr_krp_rg
    use KIM_kinds, only: dp
    use plasma_parameter, only: r_prof, n_prof, ni_prof, Te_prof, Ti_prof, iprof_length
    use constants, only: pi, ev, e_charge

    implicit none

    real(dp) :: kr, krp, rg
    complex(dp) :: kernel_value, kernel_test_value
    real(dp) :: lambda_De, lambda_Di, lambda_D
    real(dp) :: ne_core
    real(dp) :: delta_r
    integer :: r_ind = 50
    integer :: i

    call kim_init_for_test

    kr = 1.0d0
    krp = 1.0d0

    rg = r_prof(r_ind)
    print *, "r = ", rg

    kernel_value = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)

    print *, ""
    print *, "kernel value = ", kernel_value

    lambda_De = sqrt(Te_prof(r_ind) * ev / (4.0d0 * pi * e_charge**2 * n_prof(r_ind)))
    lambda_Di = sqrt(Ti_prof(1, r_ind) * ev / (4.0d0 * pi * e_charge**2 * ni_prof(1, r_ind)))
    lambda_D = sqrt(1.0d0/(1.0d0/lambda_De**2 + 1.0d0/lambda_Di**2))

    kernel_test_value = -1.0d0 / (2.0d0**3.0d0 * pi**2.0d0 * lambda_D**2.0d0)

    print *, "lambda De = ", lambda_De
    print *, "lambda Di = ", lambda_Di
    print *, "lambda D = ", lambda_D

    if (abs(kernel_value - kernel_test_value) < 1.0d-6) then
        print *, "Test constant passed"
        print *, ""
    else
        print *, "Test failed, value: ", kernel_value, " should be: ",kernel_test_value
    end if

    ne_core = n_prof(1)
    delta_r = r_prof(iprof_length) - r_prof(1)

    do i=2, iprof_length
        n_prof(i) = ne_core * (1.0d0 - r_prof(i) / delta_r)
    end do

    call calculate_backs(.false.)

    kernel_value = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    print*, "kernel value = ", kernel_value

end program