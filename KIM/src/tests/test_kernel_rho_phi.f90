program test_kernel_rho_phi

    use kernels, only: kernel_rho_phi_of_kr_krp_rg
    use KIM_kinds, only: dp
    use plasma_parameter, only: r_prof, n_prof, ni_prof, Te_prof, Ti_prof
    use constants, only: pi, ev, e_charge

    implicit none

    real(dp) :: kr, krp, rg
    real(dp) :: kernel_value, kernel_test_value
    real(dp) :: lambda_De, lambda_Di, lambda_D
    integer :: r_ind = 50

    call kim_init_for_test

    kr = 100.0d0
    krp = 100.0d0

    rg = r_prof(r_ind)

    kernel_value = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)

    print *, ""
    print *, "kernel value = ", kernel_value

    lambda_De = sqrt(Te_prof(r_ind) * ev / (4.0d0 * pi * e_charge**2 * n_prof(r_ind)))
    lambda_Di = sqrt(Ti_prof(1, r_ind) * ev / (4.0d0 * pi * e_charge**2 * ni_prof(1, r_ind)))
    lambda_D = sqrt(1.0d0/(1.0d0/lambda_De**2 + 1.0d0/lambda_Di**2))

    kernel_test_value = -1.0d0 / (2.0d0**3.0d0 * pi**2.0d0 * lambda_D**2.0d0)

    if (abs(kernel_value - kernel_test_value) < 1.0d-6) then
        print *, "Test passed"
    else
        print *, "Test failed, value: ", kernel_value, " should be: ",kernel_test_value
    end if

    


end program