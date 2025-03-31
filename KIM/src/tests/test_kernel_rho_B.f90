program test_kernel_rho_B

    use kernels, only: kernel_rho_B_of_kr_krp_rg
    use KIM_kinds, only: dp
    use plasma_parameter, only: r_prof, n_prof, ni_prof, Te_prof, Ti_prof, iprof_length, Er_prof
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

    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)

    if (abs(kernel_value - 0.0d0) < 1.0d-6) then
        print *, "Test constant passed"
        print *, ""
    else
        print *, "Test failed, value: ", kernel_value, " should be: ", 0.0d0
        error stop
    end if

    ne_core = n_prof(1)
    delta_r = r_prof(iprof_length) - r_prof(1)

    do i=2, iprof_length
        n_prof(i) = ne_core * (1.0d0 - r_prof(i) / delta_r)
    end do

    call calculate_backs(.false.)

    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(abs(kernel_value) - 0.0158746) < 1.0d-3) then
        print *, "Test linear n passed"
        print *, ""
    else
        print *, "Test linear n failed, value: ", abs(kernel_value), " should be: 0.0158746"
        error stop
    end if


    do i=2, iprof_length
        n_prof(i) = ne_core
        Te_prof(i) = Te_prof(1) * (1.0d0 - r_prof(i) / delta_r)
    end do

    call calculate_backs(.false.)

    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(abs(kernel_value) - 0.0517352) < 1.0d-3) then
        print *, "Test linear Te passed"
        print *, ""
    else
        print *, "Test linear Te failed, value: ", abs(kernel_value), " should be: 0.0517352"
        error stop
    end if


    do i=2, iprof_length
        Te_prof(i) = Te_prof(1)
        Er_prof(i) = 0.3d0
    end do

    call calculate_backs(.false.)

    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(abs(kernel_value) - 0.184567) < 1.0d-1) then
        print *, "Test constant Er passed"
        print *, ""
    else
        print *, "Test constant Er failed, value: ", abs(kernel_value), " should be: 0.184567"
        error stop
    end if


    do i=2, iprof_length
        Te_prof(i) = Te_prof(1) * (1.0d0 - r_prof(i) / delta_r)
        Er_prof(i) = 0.3d0
    end do

    call calculate_backs(.false.)

    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(abs(kernel_value) - 0.244917) < 1.0d-1) then
        print *, "Test constant Er linear Te passed"
        print *, ""
    else
        print *, "Test constant Er linear Te failed, value: ", abs(kernel_value), " should be: 0.244917"
        error stop
    end if

end program
