program test_kernel_rho_phi
    ! Behavioral checks of the restored Fourier-space rho-phi kernel:
    ! adiabatic Debye screening, gyroaverage Gaussian, swap conjugacy,
    ! phase translation, Debye switch, and linearity in A1 and A2.

    use KIM_kinds_m, only: dp
    use constants_m, only: pi, com_unit
    use species_m, only: plasma
    use kernels_m, only: kernel_rho_phi_of_kr_krp_rg, kernel_debye_case
    use kernel_test_background_m, only: setup_uniform_background, &
                                        set_forces, require_close

    implicit none

    real(dp) :: kr, krp, rg, shift, debye_sum, gauss_sum, rho_l
    complex(dp) :: kernel_value, swapped, shifted, base, with_a1, with_2a1
    integer :: sp

    call setup_uniform_background()
    rg = plasma%r_grid(50)

    ! Adiabatic screening limit: uniform plasma, zero forces, kr = krp.
    kr = 1.0d0
    krp = 1.0d0
    debye_sum = 0.0d0
    do sp = 0, plasma%n_species - 1
        debye_sum = debye_sum + 1.0d0/plasma%spec(sp)%lambda_D(1)**2
    end do
    kernel_value = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    call require_close("adiabatic Debye screening -1/(8 pi^2 lambda_D^2)", &
                       kernel_value, cmplx(-debye_sum/(8.0d0*pi**2), 0.0d0, dp), &
                       1.0d-12)

    ! Gyroaverage Gaussian and Fourier phase at kr /= krp.
    krp = 3.0d0
    gauss_sum = 0.0d0
    do sp = 0, plasma%n_species - 1
        rho_l = plasma%spec(sp)%vT(1)/abs(plasma%spec(sp)%omega_c(1))
        gauss_sum = gauss_sum + exp(-rho_l**2/2.0d0*(kr - krp)**2) &
                    /plasma%spec(sp)%lambda_D(1)**2
    end do
    kernel_value = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    call require_close("gyroaverage Gaussian with exp(i (kr-krp) rg) phase", &
                       kernel_value, &
                       -gauss_sum/(8.0d0*pi**2)*exp(com_unit*(kr - krp)*rg), &
                       1.0d-12)

    ! Swapping kr and krp conjugates the uniform kernel.
    swapped = kernel_rho_phi_of_kr_krp_rg(krp, kr, rg)
    call require_close("kr <-> krp swap conjugates the uniform kernel", &
                       swapped, conjg(kernel_value), 1.0d-12)

    ! Translating rg only rotates the phase.
    shift = 0.7d0
    shifted = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg + shift)
    call require_close("rg translation rotates phase by exp(i (kr-krp) dr)", &
                       shifted, kernel_value*exp(com_unit*(kr - krp)*shift), &
                       1.0d-12)

    ! Debye switch: positive screening sum with the gyro Gaussian.
    kernel_debye_case = .true.
    kernel_value = kernel_rho_phi_of_kr_krp_rg(kr, kr, rg)
    call require_close("Debye-case switch gives +1/(8 pi^2 lambda_D^2)", &
                       kernel_value, cmplx(debye_sum/(8.0d0*pi**2), 0.0d0, dp), &
                       1.0d-12)
    kernel_debye_case = .false.

    ! Kinetic response is linear in each thermodynamic force.
    base = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    call set_forces(1, 1.0d-2, 0.0d0)
    with_a1 = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    call set_forces(1, 2.0d-2, 0.0d0)
    with_2a1 = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    if (abs(with_a1 - base) <= 0.0d0) then
        print *, "FAIL A1 force does not change the kernel"
        error stop
    end if
    call require_close("kernel is linear in A1", with_2a1 - base, &
                       2.0d0*(with_a1 - base), 1.0d-12)

    call set_forces(1, 0.0d0, 1.0d-2)
    with_a1 = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    call set_forces(1, 0.0d0, 2.0d-2)
    with_2a1 = kernel_rho_phi_of_kr_krp_rg(kr, krp, rg)
    call require_close("kernel is linear in A2", with_2a1 - base, &
                       2.0d0*(with_a1 - base), 1.0d-12)

end program test_kernel_rho_phi
