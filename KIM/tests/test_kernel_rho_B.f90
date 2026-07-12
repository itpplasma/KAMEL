program test_kernel_rho_B
    ! Behavioral checks of the restored Fourier-space rho-B kernel:
    ! exact zero without thermodynamic forces, linearity in A1 and A2,
    ! and the swap phase relation.

    use KIM_kinds_m, only: dp
    use constants_m, only: com_unit
    use species_m, only: plasma
    use kernels_m, only: kernel_rho_B_of_kr_krp_rg
    use kernel_test_background_m, only: setup_uniform_background, &
                                        set_forces, require_close

    implicit none

    real(dp) :: kr, krp, rg
    complex(dp) :: kernel_value, with_a1, with_2a1, swapped

    call setup_uniform_background()
    rg = plasma%r_grid(50)
    kr = 1.0d0
    krp = 3.0d0

    ! Without thermodynamic forces the magnetic drive vanishes exactly.
    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(kernel_value) > 0.0d0) then
        print *, "FAIL zero-force kernel is not exactly zero: ", kernel_value
        error stop
    end if
    print *, "PASS zero thermodynamic forces give exactly zero"

    ! Linear in A1 at fixed background.
    call set_forces(1, 1.0d-2, 0.0d0)
    with_a1 = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(with_a1) <= 0.0d0) then
        print *, "FAIL A1 force does not drive the kernel"
        error stop
    end if
    call set_forces(1, 2.0d-2, 0.0d0)
    with_2a1 = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    call require_close("kernel is linear in A1", with_2a1, 2.0d0*with_a1, &
                       1.0d-12)

    ! Linear in A2 at fixed background.
    call set_forces(1, 0.0d0, 1.0d-2)
    with_a1 = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    if (abs(with_a1) <= 0.0d0) then
        print *, "FAIL A2 force does not drive the kernel"
        error stop
    end if
    call set_forces(1, 0.0d0, 2.0d-2)
    with_2a1 = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    call require_close("kernel is linear in A2", with_2a1, 2.0d0*with_a1, &
                       1.0d-12)

    ! kr <-> krp swap only reverses the Fourier phase.
    kernel_value = kernel_rho_B_of_kr_krp_rg(kr, krp, rg)
    swapped = kernel_rho_B_of_kr_krp_rg(krp, kr, rg)
    call require_close("kr <-> krp swap reverses the phase", swapped, &
                       kernel_value*exp(-2.0d0*com_unit*(kr - krp)*rg), &
                       1.0d-12)

end program test_kernel_rho_B
