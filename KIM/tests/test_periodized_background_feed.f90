program test_periodized_background_feed
    ! Behavioral checks of the periodized-background feed: the resonant
    ! layer keeps the physical profile, the kernels see an L-periodic
    ! background, and gradients come from differentiating the periodized
    ! profile rather than periodizing the physical gradient.

    use KIM_kinds_m, only: dp
    use constants_m, only: e_charge, ev, com_unit
    use species_m, only: plasma, plasma_t, init_electron_species, &
                        interpolate_plasma_backs
    use kernels_m, only: kernel_rho_phi_of_kr_krp_rg
    use rt_fourier_periodic_m, only: build_periodized_background
    use kernel_test_background_m, only: setup_uniform_background, n0, T0, b0, &
                                        ks0, nu0
    use config_m, only: number_of_ion_species, rescale_density, &
                        ion_flr_scale_factor, collision_frequency_scale
    use setup_m, only: collisions_off, omega

    implicit none

    real(dp), parameter :: r_res = 36.5d0
    real(dp), parameter :: dr_layer = 3.0d0
    real(dp), parameter :: dr_transition = 4.0d0
    real(dp), parameter :: period = 2.0d0*(dr_layer + dr_transition)
    ! 1470 cells over three periods: the spacing divides L exactly, so
    ! interpolation at r and r + L uses identical weights.
    integer, parameter :: n_points = 1471
    real(dp), parameter :: slope = -1.0d0/(4.0d0*67.0d0)

    real(dp) :: r, physical, got, centered, physical_slope, deviation
    real(dp) :: lee, lei, expected_nue
    complex(dp) :: kernel_a, kernel_b, expected_z0
    integer :: i, j, sp

    number_of_ion_species = 1
    rescale_density = .false.
    ion_flr_scale_factor = 1.0d0
    collision_frequency_scale = 3.0d0
    collisions_off = .false.
    omega = 0.0d0

    call check_recomputed_z0

    call setup_uniform_background()
    allocate (plasma%Er(plasma%grid_size), plasma%B0(plasma%grid_size), &
              plasma%q(plasma%grid_size))
    plasma%Er = 0.0d0
    plasma%B0 = b0
    plasma%q = 1.5d0
    do sp = 0, plasma%n_species - 1
        plasma%spec(sp)%n = n0*(1.0d0 + slope*(plasma%r_grid - 3.0d0))
        allocate (plasma%spec(sp)%T(plasma%grid_size))
        plasma%spec(sp)%T = T0
    end do

    call build_periodized_background(r_res, dr_layer, dr_transition, n_points)

    j = point_index(r_res)
    associate (ne => plasma%spec(0)%n(j), ni => plasma%spec(1)%n(j), &
            te => plasma%spec(0)%T(j))
        lee = 23.5d0 - log(sqrt(ne)/te**1.25d0) &
            - sqrt(1d-5 + (log(te) - 2.0d0)**2/16.0d0)
        lei = 24.0d0 - log(sqrt(ne)/te)
        expected_nue = 3.0d0*(5.8e-6*ne*lee/te**1.5d0 &
            + 7.7d-6*ni*lei/te**1.5d0)
    end associate
    call require_close("collision-frequency scale reaches the background", &
                       plasma%spec(0)%nu(j), expected_nue, 1.0d-12)

    j = point_index(r_res - 2.0d0)
    expected_z0 = -(plasma%om_E(j) - omega &
        - com_unit*plasma%spec(0)%nu(j)) &
        /(abs(plasma%kp(j))*sqrt(2.0d0)*plasma%spec(0)%vT(j))
    if (abs(plasma%spec(0)%z0(j) - expected_z0) &
            > 1.0d-12*max(abs(expected_z0), 1.0d0)) &
        error stop 'interpolated z0 violates its defining formula'

    ! Inside the resonant layer the linear physical profile survives
    ! exactly (four-point Lagrange interpolation is exact on it).
    do i = 0, 10
        r = r_res - dr_layer + 2.0d0*dr_layer*i/10.0d0
        j = point_index(r)
        physical = n0*(1.0d0 + slope*(plasma%r_grid(j) - 3.0d0))
        call require_close("resonant layer keeps the physical density", &
                           plasma%spec(0)%n(j), physical, 1.0d-10)
    end do

    ! The kernels see an L-periodic background.
    do i = 0, 2
        r = r_res - 0.4d0*period + 0.4d0*period*i
        kernel_a = kernel_rho_phi_of_kr_krp_rg(1.0d0, 1.0d0, r)
        kernel_b = kernel_rho_phi_of_kr_krp_rg(1.0d0, 1.0d0, r + period)
        if (abs(kernel_a - kernel_b) > 1.0d-10*abs(kernel_a)) then
            print *, "FAIL kernel is not periodic at r = ", r
            print *, "  K(r)   = ", kernel_a
            print *, "  K(r+L) = ", kernel_b
            error stop
        end if
    end do
    print *, "PASS kernels see an L-periodic background"

    ! Gradients are differentiated from the periodized profile: on the
    ! uniform feed grid the pipeline stencil is the centered difference.
    j = point_index(r_res + dr_layer + 0.5d0*dr_transition)
    centered = (plasma%spec(0)%n(j + 1) - plasma%spec(0)%n(j - 1)) &
               /(plasma%r_grid(j + 1) - plasma%r_grid(j - 1))
    call require_close("transition gradient is the centered difference &
                       &of the periodized density", &
                       plasma%spec(0)%dndr(j), centered, 1.0d-10)

    physical_slope = n0*slope
    deviation = abs(plasma%spec(0)%dndr(j) - physical_slope) &
                /abs(physical_slope)
    if (deviation < 0.5d0) then
        print *, "FAIL transition gradient looks like the physical slope"
        print *, "  dndr = ", plasma%spec(0)%dndr(j), &
            " physical = ", physical_slope
        error stop
    end if
    print *, "PASS transition gradient differs from the physical slope"

    ! A1 is assembled from the same periodized arrays.
    associate (spec => plasma%spec(0))
        call require_close("A1 matches the stored periodized arrays", &
                           spec%A1(j), &
                           spec%dndr(j)/spec%n(j) &
                           - spec%Zspec*e_charge/(spec%T(j)*ev)*plasma%Er(j) &
                           - 3.0d0/(2.0d0*spec%T(j))*spec%dTdr(j), 1.0d-12)
    end associate

    print *, "PASS all periodized feed checks"

contains

    subroutine check_recomputed_z0

        type(plasma_t) :: coarse
        real(dp) :: target_grid(3)
        complex(dp) :: want
        integer :: k

        coarse%n_species = 1
        coarse%grid_size = 4
        allocate (coarse%spec(0:0))
        call init_electron_species(coarse%spec(0))
        allocate (coarse%r_grid(4), coarse%ks(4), coarse%kp(4), &
                  coarse%om_E(4), coarse%B0(4), coarse%q(4), &
                  coarse%dqdr(4), coarse%Er(4))
        allocate (coarse%spec(0)%n(4), coarse%spec(0)%dndr(4), &
                  coarse%spec(0)%T(4), coarse%spec(0)%dTdr(4), &
                  coarse%spec(0)%nu(4), coarse%spec(0)%vT(4), &
                  coarse%spec(0)%omega_c(4), coarse%spec(0)%lambda_D(4), &
                  coarse%spec(0)%rho_L(4), coarse%spec(0)%z0(4))
        coarse%r_grid = [10.0d0, 20.0d0, 30.0d0, 40.0d0]
        coarse%kp = [3.0d-3, 1.0d-3, -1.0d-3, -3.0d-3]
        coarse%om_E = [9.3d4, 8.0d4, 7.0d4, 6.0d4]
        coarse%ks = ks0
        coarse%B0 = b0
        coarse%q = 1.5d0
        coarse%dqdr = 0.0d0
        coarse%Er = 0.0d0
        coarse%spec(0)%n = n0
        coarse%spec(0)%dndr = 0.0d0
        coarse%spec(0)%T = T0
        coarse%spec(0)%dTdr = 0.0d0
        coarse%spec(0)%nu = nu0
        coarse%spec(0)%vT = 2.0d9
        coarse%spec(0)%omega_c = 1.0d8
        coarse%spec(0)%lambda_D = 1.0d-2
        coarse%spec(0)%rho_L = 1.0d-1
        coarse%spec(0)%z0 = -(coarse%om_E - omega - com_unit*nu0) &
            /(abs(coarse%kp)*sqrt(2.0d0)*coarse%spec(0)%vT)

        target_grid = [15.0d0, 25.0d0, 35.0d0]
        call interpolate_plasma_backs(coarse, target_grid)
        do k = 1, size(target_grid)
            want = -(coarse%om_E(k) - omega - com_unit*coarse%spec(0)%nu(k)) &
                /(abs(coarse%kp(k))*sqrt(2.0d0)*coarse%spec(0)%vT(k))
            if (abs(coarse%spec(0)%z0(k) - want) &
                    > 1.0d-12*max(abs(want), 1.0d0)) &
                error stop 'interpolated z0 violates its defining formula'
        end do

    end subroutine check_recomputed_z0

    integer function point_index(r_target)

        implicit none

        real(dp), intent(in) :: r_target

        point_index = 1 + nint((r_target - plasma%r_grid(1)) &
                               /(plasma%r_grid(2) - plasma%r_grid(1)))

    end function point_index

    subroutine require_close(name, got_value, want_value, rel_tol)

        implicit none

        character(*), intent(in) :: name
        real(dp), intent(in) :: got_value, want_value, rel_tol

        if (abs(got_value - want_value) > rel_tol*max(abs(want_value), 1.0d0)) then
            print *, "FAIL ", name
            print *, "  got  = ", got_value
            print *, "  want = ", want_value
            error stop
        end if
        print *, "PASS ", name

    end subroutine require_close

end program test_periodized_background_feed
