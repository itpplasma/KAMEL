program test_gauss_segment_transforms

    use constants_m, only: pi
    use grid_m, only: rg_grid
    use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, &
        gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, &
        gauss_int_F3_rho_phi_t, integration_point_t
    use KIM_kinds_m, only: dp
    use species_m, only: plasma

    implicit none

    integer, parameter :: radial_nodes = 100000
    real(dp), parameter :: lower = -0.45_dp
    real(dp), parameter :: upper = 0.8_dp
    real(dp), parameter :: rho_t = 0.7_dp
    real(dp), parameter :: k_s = 0.4_dp
    real(dp), parameter :: theta = 1.1_dp
    real(dp), parameter :: x = -0.3_dp
    real(dp), parameter :: x_prime = 0.5_dp
    real(dp), parameter :: derivative_step = 2.0e-4_dp
    real(dp), parameter :: relative_tolerance = 2.0e-6_dp
    type(integration_point_t) :: point
    type(gauss_int_F0_rho_phi_t) :: integrand_f0
    type(gauss_int_F1_rho_phi_t) :: integrand_f1
    type(gauss_int_F2_rho_phi_t) :: integrand_f2
    type(gauss_int_F3_rho_phi_t) :: integrand_f3
    real(dp) :: actual_f0, actual_f1, actual_f2, actual_f3
    real(dp) :: direct_f0, direct_f1, direct_f2, direct_f3
    integer :: failures

    failures = 0
    call initialize_segment(point)

    integrand_f0%int_point = point
    integrand_f0%int_point%xlpm1 = x - 1.0_dp
    integrand_f0%int_point%xlp = x
    integrand_f0%int_point%xlpp1 = x + 1.0_dp
    integrand_f1%int_point = point

    call point%calc_Jrg1(theta, x, x_prime)
    call point%calc_Jrg23(theta, x, x_prime)
    call point%calc_Jrg4(theta, x, x_prime)
    integrand_f2%int_point = point
    integrand_f3%int_point = point

    actual_f0 = integrand_f0%f(x)
    actual_f1 = integrand_f1%f(x, x_prime, theta)
    actual_f2 = integrand_f2%f(x, x_prime, theta)
    actual_f3 = integrand_f3%f(x, x_prime, theta)
    call direct_segment_values(direct_f0, direct_f1, direct_f2, direct_f3)

    call check_close('F0 compensated normalization', actual_f0, direct_f0, failures)
    call check_close('F1 direct radial definition', actual_f1, direct_f1, failures)
    call check_close('F2 direct derivative definition', actual_f2, direct_f2, failures)
    call check_close('F3 direct derivative definition', actual_f3, direct_f3, failures)

    if (failures /= 0) error stop 'Gauss segment-transform regression failed'
    print *, 'PASS: Gauss F0--F3 segment transforms match direct radial definitions'

contains

    subroutine initialize_segment(point_out)

        type(integration_point_t), intent(out) :: point_out

        rg_grid%npts_b = 2
        rg_grid%npts_c = 1
        allocate(rg_grid%xb(2))
        rg_grid%xb = [lower, upper]

        allocate(plasma%ks(2))
        plasma%ks = k_s

        point_out%j = 1
        point_out%mphi = 0
        point_out%rhoT = rho_t
        point_out%xlm1 = x - 1.0_dp
        point_out%xl = x
        point_out%xlp1 = x + 1.0_dp
        point_out%xlpm1 = x_prime - 1.0_dp
        point_out%xlp = x_prime
        point_out%xlpp1 = x_prime + 1.0_dp

    end subroutine initialize_segment

    subroutine direct_segment_values(f0_value, f1_value, f2_value, f3_value)

        real(dp), intent(out) :: f0_value, f1_value, f2_value, f3_value
        real(dp) :: dr, radius
        integer :: radial_index

        dr = (upper - lower) / real(radial_nodes, dp)
        f0_value = 0.0_dp
        f1_value = 0.0_dp
        f2_value = 0.0_dp
        f3_value = 0.0_dp

        do radial_index = 1, radial_nodes
            radius = lower + (real(radial_index, dp) - 0.5_dp) * dr
            f0_value = f0_value + sqrt(2.0_dp * pi) / rho_t * &
                exp(-(radius - x)**2 / (2.0_dp * rho_t**2))
            f1_value = f1_value + q_density(radius, x, x_prime)
            f2_value = f2_value + f2_density(radius)
            f3_value = f3_value + f3_density(radius)
        end do

        ! The standard global-FEM F0 path intentionally retains this compensated
        ! factor. F1--F3 use the direct segment normalization.
        f0_value = pi * f0_value * dr
        f1_value = f1_value * dr
        f2_value = f2_value * dr
        f3_value = f3_value * dr

    end subroutine direct_segment_values

    pure function q_density(radius, x_value, xp_value) result(value)

        real(dp), intent(in) :: radius, x_value, xp_value
        real(dp) :: value
        real(dp) :: cosine

        cosine = cos(theta)
        value = exp(-k_s**2 * rho_t**2 - &
            (radius - (x_value + xp_value) / 2.0_dp)**2 / &
            (rho_t**2 * (1.0_dp + cosine)) - &
            (x_value - xp_value)**2 / &
            (4.0_dp * rho_t**2 * (1.0_dp - cosine))) / &
            (rho_t**2 * sin(theta))

    end function q_density

    pure function f2_density(radius) result(value)

        real(dp), intent(in) :: radius
        real(dp) :: value
        real(dp) :: d2q_dx2, d2q_dxp2

        d2q_dx2 = second_x_derivative(radius, x, x_prime)
        d2q_dxp2 = second_xp_derivative(radius, x, x_prime)
        value = rho_t**2 * k_s**2 * q_density(radius, x, x_prime) - &
            rho_t**2 * (d2q_dx2 + d2q_dxp2) / 2.0_dp

    end function f2_density

    pure function f3_density(radius) result(value)

        real(dp), intent(in) :: radius
        real(dp) :: value
        real(dp) :: mixed_derivative

        mixed_derivative = (&
            q_density(radius, x + derivative_step, x_prime + derivative_step) - &
            q_density(radius, x + derivative_step, x_prime - derivative_step) - &
            q_density(radius, x - derivative_step, x_prime + derivative_step) + &
            q_density(radius, x - derivative_step, x_prime - derivative_step)) / &
            (4.0_dp * derivative_step**2)
        value = rho_t**2 * cos(theta) * mixed_derivative

    end function f3_density

    pure function second_x_derivative(radius, x_value, xp_value) result(value)

        real(dp), intent(in) :: radius, x_value, xp_value
        real(dp) :: value

        value = (-q_density(radius, x_value + 2.0_dp * derivative_step, xp_value) + &
            16.0_dp * q_density(radius, x_value + derivative_step, xp_value) - &
            30.0_dp * q_density(radius, x_value, xp_value) + &
            16.0_dp * q_density(radius, x_value - derivative_step, xp_value) - &
            q_density(radius, x_value - 2.0_dp * derivative_step, xp_value)) / &
            (12.0_dp * derivative_step**2)

    end function second_x_derivative

    pure function second_xp_derivative(radius, x_value, xp_value) result(value)

        real(dp), intent(in) :: radius, x_value, xp_value
        real(dp) :: value

        value = (-q_density(radius, x_value, xp_value + 2.0_dp * derivative_step) + &
            16.0_dp * q_density(radius, x_value, xp_value + derivative_step) - &
            30.0_dp * q_density(radius, x_value, xp_value) + &
            16.0_dp * q_density(radius, x_value, xp_value - derivative_step) - &
            q_density(radius, x_value, xp_value - 2.0_dp * derivative_step)) / &
            (12.0_dp * derivative_step**2)

    end function second_xp_derivative

    subroutine check_close(label, actual, expected, failures_inout)

        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        integer, intent(inout) :: failures_inout
        real(dp) :: error, scale

        error = abs(actual - expected)
        scale = max(1.0_dp, abs(expected))
        if (error > relative_tolerance * scale) then
            failures_inout = failures_inout + 1
            print '(A,A)', 'FAIL: ', trim(label)
            print '(A,ES24.16)', '  actual:  ', actual
            print '(A,ES24.16)', '  expected:', expected
            print '(A,ES12.4)', '  error:   ', error
        else
            print '(A,A,A,ES12.4)', 'PASS: ', trim(label), ', error = ', error
        end if

    end subroutine check_close

end program test_gauss_segment_transforms
