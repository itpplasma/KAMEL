module radial_current_fourier_kernel_m
    use KIM_kinds_m, only: dp

    implicit none
    private

    type :: flr_jet_t
        real(dp) :: value = 0.0_dp
        real(dp) :: source = 0.0_dp
        real(dp) :: response = 0.0_dp
        real(dp) :: mixed = 0.0_dp
    end type flr_jet_t

    public :: radial_flr_coefficients
    public :: scaled_bessel_harmonic

contains

    real(dp) function scaled_bessel_harmonic(ell, bplus, bcross) result(value)
        !> Return exp(-bplus) I_ell(bcross) without overflowing the unscaled
        !> modified Bessel function. Integer-order symmetry is applied before
        !> the large-argument expansion.
        use constants_m, only: pi
        use fortnum_special, only: bessel_in

        integer, intent(in) :: ell
        real(dp), intent(in) :: bplus, bcross
        real(dp), parameter :: asymptotic_threshold = 500.0_dp
        real(dp) :: mu, inverse_eight_x, series
        integer :: order

        order = abs(ell)
        if (bcross == 0.0_dp) then
            if (order == 0) then
                value = exp(-bplus)
            else
                value = 0.0_dp
            end if
            return
        end if

        if (bcross <= asymptotic_threshold) then
            value = exp(-bplus) * bessel_in(order, bcross)
            return
        end if

        mu = 4.0_dp * real(order * order, dp)
        inverse_eight_x = 1.0_dp / (8.0_dp * bcross)
        series = 1.0_dp &
            - (mu - 1.0_dp) * inverse_eight_x &
            + 0.5_dp * (mu - 1.0_dp) * (mu - 9.0_dp) * inverse_eight_x**2 &
            - (mu - 1.0_dp) * (mu - 9.0_dp) * (mu - 25.0_dp) * &
            inverse_eight_x**3 / 6.0_dp
        value = exp(bcross - bplus) * series / sqrt(2.0_dp * pi * bcross)
    end function scaled_bessel_harmonic

    subroutine radial_flr_coefficients(ell, ks, kr, ksp, krp, rho_l, a1, a2, &
        o0, o2, w0, w2)
        !> Evaluate the phase-stripped coefficients with one (O) and two (W)
        !> radial-velocity insertions. The arguments before rho_l are the
        !> source and response perpendicular wave vectors, respectively.
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks, kr, ksp, krp, rho_l, a1, a2
        complex(dp), intent(out) :: o0, o2, w0, w2
        type(flr_jet_t) :: a0_jet, a2_jet
        real(dp) :: kperp, kperp_p

        kperp = hypot(ks, kr)
        kperp_p = hypot(ksp, krp)

        if (kperp <= tiny(1.0_dp) .and. kperp_p <= tiny(1.0_dp)) then
            call zero_wavenumber_coefficients(ell, rho_l, a1, a2, o0, o2, w0, w2)
            return
        end if
        if (kperp <= tiny(1.0_dp)) then
            call zero_source_coefficients(ell, ksp, krp, rho_l, a1, a2, &
                o0, o2, w0, w2)
            return
        end if
        if (kperp_p <= tiny(1.0_dp)) then
            call zero_response_coefficients(ell, ks, kr, rho_l, a1, a2, &
                o0, o2, w0, w2)
            return
        end if

        call ordinary_coefficient_jets(ell, ks, kr, ksp, krp, rho_l, a1, a2, &
            a0_jet, a2_jet)
        call apply_radial_operators(ell, kr, krp, kperp, kperp_p, a0_jet, o0, w0)
        call apply_radial_operators(ell, kr, krp, kperp, kperp_p, a2_jet, o2, w2)
    end subroutine radial_flr_coefficients

    subroutine ordinary_coefficient_jets(ell, ks, kr, ksp, krp, rho_l, a1, a2, &
        a0_jet, a2_jet)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks, kr, ksp, krp, rho_l, a1, a2
        type(flr_jet_t), intent(out) :: a0_jet, a2_jet
        type(flr_jet_t) :: bplus, bcross, gamma, lambda, gradient_factor
        real(dp) :: kperp, kperp_p, rho_l_squared

        rho_l_squared = rho_l**2
        kperp = hypot(ks, kr)
        kperp_p = hypot(ksp, krp)

        bplus%value = 0.5_dp * rho_l_squared * (kperp**2 + kperp_p**2)
        bplus%source = rho_l_squared * ks
        bplus%response = rho_l_squared * ksp

        bcross%value = rho_l_squared * kperp * kperp_p
        bcross%source = rho_l_squared * ks * kperp_p / kperp
        bcross%response = rho_l_squared * ksp * kperp / kperp_p
        bcross%mixed = rho_l_squared * ks * ksp / (kperp * kperp_p)

        gamma = scaled_bessel_jet(ell, bplus, bcross)
        lambda = multiply_jets(bcross, scaled_bessel_jet(ell - 1, bplus, bcross))

        gradient_factor%value = a1 + (1.0_dp - bplus%value - real(ell, dp)) * a2
        gradient_factor%source = -a2 * bplus%source
        gradient_factor%response = -a2 * bplus%response

        a0_jet = add_jets(multiply_jets(gamma, gradient_factor), &
            scale_jet(lambda, a2))
        a2_jet = scale_jet(gamma, 0.5_dp * a2)
    end subroutine ordinary_coefficient_jets

    function scaled_bessel_jet(ell, bplus, bcross) result(jet)
        integer, intent(in) :: ell
        type(flr_jet_t), intent(in) :: bplus, bcross
        type(flr_jet_t) :: jet
        real(dp) :: centre, first_derivative, second_derivative

        centre = scaled_bessel_harmonic(ell, bplus%value, bcross%value)
        first_derivative = 0.5_dp * &
            (scaled_bessel_harmonic(ell - 1, bplus%value, bcross%value) &
            + scaled_bessel_harmonic(ell + 1, bplus%value, bcross%value))
        second_derivative = 0.25_dp * &
            (scaled_bessel_harmonic(ell - 2, bplus%value, bcross%value) &
            + 2.0_dp * centre &
            + scaled_bessel_harmonic(ell + 2, bplus%value, bcross%value))

        jet%value = centre
        jet%source = -bplus%source * centre + bcross%source * first_derivative
        jet%response = -bplus%response * centre + bcross%response * first_derivative
        jet%mixed = bplus%source * bplus%response * centre &
            - (bplus%source * bcross%response + bplus%response * bcross%source) &
            * first_derivative + bcross%source * bcross%response * second_derivative &
            + bcross%mixed * first_derivative
    end function scaled_bessel_jet

    pure function multiply_jets(left, right) result(product)
        type(flr_jet_t), intent(in) :: left, right
        type(flr_jet_t) :: product

        product%value = left%value * right%value
        product%source = left%source * right%value + left%value * right%source
        product%response = left%response * right%value + left%value * right%response
        product%mixed = left%mixed * right%value + left%source * right%response &
            + left%response * right%source + left%value * right%mixed
    end function multiply_jets

    pure function add_jets(left, right) result(sum)
        type(flr_jet_t), intent(in) :: left, right
        type(flr_jet_t) :: sum

        sum%value = left%value + right%value
        sum%source = left%source + right%source
        sum%response = left%response + right%response
        sum%mixed = left%mixed + right%mixed
    end function add_jets

    pure function scale_jet(jet, scale) result(scaled)
        type(flr_jet_t), intent(in) :: jet
        real(dp), intent(in) :: scale
        type(flr_jet_t) :: scaled

        scaled%value = scale * jet%value
        scaled%source = scale * jet%source
        scaled%response = scale * jet%response
        scaled%mixed = scale * jet%mixed
    end function scale_jet

    pure subroutine apply_radial_operators(ell, kr, krp, kperp, kperp_p, &
        coefficient, observation, double_insertion)
        integer, intent(in) :: ell
        real(dp), intent(in) :: kr, krp, kperp, kperp_p
        type(flr_jet_t), intent(in) :: coefficient
        complex(dp), intent(out) :: observation, double_insertion
        complex(dp), parameter :: imaginary_unit = (0.0_dp, 1.0_dp)
        real(dp) :: source_phase_derivative, response_phase_derivative

        source_phase_derivative = real(ell, dp) * kr / kperp**2
        response_phase_derivative = real(ell, dp) * krp / kperp_p**2
        observation = coefficient%response &
            - imaginary_unit * response_phase_derivative * coefficient%value
        double_insertion = coefficient%mixed &
            + imaginary_unit * source_phase_derivative * coefficient%response &
            - imaginary_unit * response_phase_derivative * coefficient%source &
            + source_phase_derivative * response_phase_derivative * coefficient%value
    end subroutine apply_radial_operators

    pure subroutine zero_wavenumber_coefficients(ell, rho_l, a1, a2, &
        o0, o2, w0, w2)
        integer, intent(in) :: ell
        real(dp), intent(in) :: rho_l, a1, a2
        complex(dp), intent(out) :: o0, o2, w0, w2

        o0 = (0.0_dp, 0.0_dp)
        o2 = (0.0_dp, 0.0_dp)
        w0 = (0.0_dp, 0.0_dp)
        w2 = (0.0_dp, 0.0_dp)
        if (abs(ell) == 1) then
            w0 = cmplx(rho_l**2 * (0.5_dp * a1 + a2), 0.0_dp, dp)
            w2 = cmplx(0.25_dp * rho_l**2 * a2, 0.0_dp, dp)
        end if
    end subroutine zero_wavenumber_coefficients

    pure subroutine zero_source_coefficients(ell, ksp, krp, rho_l, a1, a2, &
        o0, o2, w0, w2)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ksp, krp, rho_l, a1, a2
        complex(dp), intent(out) :: o0, o2, w0, w2
        complex(dp) :: inverse_phase
        real(dp) :: bplus, coefficient0, coefficient2, derivative0, derivative2
        real(dp) :: kperp_p, rho_l_squared

        rho_l_squared = rho_l**2
        kperp_p = hypot(ksp, krp)
        bplus = 0.5_dp * rho_l_squared * kperp_p**2
        o0 = (0.0_dp, 0.0_dp)
        o2 = (0.0_dp, 0.0_dp)
        w0 = (0.0_dp, 0.0_dp)
        w2 = (0.0_dp, 0.0_dp)

        if (ell == 0) then
            o0 = cmplx(-rho_l_squared * ksp * exp(-bplus) &
                * (a1 + (2.0_dp - bplus) * a2), 0.0_dp, dp)
            o2 = cmplx(-0.5_dp * rho_l_squared * ksp * exp(-bplus) * a2, &
                0.0_dp, dp)
        end if
        if (abs(ell) /= 1) return

        coefficient0 = 0.5_dp * rho_l_squared * exp(-bplus) &
            * (a1 + (2.0_dp - bplus) * a2)
        coefficient2 = 0.25_dp * rho_l_squared * exp(-bplus) * a2
        derivative0 = -0.5_dp * rho_l_squared * exp(-bplus) &
            * (a1 + (3.0_dp - bplus) * a2)
        derivative2 = -coefficient2
        inverse_phase = cmplx(ksp, -real(ell, dp) * krp, dp) / kperp_p
        w0 = coefficient0 * inverse_phase &
            + rho_l_squared * ksp * kperp_p * derivative0
        w2 = coefficient2 * inverse_phase &
            + rho_l_squared * ksp * kperp_p * derivative2
    end subroutine zero_source_coefficients

    pure subroutine zero_response_coefficients(ell, ks, kr, rho_l, a1, a2, &
        o0, o2, w0, w2)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks, kr, rho_l, a1, a2
        complex(dp), intent(out) :: o0, o2, w0, w2
        complex(dp) :: inverse_phase
        real(dp) :: bplus, coefficient0, coefficient2, derivative0, derivative2
        real(dp) :: kperp, rho_l_squared

        rho_l_squared = rho_l**2
        kperp = hypot(ks, kr)
        bplus = 0.5_dp * rho_l_squared * kperp**2
        o0 = (0.0_dp, 0.0_dp)
        o2 = (0.0_dp, 0.0_dp)
        w0 = (0.0_dp, 0.0_dp)
        w2 = (0.0_dp, 0.0_dp)
        if (abs(ell) /= 1) return

        coefficient0 = 0.5_dp * rho_l_squared * exp(-bplus) &
            * (a1 + (2.0_dp - bplus) * a2)
        coefficient2 = 0.25_dp * rho_l_squared * exp(-bplus) * a2
        derivative0 = -0.5_dp * rho_l_squared * exp(-bplus) &
            * (a1 + (3.0_dp - bplus) * a2)
        derivative2 = -coefficient2
        o0 = cmplx(coefficient0 * kperp, 0.0_dp, dp)
        o2 = cmplx(coefficient2 * kperp, 0.0_dp, dp)
        inverse_phase = cmplx(ks, real(ell, dp) * kr, dp) / kperp
        w0 = coefficient0 * inverse_phase &
            + rho_l_squared * ks * kperp * derivative0
        w2 = coefficient2 * inverse_phase &
            + rho_l_squared * ks * kperp * derivative2
    end subroutine zero_response_coefficients

end module radial_current_fourier_kernel_m
