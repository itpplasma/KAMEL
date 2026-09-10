module radial_current_fourier_kernel_m
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t

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
    public :: hatG_jrad_phi, hatG_jrad_br, hatG_jrad_bparallel
    public :: hatG_jrad_phi_harmonic, hatG_jrad_br_harmonic
    public :: hatG_jrad_bparallel_harmonic
    public :: hatG_jrad_phi_harmonic_sp, hatG_jrad_br_harmonic_sp
    public :: hatG_jrad_bparallel_harmonic_sp

contains

    complex(dp) function hatG_jrad_phi(plasma_in, kr_response, kr_source, j, &
        cutoff) result(kernel)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr_response, kr_source
        integer, intent(in) :: j
        integer, intent(in), optional :: cutoff
        integer :: ell, harmonic_cutoff

        harmonic_cutoff = checked_harmonic_cutoff(cutoff)
        kernel = (0.0_dp, 0.0_dp)
        do ell = -harmonic_cutoff, harmonic_cutoff
            kernel = kernel + hatG_jrad_phi_harmonic(plasma_in, ell, &
                kr_response, kr_source, j)
        end do
    end function hatG_jrad_phi

    complex(dp) function hatG_jrad_br(plasma_in, kr_response, kr_source, j, &
        cutoff) result(kernel)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr_response, kr_source
        integer, intent(in) :: j
        integer, intent(in), optional :: cutoff
        integer :: ell, harmonic_cutoff

        harmonic_cutoff = checked_harmonic_cutoff(cutoff)
        kernel = (0.0_dp, 0.0_dp)
        do ell = -harmonic_cutoff, harmonic_cutoff
            kernel = kernel + hatG_jrad_br_harmonic(plasma_in, ell, &
                kr_response, kr_source, j)
        end do
    end function hatG_jrad_br

    complex(dp) function hatG_jrad_bparallel(plasma_in, kr_response, kr_source, j, &
        cutoff) result(kernel)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr_response, kr_source
        integer, intent(in) :: j
        integer, intent(in), optional :: cutoff
        integer :: ell, harmonic_cutoff

        harmonic_cutoff = checked_harmonic_cutoff(cutoff)
        kernel = (0.0_dp, 0.0_dp)
        do ell = -harmonic_cutoff, harmonic_cutoff
            kernel = kernel + hatG_jrad_bparallel_harmonic(plasma_in, ell, &
                kr_response, kr_source, j)
        end do
    end function hatG_jrad_bparallel

    complex(dp) function hatG_jrad_phi_harmonic(plasma_in, ell, kr_response, &
            kr_source, j) result(kernel)
        use config_m, only: turn_off_electrons, turn_off_ions

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: ell, j
        real(dp), intent(in) :: kr_response, kr_source
        integer :: sp

        kernel = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_electrons .and. sp == 0) cycle
            if (turn_off_ions .and. sp >= 1) cycle
            if (use_collisionless_ion(sp)) then
                kernel = kernel + collisionless_jrad_phi_harmonic_sp(plasma_in, &
                    sp, ell, kr_response, kr_source, j)
            else
                kernel = kernel + hatG_jrad_phi_harmonic_sp(plasma_in, sp, ell, &
                    kr_response, kr_source, j)
            end if
        end do
    end function hatG_jrad_phi_harmonic

    complex(dp) function hatG_jrad_br_harmonic(plasma_in, ell, kr_response, &
            kr_source, j) result(kernel)
        use config_m, only: turn_off_electrons, turn_off_ions

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: ell, j
        real(dp), intent(in) :: kr_response, kr_source
        integer :: sp

        kernel = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_electrons .and. sp == 0) cycle
            if (turn_off_ions .and. sp >= 1) cycle
            if (use_collisionless_ion(sp)) then
                kernel = kernel + collisionless_jrad_br_harmonic_sp(plasma_in, &
                    sp, ell, kr_response, kr_source, j)
            else
                kernel = kernel + hatG_jrad_br_harmonic_sp(plasma_in, sp, ell, &
                    kr_response, kr_source, j)
            end if
        end do
    end function hatG_jrad_br_harmonic

    complex(dp) function hatG_jrad_bparallel_harmonic(plasma_in, ell, kr_response, &
            kr_source, j) result(kernel)
        use config_m, only: turn_off_electrons, turn_off_ions

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: ell, j
        real(dp), intent(in) :: kr_response, kr_source
        integer :: sp

        kernel = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_electrons .and. sp == 0) cycle
            if (turn_off_ions .and. sp >= 1) cycle
            if (use_collisionless_ion(sp)) then
                kernel = kernel + collisionless_jrad_bparallel_harmonic_sp( &
                    plasma_in, sp, ell, kr_response, kr_source, j)
            else
                kernel = kernel + hatG_jrad_bparallel_harmonic_sp(plasma_in, &
                    sp, ell, kr_response, kr_source, j)
            end if
        end do
    end function hatG_jrad_bparallel_harmonic

    complex(dp) function hatG_jrad_phi_harmonic_sp(plasma_in, sp, ell, &
        kr_response, kr_source, j) result(kernel)
        use config_m, only: artificial_debye_case

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp) :: o0, o2, w0, w2

        kernel = (0.0_dp, 0.0_dp)
        if (artificial_debye_case /= 0 .and. artificial_debye_case /= 2) return
        call check_harmonic_is_available(plasma_in, sp, ell)
        call radial_flr_coefficients(ell, plasma_in%ks(j), kr_source, &
            plasma_in%ks(j), kr_response, plasma_in%spec(sp)%rho_L(j), &
            plasma_in%spec(sp)%A1(j), plasma_in%spec(sp)%A2(j), o0, o2, w0, w2)
        kernel = -radial_kernel_phase(ell, kr_response, kr_source, plasma_in%ks(j), j) &
            * plasma_in%spec(sp)%vT(j)**2 * plasma_in%ks(j) &
            / (plasma_in%spec(sp)%lambda_D(j)**2 * plasma_in%spec(sp)%nu(j)) &
            * (o0 * plasma_in%spec(sp)%I00(j, ell) &
            + o2 * plasma_in%spec(sp)%I02(j, ell))
    end function hatG_jrad_phi_harmonic_sp

    complex(dp) function hatG_jrad_br_harmonic_sp(plasma_in, sp, ell, &
        kr_response, kr_source, j) result(kernel)
        use config_m, only: artificial_debye_case
        use constants_m, only: com_unit, sol

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp) :: o0, o2, w0, w2

        kernel = (0.0_dp, 0.0_dp)
        if (artificial_debye_case /= 0 .and. artificial_debye_case /= 2) return
        call check_harmonic_is_available(plasma_in, sp, ell)
        call radial_flr_coefficients(ell, plasma_in%ks(j), kr_source, &
            plasma_in%ks(j), kr_response, plasma_in%spec(sp)%rho_L(j), &
            plasma_in%spec(sp)%A1(j), plasma_in%spec(sp)%A2(j), o0, o2, w0, w2)
        kernel = -com_unit &
            * radial_kernel_phase(ell, kr_response, kr_source, plasma_in%ks(j), j) &
            * plasma_in%spec(sp)%vT(j)**3 &
            / (plasma_in%spec(sp)%lambda_D(j)**2 * plasma_in%spec(sp)%nu(j) * sol) &
            * (o0 * plasma_in%spec(sp)%I01(j, ell) &
            + o2 * plasma_in%spec(sp)%I03(j, ell))
    end function hatG_jrad_br_harmonic_sp

    complex(dp) function hatG_jrad_bparallel_harmonic_sp(plasma_in, sp, ell, &
        kr_response, kr_source, j) result(kernel)
        use config_m, only: artificial_debye_case
        use constants_m, only: sol

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp) :: o0, o2, w0, w2

        kernel = (0.0_dp, 0.0_dp)
        if (artificial_debye_case /= 0 .and. artificial_debye_case /= 2) return
        call check_harmonic_is_available(plasma_in, sp, ell)
        call radial_flr_coefficients(ell, plasma_in%ks(j), kr_source, &
            plasma_in%ks(j), kr_response, plasma_in%spec(sp)%rho_L(j), &
            plasma_in%spec(sp)%A1(j), plasma_in%spec(sp)%A2(j), o0, o2, w0, w2)
        kernel = radial_kernel_phase(ell, kr_response, kr_source, plasma_in%ks(j), j) &
            * plasma_in%spec(sp)%vT(j)**2 * plasma_in%spec(sp)%omega_c(j) &
            / (plasma_in%spec(sp)%lambda_D(j)**2 * plasma_in%spec(sp)%nu(j) * sol) &
            * (w0 * plasma_in%spec(sp)%I00(j, ell) &
            + w2 * plasma_in%spec(sp)%I02(j, ell))
    end function hatG_jrad_bparallel_harmonic_sp

    complex(dp) function collisionless_jrad_phi_harmonic_sp(plasma_in, sp, ell, &
            kr_response, kr_source, j) result(kernel)
        use constants_m, only: com_unit, pi

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp) :: even_o, even_w, odd_o

        call collisionless_radial_moments(plasma_in, sp, ell, kr_response, &
            kr_source, j, even_o, odd_o, even_w)
        kernel = com_unit &
            * radial_kernel_phase(ell, kr_response, kr_source, plasma_in%ks(j), j) &
            * plasma_in%ks(j) * plasma_in%spec(sp)%vT(j) &
            / (sqrt(2.0_dp * pi) * plasma_in%spec(sp)%lambda_D(j)**2) * even_o
    end function collisionless_jrad_phi_harmonic_sp

    complex(dp) function collisionless_jrad_br_harmonic_sp(plasma_in, sp, ell, &
            kr_response, kr_source, j) result(kernel)
        use constants_m, only: pi, sol

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp) :: even_o, even_w, odd_o

        call collisionless_radial_moments(plasma_in, sp, ell, kr_response, &
            kr_source, j, even_o, odd_o, even_w)
        kernel = -radial_kernel_phase(ell, kr_response, kr_source, plasma_in%ks(j), j) &
            * plasma_in%spec(sp)%vT(j) &
            / (sqrt(2.0_dp * pi) * plasma_in%spec(sp)%lambda_D(j)**2 * sol) * odd_o
    end function collisionless_jrad_br_harmonic_sp

    complex(dp) function collisionless_jrad_bparallel_harmonic_sp(plasma_in, sp, ell, &
            kr_response, kr_source, j) result(kernel)
        use constants_m, only: com_unit, pi, sol

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp) :: even_o, even_w, odd_o

        call collisionless_radial_moments(plasma_in, sp, ell, kr_response, &
            kr_source, j, even_o, odd_o, even_w)
        kernel = -com_unit &
            * radial_kernel_phase(ell, kr_response, kr_source, plasma_in%ks(j), j) &
            * plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%vT(j) &
            / (sqrt(2.0_dp * pi) * plasma_in%spec(sp)%lambda_D(j)**2 * sol) * even_w
    end function collisionless_jrad_bparallel_harmonic_sp

    subroutine collisionless_radial_moments(plasma_in, sp, ell, kr_response, &
            kr_source, j, even_o, odd_o, even_w)
        use config_m, only: artificial_debye_case, collisionless_kpar_epsilon
        use constants_m, only: pi
        use Krook_kernel_plasma_prefacs_m, only: Krook_collisionless_kpar, &
            Krook_collisionless_kpar_magnitude
        use setup_m, only: omega

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell, j
        real(dp), intent(in) :: kr_response, kr_source
        complex(dp), intent(out) :: even_o, odd_o, even_w

        complex(dp) :: o0, o2, w0, w2, k_pole, response_z, zeta
        complex(dp) :: plasma_Z
        real(dp) :: k_abs, resonance, v_t

        call validate_collisionless_radial_inputs(plasma_in, sp, j, &
            collisionless_kpar_epsilon)
        if (artificial_debye_case /= 0) then
            error stop 'Collisionless radial-current ions require artificial_debye_case=0'
        end if
        call check_harmonic_is_available(plasma_in, sp, ell)
        call radial_flr_coefficients(ell, plasma_in%ks(j), kr_source, &
            plasma_in%ks(j), kr_response, plasma_in%spec(sp)%rho_L(j), &
            plasma_in%spec(sp)%A1(j), plasma_in%spec(sp)%A2(j), o0, o2, w0, w2)

        v_t = plasma_in%spec(sp)%vT(j)
        k_abs = Krook_collisionless_kpar_magnitude(plasma_in%kp(j), &
            collisionless_kpar_epsilon)
        k_pole = Krook_collisionless_kpar(plasma_in%kp(j), &
            collisionless_kpar_epsilon)
        resonance = plasma_in%om_E(j) &
            + real(ell, dp) * plasma_in%spec(sp)%omega_c(j) - omega
        zeta = cmplx(-resonance / (sqrt(2.0_dp) * k_abs * v_t), 0.0_dp, dp)
        response_z = plasma_Z(zeta)

        even_o = sqrt(pi) / k_abs * (o0 * response_z &
            + 2.0_dp * o2 * zeta * (1.0_dp + zeta * response_z))
        even_w = sqrt(pi) / k_abs * (w0 * response_z &
            + 2.0_dp * w2 * zeta * (1.0_dp + zeta * response_z))
        odd_o = sqrt(2.0_dp * pi) * v_t / k_pole &
            * (o0 * (1.0_dp + zeta * response_z) &
            + 2.0_dp * o2 * (0.5_dp + zeta**2 * (1.0_dp + zeta * response_z)))
    end subroutine collisionless_radial_moments

    logical function use_collisionless_ion(sp) result(use_collisionless)
        use config_m, only: ion_collision_model

        integer, intent(in) :: sp

        select case (trim(ion_collision_model))
        case ('FokkerPlanck')
            use_collisionless = .false.
        case ('collisionless')
            use_collisionless = sp >= 1
        case default
            error stop 'Periodic ion_collision_model must be FokkerPlanck or collisionless'
        end select
    end function use_collisionless_ion

    subroutine validate_collisionless_radial_inputs(plasma_in, sp, j, epsilon)
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: epsilon

        if (sp < 1 .or. sp >= plasma_in%n_species) then
            error stop 'Collisionless radial-current kernel requires an ion species index'
        end if
        if (j < 1 .or. j > plasma_in%grid_size) then
            error stop 'Collisionless radial-current kernel radial index is out of range'
        end if
        if (epsilon <= 0.0_dp) then
            error stop 'Collisionless radial-current kernel requires epsilon > 0'
        end if
        if (plasma_in%spec(sp)%vT(j) <= 0.0_dp) then
            error stop 'Collisionless radial-current kernel requires vT > 0'
        end if
        if (plasma_in%spec(sp)%lambda_D(j) <= 0.0_dp) then
            error stop 'Collisionless radial-current kernel requires lambda_D > 0'
        end if
        if (abs(plasma_in%spec(sp)%omega_c(j)) <= tiny(1.0_dp)) then
            error stop 'Collisionless radial-current kernel requires nonzero omega_c'
        end if
    end subroutine validate_collisionless_radial_inputs

    complex(dp) function radial_kernel_phase(ell, kr_response, kr_source, ks, j) &
        result(phase)
        use constants_m, only: com_unit, pi
        use grid_m, only: rg_grid

        integer, intent(in) :: ell, j
        real(dp), intent(in) :: kr_response, kr_source, ks
        real(dp) :: alpha_response, alpha_source

        alpha_response = atan2(kr_response, ks)
        alpha_source = atan2(kr_source, ks)
        phase = exp(-com_unit * (kr_response - kr_source) * rg_grid%xb(j) &
            + com_unit * real(ell, dp) * (alpha_response - alpha_source)) &
            / (8.0_dp * pi**2)
    end function radial_kernel_phase

    integer function checked_harmonic_cutoff(cutoff) result(harmonic_cutoff)
        use setup_m, only: mphi_max

        integer, intent(in), optional :: cutoff

        harmonic_cutoff = mphi_max
        if (present(cutoff)) harmonic_cutoff = cutoff
        if (harmonic_cutoff < 0 .or. harmonic_cutoff > mphi_max) then
            error stop 'Radial-current harmonic cutoff must satisfy 0 <= cutoff <= mphi_max'
        end if
    end function checked_harmonic_cutoff

    subroutine check_harmonic_is_available(plasma_in, sp, ell)
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, ell

        if (ell < lbound(plasma_in%spec(sp)%I00, 2) &
            .or. ell > ubound(plasma_in%spec(sp)%I00, 2)) then
            error stop 'Requested radial-current harmonic is unavailable'
        end if
    end subroutine check_harmonic_is_available

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
