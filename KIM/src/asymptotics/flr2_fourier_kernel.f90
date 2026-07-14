module flr2_fourier_kernel_m
    ! Oracle-verified Fokker-Planck off-diagonal Fourier-space kernel integrand
    ! G(k_r, k'_r, r_g) for the forced-periodicity solver. Production uses the
    ! full finite-radius k_perp model; the historical radial-only diagonal is
    ! retained as an explicitly named approximation for compatibility checks.
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    implicit none
    private
    public :: hatG_rho_phi, hatG_rho_B
    public :: hatG_rho_phi_diag_sp, hatG_rho_B_diag_sp
    public :: core_rho_phi_sp, core_rho_B_sp
    public :: scaled_bessel_pair, flr_arg_pair
    public :: fp_rho_phi_harmonic, fp_rho_B_harmonic, fp_rho_phi_debye
    public :: fp_ks_model_from_name, fp_mphi_truncation_check

    integer, parameter, public :: FP_KPERP_FULL = 1
    integer, parameter, public :: FP_KR_ONLY_APPROX = 2
    ! Selected via the KIM_CONFIG namelist key fp_ks_model_name
    ! ('kperp_full' | 'kr_only'); kim_read_config applies the mapping.
    integer, save, public :: fp_ks_model = FP_KPERP_FULL
contains
    pure subroutine fp_ks_model_from_name(name, model, status)
        ! Map the config-facing model name onto the module constants.
        ! status /= 0 marks an unknown name; model then stays at the default.
        character(len=*), intent(in) :: name
        integer, intent(out) :: model, status

        status = 0
        select case (trim(name))
        case ('kperp_full')
            model = FP_KPERP_FULL
        case ('kr_only')
            model = FP_KR_ONLY_APPROX
        case default
            model = FP_KPERP_FULL
            status = 1
        end select
    end subroutine fp_ks_model_from_name

    pure subroutine fp_mphi_truncation_check(harmonics, tolerance, tail_fraction, ok)
        ! Validation gate for the cyclotron-harmonic truncation. harmonics(:)
        ! are the per-m_phi kernel contributions ordered m_phi = -mphi_max ..
        ! +mphi_max (size 2*mphi_max+1). The estimated truncation tail is the
        ! outermost-shell magnitude relative to the total. Sums too short to
        ! carry a measurable tail (mphi_max = 0), even-sized inputs, and a
        ! vanishing total fail closed: hosts must not treat an unassessable
        ! truncation as verified.
        complex(dp), intent(in) :: harmonics(:)
        real(dp), intent(in) :: tolerance
        real(dp), intent(out) :: tail_fraction
        logical, intent(out) :: ok
        complex(dp) :: total
        integer :: n

        n = size(harmonics)
        total = sum(harmonics)
        if (n < 3 .or. mod(n, 2) == 0 .or. abs(total) == 0.0_dp) then
            tail_fraction = huge(1.0_dp)
            ok = .false.
            return
        end if
        tail_fraction = (abs(harmonics(1)) + abs(harmonics(n)))/abs(total)
        ok = tail_fraction < tolerance
    end subroutine fp_mphi_truncation_check

    pure subroutine flr_arg_pair(ks, kr, krp, rho_L, model, bplus, bcross, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(dp), intent(in) :: ks, kr, krp, rho_L
        integer, intent(in) :: model
        real(dp), intent(out) :: bplus, bcross
        integer, intent(out) :: status
        real(dp) :: rho_L2

        status = 0
        bplus = 0.0_dp
        bcross = 0.0_dp
        if (.not. ieee_is_finite(kr) .or. .not. ieee_is_finite(krp) .or. &
            .not. ieee_is_finite(rho_L) .or. rho_L < 0.0_dp) then
            status = 1
            return
        end if

        rho_L2 = rho_L*rho_L
        select case (model)
        case (FP_KPERP_FULL)
            if (.not. ieee_is_finite(ks)) then
                status = 2
                return
            end if
            bplus = rho_L2*(2.0_dp*ks*ks + kr*kr + krp*krp)/2.0_dp
            bcross = rho_L2*sqrt(ks*ks + kr*kr)*sqrt(ks*ks + krp*krp)
        case (FP_KR_ONLY_APPROX)
            bplus = rho_L2*(kr*kr + krp*krp)/2.0_dp
            bcross = rho_L2*abs(kr*krp)
        case default
            status = 3
        end select
    end subroutine flr_arg_pair

    subroutine flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr, krp
        real(dp), intent(out) :: bplus, bcross
        integer :: status

        call flr_arg_pair(plasma_in%ks(j), kr, krp, plasma_in%spec(sp)%rho_L(j), &
            fp_ks_model, bplus, bcross, status)
        if (status /= 0) error stop 'invalid finite-radius FLR arguments'
    end subroutine flr_arg_pair_sp

    subroutine scaled_bessel_harmonic(bplus, bcross, order, value, status)
        use fortnum_special, only: bessel_in
        real(dp), intent(in) :: bplus, bcross
        integer, intent(in) :: order
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp), parameter :: asymptotic_threshold = 500.0_dp
        real(dp), parameter :: twopi = 6.283185307179586_dp
        real(dp) :: mu, term, series
        integer :: k, n

        status = 0
        value = 0.0_dp
        if (bcross < 0.0_dp .or. bplus < 0.0_dp .or. &
            bplus + 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, bplus) < bcross) then
            status = 1
            return
        end if
        n = abs(order)
        if (bcross == 0.0_dp) then
            if (n == 0) value = exp(-bplus)
        else if (bcross > asymptotic_threshold) then
            mu = 4.0_dp*real(n*n, dp)
            series = 1.0_dp
            term = 1.0_dp
            do k = 1, 5
                term = term * (-(mu - real((2*k - 1)**2, dp))) / &
                    (real(k, dp)*8.0_dp*bcross)
                series = series + term
            end do
            value = exp(bcross - bplus)*series/sqrt(twopi*bcross)
        else
            value = exp(-bplus)*bessel_in(n, bcross)
        end if
    end subroutine scaled_bessel_harmonic

    subroutine scaled_bessel_pair(bplus, bcross, sI0, sIm1)
        real(dp), intent(in) :: bplus, bcross
        real(dp), intent(out) :: sI0, sIm1
        integer :: status

        call scaled_bessel_harmonic(bplus, bcross, 0, sI0, status)
        if (status /= 0) error stop 'invalid scaled Bessel arguments'
        call scaled_bessel_harmonic(bplus, bcross, -1, sIm1, status)
        if (status /= 0) error stop 'invalid scaled Bessel arguments'
    end subroutine scaled_bessel_pair

    subroutine fp_rho_phi_harmonic(lambda_D, vT, omega_c, nu, rho_L, ks, kr, krp, rg, &
            mphi, A1, A2, I00, I20, model, value, status)
        ! rho_L is an explicit input (not recomputed from vT/omega_c) so the
        ! caller's stored species value — including the ion_flr_scale_factor
        ! config option — reaches the FLR arguments, as in flr_arg_pair_sp.
        use constants_m, only: pi, com_unit
        real(dp), intent(in) :: lambda_D, vT, omega_c, nu, rho_L, ks, kr, krp, rg
        integer, intent(in) :: mphi, model
        real(dp), intent(in) :: A1, A2
        complex(dp), intent(in) :: I00, I20
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp) :: bplus, bcross, sIm, sImm1, energy_moment, alpha, alphap
        complex(dp) :: phase

        value = (0.0_dp, 0.0_dp)
        if (lambda_D <= 0.0_dp .or. omega_c == 0.0_dp .or. nu <= 0.0_dp) then
            status = 4
            return
        end if
        call flr_arg_pair(ks, kr, krp, rho_L, model, bplus, bcross, status)
        if (status /= 0) return
        call scaled_bessel_harmonic(bplus, bcross, mphi, sIm, status)
        if (status /= 0) return
        call scaled_bessel_harmonic(bplus, bcross, mphi - 1, sImm1, status)
        if (status /= 0) return

        energy_moment = (1.0_dp - bplus - real(mphi, dp))*sIm + bcross*sImm1
        alpha = atan2(kr, ks)
        alphap = atan2(krp, ks)
        phase = exp(com_unit*((kr - krp)*rg + real(mphi, dp)*(alphap - alpha)))
        value = phase/(4.0_dp*pi*lambda_D**2) * com_unit*vT**2*ks/(omega_c*nu) * &
            (I00*(A1*sIm + A2*energy_moment) + 0.5_dp*A2*I20*sIm)
    end subroutine fp_rho_phi_harmonic

    subroutine fp_rho_B_harmonic(lambda_D, vT, omega_c, nu, rho_L, ks, kr, krp, rg, &
            mphi, A1, A2, I01, I21, model, value, status)
        ! rho_L is an explicit input for the same reason as in
        ! fp_rho_phi_harmonic.
        use constants_m, only: pi, com_unit, sol
        real(dp), intent(in) :: lambda_D, vT, omega_c, nu, rho_L, ks, kr, krp, rg
        integer, intent(in) :: mphi, model
        real(dp), intent(in) :: A1, A2
        complex(dp), intent(in) :: I01, I21
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp) :: bplus, bcross, sIm, sImm1, energy_moment, alpha, alphap
        complex(dp) :: phase

        value = (0.0_dp, 0.0_dp)
        if (lambda_D <= 0.0_dp .or. omega_c == 0.0_dp .or. nu <= 0.0_dp) then
            status = 4
            return
        end if
        call flr_arg_pair(ks, kr, krp, rho_L, model, bplus, bcross, status)
        if (status /= 0) return
        call scaled_bessel_harmonic(bplus, bcross, mphi, sIm, status)
        if (status /= 0) return
        call scaled_bessel_harmonic(bplus, bcross, mphi - 1, sImm1, status)
        if (status /= 0) return

        energy_moment = (1.0_dp - bplus - real(mphi, dp))*sIm + bcross*sImm1
        alpha = atan2(kr, ks)
        alphap = atan2(krp, ks)
        phase = exp(com_unit*((kr - krp)*rg + real(mphi, dp)*(alphap - alpha)))
        value = -phase/(4.0_dp*pi*lambda_D**2) * vT**3/(omega_c*nu*sol) * &
            (I01*(A1*sIm + A2*energy_moment) + 0.5_dp*A2*I21*sIm)
    end subroutine fp_rho_B_harmonic

    pure subroutine fp_rho_phi_debye(lambda_D, kr, krp, rg, include_debye, value)
        use constants_m, only: pi, com_unit
        real(dp), intent(in) :: lambda_D, kr, krp, rg
        logical, intent(in) :: include_debye
        complex(dp), intent(out) :: value

        value = (0.0_dp, 0.0_dp)
        if (include_debye) then
            value = -exp(com_unit*(kr - krp)*rg)/(4.0_dp*pi*lambda_D**2)
        end if
    end subroutine fp_rho_phi_debye

    complex(dp) function hatG_rho_phi(plasma_in, kr, krp, j) result(G)
        ! Full two-wavenumber rho-Phi kernel Ghat^{rho Phi}(k_r,k'_r,r_g).
        ! Each retained harmonic carries its gyro-angle/radial Fourier phase,
        ! finite-radius moments, and 1/(4*pi) normalization. The Debye term is
        ! independent of m_phi and is therefore included exactly once/species.
        use grid_m, only: rg_grid
        use config_m, only: turn_off_ions, turn_off_electrons
        use config_m, only: artificial_debye_case
        use setup_m, only: mphi_max
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp, mphi, status
        complex(dp) :: term

        G = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            call fp_rho_phi_debye(plasma_in%spec(sp)%lambda_D(j), kr, krp, &
                rg_grid%xb(j), artificial_debye_case <= 1, term)
            G = G + term
            if (artificial_debye_case == 1) cycle

            do mphi = -mphi_max, mphi_max
                call fp_rho_phi_harmonic(plasma_in%spec(sp)%lambda_D(j), &
                    plasma_in%spec(sp)%vT(j), plasma_in%spec(sp)%omega_c(j), &
                    plasma_in%spec(sp)%nu(j), plasma_in%spec(sp)%rho_L(j), &
                    plasma_in%ks(j), kr, krp, rg_grid%xb(j), &
                    mphi, plasma_in%spec(sp)%A1(j), plasma_in%spec(sp)%A2(j), &
                    plasma_in%spec(sp)%I00(j, mphi), plasma_in%spec(sp)%I20(j, mphi), &
                    fp_ks_model, term, status)
                if (status /= 0) error stop 'invalid rho-Phi harmonic arguments'
                G = G + term
            end do
        end do
    end function hatG_rho_phi

    complex(dp) function hatG_rho_B(plasma_in, kr, krp, j) result(G)
        ! Full two-wavenumber rho-B kernel Ghat^{rho B}(k_r,k'_r,r_g).
        ! Each retained harmonic carries its gyro-angle/radial Fourier phase,
        ! signed current moment, finite-radius factors, and 1/(4*pi)
        ! normalization. rho-B has no Debye contribution.
        use grid_m, only: rg_grid
        use config_m, only: turn_off_ions, turn_off_electrons
        use config_m, only: artificial_debye_case
        use setup_m, only: mphi_max
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp, mphi, status
        complex(dp) :: term

        G = (0.0_dp, 0.0_dp)
        if (artificial_debye_case == 1) return
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            do mphi = -mphi_max, mphi_max
                call fp_rho_B_harmonic(plasma_in%spec(sp)%lambda_D(j), &
                    plasma_in%spec(sp)%vT(j), plasma_in%spec(sp)%omega_c(j), &
                    plasma_in%spec(sp)%nu(j), plasma_in%spec(sp)%rho_L(j), &
                    plasma_in%ks(j), kr, krp, rg_grid%xb(j), &
                    mphi, plasma_in%spec(sp)%A1(j), plasma_in%spec(sp)%A2(j), &
                    plasma_in%spec(sp)%I01(j, mphi), plasma_in%spec(sp)%I21(j, mphi), &
                    fp_ks_model, term, status)
                if (status /= 0) error stop 'invalid rho-B harmonic arguments'
                G = G + term
            end do
        end do
    end function hatG_rho_B

    complex(dp) function core_rho_phi_sp(plasma_in, sp, bplus, bcross, j) result(G)
        ! Historical m_phi=0 per-species rho-Phi core used by diagonal
        ! compatibility tests and diagnostics. The production off-diagonal
        ! path uses fp_rho_phi_harmonic instead. No phase or 1/(4*pi).
        use config_m, only: artificial_debye_case
        use constants_m, only: com_unit
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: bplus, bcross
        real(dp) :: ks, sI0, sIm1

        G = (0.0_dp, 0.0_dp)
        ks = plasma_in%ks(j)

        if (artificial_debye_case <= 1) then
            G = - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0
        end if

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            call scaled_bessel_pair(bplus, bcross, sI0, sIm1)
            G = G + 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 &
                * com_unit * plasma_in%spec(sp)%vT(j)**2.0d0 * ks &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j)) * &
                (&
                    plasma_in%spec(sp)%I00(j, 0) * (&
                        sI0 * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-bplus)) &
                        + plasma_in%spec(sp)%A2(j) * bcross * sIm1 &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I20(j, 0) * plasma_in%spec(sp)%A2(j) * sI0 &
                )
        end if
    end function core_rho_phi_sp

    complex(dp) function core_rho_B_sp(plasma_in, sp, bplus, bcross, j) result(G)
        ! Historical m_phi=0 per-species signed rho-B core used by diagonal
        ! compatibility tests and diagnostics. The production off-diagonal
        ! path uses fp_rho_B_harmonic instead. No phase or 1/(4*pi).
        use config_m, only: artificial_debye_case
        use constants_m, only: sol
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: bplus, bcross
        real(dp) :: sI0, sIm1

        G = (0.0_dp, 0.0_dp)

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            call scaled_bessel_pair(bplus, bcross, sI0, sIm1)
            G = G - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 * plasma_in%spec(sp)%vT(j)**3.0d0 &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j) * sol) * &
                (&
                    plasma_in%spec(sp)%I01(j, 0) * (&
                        sI0 * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-bplus)) &
                        + plasma_in%spec(sp)%A2(j) * bcross * sIm1 &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I21(j, 0) * plasma_in%spec(sp)%A2(j) * sI0 &
                )
        end if
    end function core_rho_B_sp

    complex(dp) function hatG_rho_phi_diag_sp(plasma_in, sp, kr, j) result(G)
        ! Per-species DIAGONAL (single-wavenumber) rho-Phi integrand
        ! contribution: returns Debye term + FP term = the full per-species
        ! contribution to kernel_phi_temp. No phase, NO /(4*pi) — the caller
        ! applies those. Species selection (turn_off_*) stays in the caller.
        ! Delegates to core_rho_phi_sp with b_+ = b_x = b so the diagonal math
        ! is single-homed.
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr
        real(dp) :: b, bcross

        call flr_arg_pair_sp(plasma_in, sp, kr, kr, j, b, bcross)
        G = core_rho_phi_sp(plasma_in, sp, b, bcross, j)
    end function hatG_rho_phi_diag_sp

    complex(dp) function hatG_rho_B_diag_sp(plasma_in, sp, kr, j) result(G)
        ! Per-species DIAGONAL rho-B integrand contribution: returns the
        ! SIGNED FP term (i.e. including the leading minus, so callers
        ! accumulate with '+'). No phase, NO /(4*pi). Species selection
        ! (turn_off_*) stays in the caller. Delegates to core_rho_B_sp with
        ! b_+ = b_x = b so the diagonal math is single-homed.
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr
        real(dp) :: b, bcross

        call flr_arg_pair_sp(plasma_in, sp, kr, kr, j, b, bcross)
        G = core_rho_B_sp(plasma_in, sp, b, bcross, j)
    end function hatG_rho_B_diag_sp
end module flr2_fourier_kernel_m
