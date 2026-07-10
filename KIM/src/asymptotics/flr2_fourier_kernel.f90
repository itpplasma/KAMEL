module flr2_fourier_kernel_m
    ! Fused Fokker-Planck off-diagonal Fourier-space kernel integrand
    ! G(k_r, k'_r, r_g) for the forced-periodicity solver. Collapses to the
    ! diagonal source-of-truth calc_hatK_Phi_in_Fourier at k_r = k'_r.
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    implicit none
    private
    public :: hatG_rho_phi, hatG_rho_B
    public :: hatG_rho_phi_diag_sp, hatG_rho_B_diag_sp

    ! FLR-parameter convention for the two-wavenumber kernel. When .false.
    ! (default) the perpendicular wavenumber k_s is dropped from b_+ / b_x, so
    ! the kernel depends only on the radial wavenumbers k_r, k'_r. The .true.
    ! path (which reintroduces k_s^2) is implemented in a later task.
    logical, save, public :: kern_include_ks2 = .false.
contains
    complex(dp) function hatG_rho_phi(plasma_in, kr, krp, j) result(G)
        ! Fused two-wavenumber rho-Phi kernel Ghat^{rho Phi}(k_r, k'_r, r_g):
        ! sum the per-species core over the non-turned-off species, then apply
        ! the Fourier phase exp(i (k_r - k'_r) r_g) and the /(4*pi) source-of-
        ! truth normalization. On the diagonal k_r = k'_r the FLR pair collapses
        ! to a single b and the phase is unity, so this equals the diagonal
        ! integrand summed over species and divided by 4*pi (collapse by
        ! construction, since the diagonal shares this same core).
        use grid_m, only: rg_grid
        use constants_m, only: pi, com_unit
        use config_m, only: turn_off_ions, turn_off_electrons
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp
        real(dp) :: rho_L2, bplus, bcross
        complex(dp) :: acc

        acc = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            if (kern_include_ks2) then
                error stop 'hatG_rho_phi: kern_include_ks2=.true. path not implemented yet (Task 3b)'
            else
                rho_L2 = plasma_in%spec(sp)%rho_L(j)**2.0d0
                bplus = rho_L2 * (kr**2.0d0 + krp**2.0d0) / 2.0d0
                bcross = rho_L2 * abs(kr * krp)
            end if

            acc = acc + core_rho_phi_sp(plasma_in, sp, bplus, bcross, j)
        end do

        G = exp(com_unit * (kr - krp) * rg_grid%xb(j)) / (4.0d0 * pi) * acc
    end function hatG_rho_phi

    complex(dp) function hatG_rho_B(plasma_in, kr, krp, j) result(G)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        G = (0.0_dp, 0.0_dp)   ! stub — replaced in a later task
    end function hatG_rho_B

    complex(dp) function core_rho_phi_sp(plasma_in, sp, bplus, bcross, j) result(G)
        ! Per-species rho-Phi core (Debye + FP) parameterized by the FLR
        ! argument pair (b_+, b_x). This is the single home for the rho-Phi
        ! integrand math: the diagonal integrand calls it with bplus=bcross=b,
        ! and the fused two-wavenumber kernel calls it with the off-diagonal
        ! pair. b_+ scales exp(-b) and the (1 - b) polynomial factor; b_x is the
        ! argument (and the leading multiplier) of the modified Bessel terms.
        ! No phase, NO /(4*pi); species selection (turn_off_*) stays in callers.
        use config_m, only: artificial_debye_case
        use constants_m, only: com_unit
        use fortnum_special, only: bessel_in
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: bplus, bcross
        real(dp) :: ks

        G = (0.0_dp, 0.0_dp)
        ks = plasma_in%ks(j)

        if (artificial_debye_case <= 1) then
            G = - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0
        end if

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            G = G + 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 &
                * com_unit * plasma_in%spec(sp)%vT(j)**2.0d0 * ks &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j)) * exp(-bplus) * &
                (&
                    plasma_in%spec(sp)%I00(j, 0) * (&
                        bessel_in(0, bcross) * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-bplus)) &
                        + plasma_in%spec(sp)%A2(j) * bcross * bessel_in(-1, bcross) &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I20(j, 0) * plasma_in%spec(sp)%A2(j) * bessel_in(0, bcross) &
                )
        end if
    end function core_rho_phi_sp

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
        real(dp) :: b

        b = (kr**2.0d0) * plasma_in%spec(sp)%rho_L(j)**2.0d0
        G = core_rho_phi_sp(plasma_in, sp, b, b, j)
    end function hatG_rho_phi_diag_sp

    complex(dp) function hatG_rho_B_diag_sp(plasma_in, sp, kr, j) result(G)
        ! Per-species DIAGONAL rho-B integrand contribution: returns the
        ! SIGNED FP term (i.e. including the leading minus, so callers
        ! accumulate with '+'). No phase, NO /(4*pi). Species selection
        ! (turn_off_*) stays in the caller.
        use config_m, only: artificial_debye_case
        use constants_m, only: sol
        use fortnum_special, only: bessel_in
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr
        real(dp) :: b

        G = (0.0_dp, 0.0_dp)
        b = (kr**2.0d0) * plasma_in%spec(sp)%rho_L(j)**2.0d0

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            G = G - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 * plasma_in%spec(sp)%vT(j)**3.0d0 &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j) * sol) * exp(-b) * &
                (&
                    plasma_in%spec(sp)%I01(j, 0) * (&
                        bessel_in(0, b) * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-b)) &
                        + plasma_in%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I21(j, 0) * plasma_in%spec(sp)%A2(j) * bessel_in(0, b) &
                )
        end if
    end function hatG_rho_B_diag_sp
end module flr2_fourier_kernel_m
