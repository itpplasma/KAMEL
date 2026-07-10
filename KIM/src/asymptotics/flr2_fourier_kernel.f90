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
contains
    complex(dp) function hatG_rho_phi(plasma_in, kr, krp, j) result(G)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        G = (0.0_dp, 0.0_dp)   ! stub — replaced in a later task
    end function hatG_rho_phi

    complex(dp) function hatG_rho_B(plasma_in, kr, krp, j) result(G)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        G = (0.0_dp, 0.0_dp)   ! stub — replaced in a later task
    end function hatG_rho_B

    complex(dp) function hatG_rho_phi_diag_sp(plasma_in, sp, kr, j) result(G)
        ! Per-species DIAGONAL (single-wavenumber) rho-Phi integrand
        ! contribution: returns Debye term + FP term = the full per-species
        ! contribution to kernel_phi_temp. No phase, NO /(4*pi) — the caller
        ! applies those. Species selection (turn_off_*) stays in the caller.
        use config_m, only: artificial_debye_case
        use constants_m, only: com_unit
        use fortnum_special, only: bessel_in
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr
        real(dp) :: b

        G = (0.0_dp, 0.0_dp)
        b = (kr**2.0d0) * plasma_in%spec(sp)%rho_L(j)**2.0d0

        if (artificial_debye_case <= 1) then
            G = - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0
        end if

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            G = G + 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 &
                * com_unit * plasma_in%spec(sp)%vT(j)**2.0d0 * plasma_in%ks(j) &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j)) * exp(-b) * &
                (&
                    plasma_in%spec(sp)%I00(j, 0) * (&
                        bessel_in(0, b) * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-b)) &
                        + plasma_in%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I20(j, 0) * plasma_in%spec(sp)%A2(j) * bessel_in(0, b) &
                )
        end if
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
