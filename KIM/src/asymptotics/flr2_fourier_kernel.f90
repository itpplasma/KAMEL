module flr2_fourier_kernel_m
    ! Fused Fokker-Planck off-diagonal Fourier-space kernel integrand
    ! G(k_r, k'_r, r_g) for the forced-periodicity solver. Collapses to the
    ! diagonal source-of-truth calc_hatK_Phi_in_Fourier at k_r = k'_r.
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    implicit none
    private
    public :: hatG_rho_phi, hatG_rho_B
    public :: hatG_j_phi, hatG_j_B
    public :: hatG_rho_phi_diag_sp, hatG_rho_B_diag_sp
    public :: core_rho_phi_sp, core_rho_B_sp
    public :: core_j_phi_sp, core_j_B_sp
    public :: scaled_bessel_pair
    public :: flr_arg_pair_sp
    public :: set_global_kernel_approximations

    ! FLR-parameter convention for the two-wavenumber kernel. When .false.
    ! the perpendicular wavenumber k_s is dropped from b_+ / b_x, so
    ! the kernel depends only on the radial wavenumbers k_r, k'_r and collapses
    ! exactly to the diagonal source-of-truth calc_hatK_Phi_in_Fourier. The
    ! optional full form reintroduces k_s^2 via kperp = sqrt(k_s^2 + k_r^2).
    ! The global FEM kernel drops this contribution before its analytical
    ! inverse Fourier transform. The full periodic kernel retains it by default.
    logical, save, public :: kern_include_ks2 = .true.

    ! The global FEM implementation evaluates electron kernels in the zero-FLR
    ! limit but retains ion FLR. Keep the same approximation in the periodic
    ! solver only when the optional global-comparison approximation is enabled.
    logical, save, public :: kern_zero_flr_electrons = .false.
contains
    subroutine set_global_kernel_approximations(enabled)
        ! Select between the full periodic Fourier kernel and the simplified
        ! kernel used by the global FEM implementation.
        logical, intent(in) :: enabled

        kern_include_ks2 = .not. enabled
        kern_zero_flr_electrons = enabled
    end subroutine set_global_kernel_approximations

    subroutine flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
        ! Shared FLR argument pair (b_+, b_x) for the two-wavenumber kernels.
        ! Single home for the kern_include_ks2 branch used by both hatG_rho_phi
        ! and hatG_rho_B. When .false. the perpendicular wavenumber
        ! k_s is dropped, so the pair depends only on the radial wavenumbers
        ! k_r, k'_r and collapses exactly to the diagonal source-of-truth. The
        ! .true. path reintroduces k_s^2 via kperp = sqrt(k_s^2 + k_r^2), giving
        ! b_+ = rho_L^2/2 (kperp^2 + kperpp^2) and b_x = rho_L^2 kperp kperpp.
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr, krp
        real(dp), intent(out) :: bplus, bcross
        real(dp) :: rho_L2, ks

        if (kern_zero_flr_electrons .and. sp == 0) then
            bplus = 0.0_dp
            bcross = 0.0_dp
            return
        end if

        rho_L2 = plasma_in%spec(sp)%rho_L(j)**2.0d0
        if (kern_include_ks2) then
            ks = plasma_in%ks(j)
            bplus = rho_L2 * (2.0d0 * ks**2 + kr**2 + krp**2) / 2.0d0
            bcross = rho_L2 * sqrt(ks**2 + kr**2) * sqrt(ks**2 + krp**2)
        else
            bplus = rho_L2 * (kr**2 + krp**2) / 2.0d0
            bcross = rho_L2 * abs(kr * krp)
        end if
    end subroutine flr_arg_pair_sp

    subroutine scaled_bessel_pair(bplus, bcross, sI0, sIm1)
        ! Overflow-safe scaled modified Bessel products:
        !   sI0  = exp(-bplus) * I_0(bcross)
        !   sIm1 = exp(-bplus) * I_{-1}(bcross)
        ! For large bcross the unscaled I_n(bcross) would overflow before the
        ! exp(-bplus) damping (b_+ >= b_x by construction), so use the large-
        ! argument asymptotics (lifted from the recovered kernel_mod.f90).
        use fortnum_special, only: bessel_in
        real(dp), intent(in) :: bplus, bcross
        real(dp), intent(out) :: sI0, sIm1
        real(dp), parameter :: asymptotic_threshold = 500.0d0  ! safely below the ~710 I_0 overflow
        real(dp), parameter :: twopi = 6.283185307179586d0

        if (bcross > asymptotic_threshold) then
            sI0  = exp(bcross - bplus) / sqrt(twopi * bcross)
            sIm1 = exp(-bplus + asinh(-1.0d0/bcross) + bcross*sqrt(1.0d0 + 1.0d0/bcross**2)) &
                   / sqrt(twopi * bcross * sqrt(1.0d0 + 1.0d0/bcross**2))
        else
            sI0  = exp(-bplus) * bessel_in(0, bcross)
            sIm1 = exp(-bplus) * bessel_in(-1, bcross)
        end if
    end subroutine scaled_bessel_pair

    complex(dp) function hatG_rho_phi(plasma_in, kr, krp, j) result(G)
        ! Fused two-wavenumber rho-Phi kernel Ghat^{rho Phi}(k_r, k'_r, r_g):
        ! sum the per-species core over the non-turned-off species, then apply
        ! the Fourier phase exp(-i (k_r - k'_r) r_g) and the continuous-kernel
        ! normalization 1/(8*pi^2). The latter follows directly from the
        ! Fourier definition K_kk'=(1/2*pi) int dr dr' exp(-ikr) K(r,r')
        ! exp(ik'r'): a local real-space core/(4*pi) becomes core/(8*pi^2).
        use grid_m, only: rg_grid
        use constants_m, only: pi, com_unit
        use config_m, only: turn_off_ions, turn_off_electrons
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp
        real(dp) :: bplus, bcross
        complex(dp) :: acc

        acc = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            call flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
            acc = acc + core_rho_phi_sp(plasma_in, sp, bplus, bcross, j)
        end do

        G = exp(-com_unit * (kr - krp) * rg_grid%xb(j)) / (8.0_dp * pi**2) * acc
    end function hatG_rho_phi

    complex(dp) function hatG_rho_B(plasma_in, kr, krp, j) result(G)
        ! Fused two-wavenumber rho-B kernel Ghat^{rho B}(k_r, k'_r, r_g): sum the
        ! signed per-species core over the non-turned-off species, then apply the
        ! Fourier phase exp(-i (k_r - k'_r) r_g) and the same 1/(8*pi^2)
        ! continuous-kernel normalization as rho-Phi.
        ! rho-B has no Debye term (handled inside core_rho_B_sp).
        use grid_m, only: rg_grid
        use constants_m, only: pi, com_unit
        use config_m, only: turn_off_ions, turn_off_electrons
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp
        real(dp) :: bplus, bcross
        complex(dp) :: acc

        acc = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            call flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
            acc = acc + core_rho_B_sp(plasma_in, sp, bplus, bcross, j)
        end do

        G = exp(-com_unit * (kr - krp) * rg_grid%xb(j)) / (8.0_dp * pi**2) * acc
    end function hatG_rho_B

    complex(dp) function hatG_j_phi(plasma_in, kr, krp, j) result(G)
        ! Fused two-wavenumber j-Phi kernel Ghat^{jPhi}(k_r, k'_r, r_g): sum the
        ! per-species core over the non-turned-off species, then apply the
        ! Fourier phase exp(-i (k_r - k'_r) r_g). On the diagonal k_r = k'_r the
        ! FLR pair collapses to a single b and the phase is unity, so this equals
        ! the diagonal j-Phi expression of calc_hatK_Phi_in_Fourier summed over
        ! species (thesis 14.5).
        !
        ! Equations (14.5) and (14.6) give the current kernels the Fourier
        ! normalization 1/(2^3*pi^2). This factor belongs to the kinetic kernel,
        ! independently of the 4*pi appearing in Poisson's equation.
        use grid_m, only: rg_grid
        use constants_m, only: pi, com_unit
        use config_m, only: turn_off_ions, turn_off_electrons
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp
        real(dp) :: bplus, bcross
        complex(dp) :: acc

        acc = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            call flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
            acc = acc + core_j_phi_sp(plasma_in, sp, bplus, bcross, j)
        end do

        G = exp(-com_unit * (kr - krp) * rg_grid%xb(j)) &
            / (8.0_dp * pi**2) * acc
    end function hatG_j_phi

    complex(dp) function hatG_j_B(plasma_in, kr, krp, j) result(G)
        ! Fused two-wavenumber j-B kernel Ghat^{jB}(k_r, k'_r, r_g): sum the
        ! signed per-species core over the non-turned-off species, then apply the
        ! Fourier phase exp(-i (k_r - k'_r) r_g). On the diagonal k_r = k'_r this
        ! equals the diagonal j-B expression of calc_hatK_Phi_in_Fourier summed
        ! over species. Direct output/source moment counting requires I11 and I13;
        ! the I02 and I22 indices printed in thesis equation 14.6 are incorrect.
        !
        ! The same 1/(8*pi^2) normalization applies as in hatG_j_phi.
        use grid_m, only: rg_grid
        use constants_m, only: pi, com_unit
        use config_m, only: turn_off_ions, turn_off_electrons
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        integer :: sp
        real(dp) :: bplus, bcross
        complex(dp) :: acc

        acc = (0.0_dp, 0.0_dp)
        do sp = 0, plasma_in%n_species - 1
            if (turn_off_ions .and. sp >= 1) cycle
            if (turn_off_electrons .and. sp == 0) cycle

            call flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
            acc = acc + core_j_B_sp(plasma_in, sp, bplus, bcross, j)
        end do

        G = exp(-com_unit * (kr - krp) * rg_grid%xb(j)) &
            / (8.0_dp * pi**2) * acc
    end function hatG_j_B

    complex(dp) function core_j_phi_sp(plasma_in, sp, bplus, bcross, j) result(G)
        ! Per-species j-Phi core (FP only, NO Debye term) parameterized by the
        ! FLR argument pair (b_+, b_x). Shares the skeleton of core_rho_phi_sp;
        ! only the moment pair and the prefactor differ.
        !
        ! Moments (I01, I21) per thesis (14.5). j-Phi and rho-B legitimately
        ! SHARE moments: the pair is (I^{0s}, I^{2s}) with s = moment power +
        ! drive power, and both have s = 1. This is not a copy-paste artifact.
        !
        ! The prefactor carries vT^3 = vT^(2+s), not the vT^2 printed in (14.5):
        ! (14.3)-(14.6) all print vT^2, and (14.4)'s was confirmed dimensionally
        ! to be a misprint. The code's powers follow vT^(2+s) throughout.
        !
        ! No phase or Fourier normalization; callers supply both after summing
        ! species contributions.
        use config_m, only: artificial_debye_case
        use constants_m, only: com_unit
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: bplus, bcross
        real(dp) :: ks, sI0, sIm1

        G = (0.0_dp, 0.0_dp)
        ks = plasma_in%ks(j)

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            call scaled_bessel_pair(bplus, bcross, sI0, sIm1)
            G = 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 &
                * com_unit * plasma_in%spec(sp)%vT(j)**3.0d0 * ks &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j)) * &
                (&
                    plasma_in%spec(sp)%I01(j, 0) * (&
                        sI0 * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-bplus)) &
                        + plasma_in%spec(sp)%A2(j) * bcross * sIm1 &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I21(j, 0) * plasma_in%spec(sp)%A2(j) * sI0 &
                )
        end if
    end function core_j_phi_sp

    complex(dp) function core_j_B_sp(plasma_in, sp, bplus, bcross, j) result(G)
        ! Per-species j-B core (SIGNED FP term, NO Debye) parameterized by the
        ! FLR argument pair (b_+, b_x). Carries the leading minus, so callers
        ! accumulate with '+'.
        !
        ! Direct output/source power counting gives moments I11 and I13: the
        ! current output supplies one power of u_par, while the B drive supplies
        ! one or three source powers. The thesis's I02 and I22 indices are
        ! incorrect. With s = 2, the prefactor is vT^(2+s) = vT^4, not the vT^2
        ! printed in thesis equation 14.6. See core_j_phi_sp.
        !
        ! No phase or Fourier normalization; callers supply both after summing
        ! species contributions.
        use config_m, only: artificial_debye_case
        use constants_m, only: sol
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: bplus, bcross
        real(dp) :: sI0, sIm1

        G = (0.0_dp, 0.0_dp)

        if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
            call scaled_bessel_pair(bplus, bcross, sI0, sIm1)
            G = - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 &
                * plasma_in%spec(sp)%vT(j)**4.0d0 &
                / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j) * sol) * &
                (&
                    plasma_in%spec(sp)%I11(j, 0) * (&
                        sI0 * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-bplus)) &
                        + plasma_in%spec(sp)%A2(j) * bcross * sIm1 &
                    )&
                    + 0.5d0 * plasma_in%spec(sp)%I13(j, 0) * plasma_in%spec(sp)%A2(j) * sI0 &
                )
        end if
    end function core_j_B_sp

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
        ! Per-species rho-B core (SIGNED FP term) parameterized by the FLR
        ! argument pair (b_+, b_x). Single home for the rho-B integrand math:
        ! the diagonal integrand calls it with bplus=bcross=b, and the fused
        ! two-wavenumber kernel calls it with the off-diagonal pair. b_+ scales
        ! exp(-b) and the (1 - b) polynomial factor; b_x is the argument (and the
        ! leading multiplier) of the modified Bessel terms. rho-B carries the
        ! leading minus, so callers accumulate with '+'. NO Debye term, no phase,
        ! NO /(4*pi); species selection (turn_off_*) stays in callers.
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
        real(dp) :: b

        b = (kr**2.0d0) * plasma_in%spec(sp)%rho_L(j)**2.0d0
        G = core_rho_phi_sp(plasma_in, sp, b, b, j)
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
        real(dp) :: b

        b = (kr**2.0d0) * plasma_in%spec(sp)%rho_L(j)**2.0d0
        G = core_rho_B_sp(plasma_in, sp, b, b, j)
    end function hatG_rho_B_diag_sp
end module flr2_fourier_kernel_m
