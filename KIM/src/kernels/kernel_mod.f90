module kernels_m
    ! Continuous-Fourier KIM kernels G(k_r, k_r', r_g) restored from
    ! kernel_mod.f90 at 1ae0aeeb~1 for the enforced-periodicity solver:
    ! the periodic matrix element is K_{m,m'} = (2 pi / L) times the
    ! one-period r_g integral of these functions (Kasilov, "Enforced
    ! periodicity", 2026-07-04). Fourier phase convention
    ! exp(i k_r r); CGS-Gaussian units; the rho_phi kernel carries
    ! 1/(8 pi^2), the rho_B kernel i/(8 pi^2 c).

    use KIM_kinds_m, only: dp

    implicit none

    ! Adiabatic-only response: drop the thermodynamic-force terms while
    ! keeping the sign, gyroaverage, and Fourier phase of the full
    ! expression. Matches the artificial_debye_case semantics of the
    ! hat-basis path; the historical branch flipped the sign and dropped
    ! the phase and is not restored.
    logical :: kernel_debye_case = .false.

    real(dp) :: bessel_large_arg_limit = 10d0

    integer, parameter :: nlagr = 4
    integer, parameter :: nder = 0

contains

    subroutine lagrange_weights(val_rg, ibeg, iend, weights)

        use species_m, only: plasma

        implicit none

        real(dp), intent(in) :: val_rg
        integer, intent(out) :: ibeg, iend
        real(dp), intent(out) :: weights(nlagr)

        real(dp) :: coef(0:nder, nlagr)
        integer :: ir

        call binsrc(plasma%r_grid, 1, plasma%grid_size, val_rg, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > plasma%grid_size) then
            iend = plasma%grid_size
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, val_rg, plasma%r_grid(ibeg:iend), coef)
        weights = coef(0, :)

    end subroutine lagrange_weights

    ! This is without the exp(i k_r(r_g - x_l)) factor
    complex(dp) function kernel_rho_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

        use constants_m, only: pi, com_unit
        use species_m, only: plasma
        use fortnum_special, only: bessel_in

        implicit none

        real(dp), intent(in) :: val_kr, val_krp, val_rg

        complex(dp), external :: plasma_Z

        real(dp) :: weights(nlagr)
        integer :: ibeg, iend, sp
        real(dp) :: eval_bp, eval_bt
        complex(dp) :: z0_interp, eval_besselI0, eval_besselIm1
        real(dp) :: vT_interp, omc_interp, ks_interp, kp_interp, A1_interp, &
                    A2_interp, lambda_D_interp, rhoL_interp, kperp, kperpp

        kernel_rho_phi_of_kr_krp_rg = (0.0d0, 0.0d0)

        call lagrange_weights(val_rg, ibeg, iend, weights)
        ks_interp = sum(weights*plasma%ks(ibeg:iend))
        kp_interp = sum(weights*plasma%kp(ibeg:iend))
        kperp = sqrt(ks_interp**2 + val_kr**2)
        kperpp = sqrt(ks_interp**2 + val_krp**2)

        do sp = 0, plasma%n_species - 1
            associate (spec => plasma%spec(sp))
                vT_interp = sum(weights*spec%vT(ibeg:iend))
                omc_interp = sum(weights*spec%omega_c(ibeg:iend))
                A1_interp = sum(weights*spec%A1(ibeg:iend))
                A2_interp = sum(weights*spec%A2(ibeg:iend))
                lambda_D_interp = sum(weights*spec%lambda_D(ibeg:iend))
                z0_interp = sum(weights*spec%z0(ibeg:iend))
            end associate

            rhoL_interp = vT_interp/abs(omc_interp)

            eval_bp = rhoL_interp**2/2.0d0*(kperp**2 + kperpp**2)
            eval_bt = rhoL_interp**2*kperp*kperpp

            if (kernel_debye_case) then
                A1_interp = 0.0d0
                A2_interp = 0.0d0
            end if

            if (eval_bt > bessel_large_arg_limit) then
                ! limit close to magnetic axis (k_s -> infinity) and
                ! large k_r and k_rp: asymptotics for Bessel I functions
                eval_besselI0 = exp(eval_bt - eval_bp) &
                                /sqrt(2.0d0*pi*eval_bt)
                eval_besselIm1 = exp(-eval_bp + asinh(-1.0d0/eval_bt) &
                                     + eval_bt*sqrt(1.0d0 + 1.0d0/eval_bt**2)) &
                                 /sqrt(2.0d0*pi*eval_bt &
                                       *sqrt(1.0d0 + 1.0d0/eval_bt**2))
            else
                eval_besselI0 = bessel_in(0, eval_bt)*exp(-eval_bp)
                eval_besselIm1 = bessel_in(-1, eval_bt)*exp(-eval_bp)
            end if

            kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg &
                + 1.0d0/lambda_D_interp**2 &
                *exp(com_unit*(val_kr - val_krp)*val_rg) &
                *(-exp(-rhoL_interp**2/2.0d0*(val_kr - val_krp)**2) &
                  + ks_interp*rhoL_interp/(kp_interp*sqrt(2.0d0)) &
                  *(A1_interp*eval_besselI0*plasma_Z(z0_interp) &
                    + A2_interp*(plasma_Z(z0_interp)*eval_besselI0 &
                                 *(1.0d0 + eval_bp + z0_interp**2) &
                                 + eval_besselIm1*eval_bt &
                                 + z0_interp*eval_besselI0)))
        end do

        kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg &
                                      /(2.0d0**3*pi**2)

    end function kernel_rho_phi_of_kr_krp_rg

    complex(dp) function kernel_rho_B_of_kr_krp_rg(val_kr, val_krp, val_rg)

        use constants_m, only: pi, com_unit, sol
        use species_m, only: plasma
        use fortnum_special, only: bessel_in

        implicit none

        real(dp), intent(in) :: val_kr, val_krp, val_rg

        complex(dp), external :: plasma_Z

        real(dp) :: weights(nlagr)
        integer :: ibeg, iend, sp
        real(dp) :: eval_bp, eval_bt
        complex(dp) :: z0_interp, eval_besselI0, eval_besselIm1
        real(dp) :: vT_interp, omc_interp, ks_interp, kp_interp, A1_interp, &
                    A2_interp, lambda_D_interp, rhoL_interp, kperp, kperpp

        kernel_rho_B_of_kr_krp_rg = (0.0d0, 0.0d0)

        call lagrange_weights(val_rg, ibeg, iend, weights)
        ks_interp = sum(weights*plasma%ks(ibeg:iend))
        kp_interp = sum(weights*plasma%kp(ibeg:iend))
        kperp = sqrt(ks_interp**2 + val_kr**2)
        kperpp = sqrt(ks_interp**2 + val_krp**2)

        do sp = 0, plasma%n_species - 1
            associate (spec => plasma%spec(sp))
                vT_interp = sum(weights*spec%vT(ibeg:iend))
                omc_interp = sum(weights*spec%omega_c(ibeg:iend))
                A1_interp = sum(weights*spec%A1(ibeg:iend))
                A2_interp = sum(weights*spec%A2(ibeg:iend))
                lambda_D_interp = sum(weights*spec%lambda_D(ibeg:iend))
                z0_interp = sum(weights*spec%z0(ibeg:iend))
            end associate

            rhoL_interp = vT_interp/abs(omc_interp)

            eval_bp = rhoL_interp**2/2.0d0*(kperp**2 + kperpp**2)
            eval_bt = rhoL_interp**2*kperp*kperpp

            if (kernel_debye_case) then
                A1_interp = 0.0d0
                A2_interp = 0.0d0
            end if

            if (eval_bt > bessel_large_arg_limit) then
                ! limit close to magnetic axis (k_s -> infinity) and
                ! large k_r and k_rp: asymptotics for Bessel I functions
                eval_besselI0 = exp(eval_bt - eval_bp) &
                                /sqrt(2.0d0*pi*eval_bt)
                eval_besselIm1 = exp(-eval_bp + asinh(-1.0d0/eval_bt) &
                                     + eval_bt*sqrt(1.0d0 + 1.0d0/eval_bt**2)) &
                                 /sqrt(2.0d0*pi*eval_bt &
                                       *sqrt(1.0d0 + 1.0d0/eval_bt**2))
            else
                eval_besselI0 = bessel_in(0, eval_bt)*exp(-eval_bp)
                eval_besselIm1 = bessel_in(-1, eval_bt)*exp(-eval_bp)
            end if

            kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg &
                + exp(com_unit*(val_kr - val_krp)*val_rg) &
                *vT_interp**2/(lambda_D_interp**2*omc_interp*kp_interp) &
                *(0.5d0*A1_interp*eval_besselI0 &
                  *(z0_interp*plasma_Z(z0_interp) + 1.0d0) &
                  + A2_interp*(0.5d0*eval_besselI0 &
                               + (z0_interp*plasma_Z(z0_interp) + 1.0d0) &
                               *((1.0d0 + eval_bp + z0_interp**2)*eval_besselI0 &
                                 + eval_bp*eval_besselIm1)))
        end do

        kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg &
                                    *com_unit/(8.0d0*pi**2*sol)

    end function kernel_rho_B_of_kr_krp_rg

end module kernels_m
