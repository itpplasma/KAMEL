module kernels

    use plasma_parameter
    use constants
    use config
    use grid
    use back_quants
    use omp_lib
    use KIM_kinds, only: dp

    implicit none

    complex(dp), dimension(:,:,:), allocatable :: K_rho_phi_of_rg
    complex(dp), dimension(:,:,:), allocatable :: K_rho_B_of_rg

    complex(dp), dimension(:,:), allocatable :: K_rho_phi_llp
    complex(dp), dimension(:,:), allocatable :: K_rho_B_llp

    complex(dp), dimension(:,:), allocatable :: K_j_phi_llp
    complex(dp), dimension(:,:), allocatable :: K_j_B_llp

    logical :: write_out

    integer :: nlagr = 4
    integer :: max_threads

    real(dp) :: bessel_large_arg_limit = 3d0
    real(dp) :: large_z0_limit = 4.5d0


    contains

        subroutine fill_rho_kernels

            use config, only: fstatus
            use loading_bar
            use grid, only: varphi_lkr, rg_grid, kr_grid, krp_grid

            implicit none
            integer :: i_kr, i_krp, i_rg
            integer :: count_loading = 0

            if (fstatus == 1) write(*,*) 'Status: Fill rho kernels'
            if (.not. allocated(K_rho_phi_of_rg)) allocate(K_rho_phi_of_rg(krp_grid%npts_b, kr_grid%npts_b, rg_grid%npts_b))
            if (.not. allocated(K_rho_B_of_rg)) allocate(K_rho_B_of_rg(kr_grid%npts_b, krp_grid%npts_b, rg_grid%npts_b))

            K_rho_phi_of_rg = 0.0d0
            K_rho_B_of_rg = 0.0d0

            !$OMP PARALLEL DO collapse(3) default(none) schedule(guided) &
            !$OMP PRIVATE(i_krp, i_kr, i_rg) &
            !$OMP SHARED(K_rho_phi_of_rg, K_rho_B_of_rg, &
            !$OMP kr_grid, krp_grid, rg_grid, count_loading)
            do i_krp = 1, krp_grid%npts_b
                do i_kr = 1, kr_grid%npts_b
                    do i_rg = 1, rg_grid%npts_b
                        K_rho_phi_of_rg(i_krp, i_kr, i_rg) = kernel_rho_phi_of_kr_krp_rg(krp_grid%xb(i_kr), kr_grid%xb(i_krp)&
                                                            , rg_grid%xb(i_rg))
                        !K_rho_B_of_rg(i_krp, i_kr, i_rg) = kernel_rho_B_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg))

                        if (isnan(real(K_rho_phi_of_rg(i_krp, i_kr, i_rg)))) then
                            write(*,*) "K_rho_phi_of_rg: NaN detected, kr = ", kr_grid%xb(i_kr), ", krp = ", krp_grid%xb(i_krp)&
                                        , ", rg = ", rg_grid%xb(i_rg)
                        end if
                    end do
                end do
            end do
            !$OMP END PARALLEL DO

            if (fstatus == 1) write(*,*) 'Status: Finished filling rho kernels'

        end subroutine fill_rho_kernels


        ! This is without the exp(i k_r(r_g - x_l)) factor
        complex(dp) function kernel_rho_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

            use setup, only: omega
            use constants, only: pi
            use KIM_kinds, only: dp
            use gsl_mod, only: gsl_sf_bessel_In

            implicit none

            real(dp), dimension(:,:), allocatable :: coef
            integer :: ibeg, iend
            real(dp), intent(in) :: val_kr, val_krp, val_rg
            integer :: ir

            ! sub functions appearing in the kernels
            complex(dp) :: eval_bp, eval_bt ! b_+ and b_\times
            complex(dp) :: z0_interp
            complex(dp) :: eval_besselI0, eval_besselIm1

            complex(dp) :: besselI ! complex bessel function from bessel.f90
            complex(dp) :: plasma_Z ! plasma dispersion function

            integer :: sigma ! for loop over species

            ! interpolated values of the parameters
            real(dp) :: vT_interp, omc_interp, ks_interp, om_E_interp, kp_interp, &
                        A1_interp, A2_interp, lambda_D_interp, nu_interp, rhoL_interp,&
                        kperp, kperpp

            kernel_rho_phi_of_kr_krp_rg = 0.0d0

            if(.not. allocated(coef)) then
                allocate(coef(0:nder, nlagr))
            end if

            call binsrc(r_prof, 1, iprof_length, val_rg, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, val_rg, r_prof(ibeg:iend), coef)

            ks_interp = sum(coef(0,:) * ks(ibeg:iend))
            kp_interp = sum(coef(0,:) * kp(ibeg:iend))
            om_E_interp = sum(coef(0,:) * om_E(ibeg:iend))
            kperp = sqrt(ks_interp**2 + val_kr**2)
            kperpp = sqrt(ks_interp**2 + val_krp**2)

            print *, ""
            print *, "ks_interp = ", ks_interp
            print *, "kp_interp = ", kp_interp
            print *, "kperp = ", kperp
            print *, "kperpp = ", kperpp
            print *, ""

            do sigma = 0, number_of_ion_species

                if (sigma == 0) then ! electrons
                    vT_interp = sum(coef(0,:) * vTe(ibeg:iend))
                    omc_interp = sum(coef(0,:) * omce(ibeg:iend))
                    A1_interp = sum(coef(0,:) * A1e(ibeg:iend))
                    A2_interp = sum(coef(0,:) * A2e(ibeg:iend))
                    lambda_D_interp = sum(coef(0,:) * lambda_De(ibeg:iend))
                    z0_interp = sum(coef(0,:) * z0e(ibeg:iend))
                    nu_interp = sum(coef(0,:) * nue(ibeg:iend))
                else
                    vT_interp = sum(coef(0,:) * vTi(sigma, ibeg:iend))
                    omc_interp = sum(coef(0,:) * omci(sigma, ibeg:iend))
                    A1_interp = sum(coef(0,:) * A1i(sigma,ibeg:iend))
                    A2_interp = sum(coef(0,:) * A2i(sigma,ibeg:iend))
                    lambda_D_interp = sum(coef(0,:) * lambda_Di(sigma,ibeg:iend))
                    z0_interp = sum(coef(0,:) * z0i(sigma, ibeg:iend))
                    nu_interp = sum(coef(0,:) * nui(sigma, ibeg:iend))
                end if

                rhoL_interp = vT_interp/abs(omc_interp)
                print *, ""

                eval_bp = rhoL_interp**2.0d0 / 2.0d0 * (kperp**2.0d0 + kperpp**2.0d0)
                eval_bt = rhoL_interp**2.0d0 * kperp * kperpp

                print *, "z0_interp = ", z0_interp
                print *, "vT_interp = ", vT_interp
                print *, "omc_interp = ", omc_interp
                print *, "A1_interp = ", A1_interp
                print *, "A2_interp = ", A2_interp
                print *, "rhoL_interp = ", rhoL_interp
                print *, "plasma_Z = ", plasma_Z(z0_interp)
                print *, "eval_bp = ", eval_bp
                print *, "eval_bt = ", eval_bt
                print *, "Lambda_D = ", lambda_D_interp


                if (kernel_debye_case .eqv. .true.)then
                    eval_besselI0 = 0.0d0
                    eval_besselIm1 = 0.0d0
                    A1_interp = 0.0d0
                    A2_interp = 0.0d0
                    z0_interp = 0.0d0

                    kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg + 1.0d0 / lambda_D_interp**2.0d0 &
                                            * exp(-vT_interp**2.0d0 / (2.0d0 * omc_interp ** 2.0d0) &
                                            * (val_kr - val_krp)**2.0d0)
                else
                    if (real(eval_bt) > bessel_large_arg_limit) then
                        ! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                        ! use asymptotics for Bessel I functions
                        eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                        eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1/eval_bt**2.0d0)) &
                                / (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                    else
                        !eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                        !eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)

                        eval_besselI0 = gsl_sf_bessel_In(0, real(eval_bt, dp)) * exp(-eval_bp)
                        eval_besselIm1 = gsl_sf_bessel_In(-1, real(eval_bt, dp)) * exp(-eval_bp)

                    end if

                    print *, ""
                    print *, "eval_besselI0 = ", eval_besselI0
                    print *, "Bessel0 = ", gsl_sf_bessel_In(0, real(eval_bt,dp))
                    print *, "eval_besselIm1 = ", eval_besselIm1
                    print *, "Bessel-1 = ", gsl_sf_bessel_In(-1, real(eval_bt,dp))

                    kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg + &
                        1.0d0/(lambda_D_interp**2.0d0) * exp(com_unit * (val_kr - val_krp) * val_rg) &
                        * (-exp(-rhoL_interp**2.0d0/ 2.0d0 * (val_kr - val_krp)**2.0d0) &
                        + ks_interp * rhoL_interp / (kp_interp * sqrt(2.0d0)) &
                            * (A1_interp * eval_besselI0 * plasma_Z(z0_interp)  &
                            + A2_interp * (plasma_Z(z0_interp) * eval_besselI0 * (1 + eval_bp + z0_interp**2.0d0) &
                            + eval_besselIm1 * eval_bt + z0_interp * eval_besselI0))&
                        )
                    print *, "kernel_rho_phi_of_kr_krp_rg = ", kernel_rho_phi_of_kr_krp_rg
                end if
            end do

            deallocate(coef)

            kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg / (2d0**3 * pi**2)

            if (isnan(real(kernel_rho_phi_of_kr_krp_rg))) then
                write(*,*) "kernel_rho_phi_of_kr_krp_rg: NaN detected"
            end if

        end function kernel_rho_phi_of_kr_krp_rg


        ! TODO: implement the following functions
        complex(dp) function kernel_rho_B_of_kr_krp_rg(val_kr, val_krp, val_rg)

            use setup, only: omega
            use constants, only: sol, com_unit, pi
            use KIM_kinds, only: dp

            implicit none

            real(dp), intent(in) :: val_kr, val_krp, val_rg

            real(dp), dimension(:,:), allocatable :: coef
            integer :: ibeg, iend
            integer :: ir

            ! sub functions appearing in the kernels
            complex(dp) :: a0, a1, a2
            complex(dp) :: eval_bp, eval_bt ! b_+ and b_\times
            complex(dp) :: z0_interp
            complex(dp) :: eval_besselI0, eval_besselIm1

            complex(dp) :: besselI ! complex bessel function from bessel.f90
            complex(dp) :: plasma_Z ! plasma dispersion function

            integer :: sigma ! for loop over species

            ! interpolated values of the parameters
            real(dp) :: vT_interp, omc_interp, ks_interp, om_E_interp, kp_interp, &
                                        A1_interp, A2_interp, lambda_D_interp, nu_interp


            kernel_rho_B_of_kr_krp_rg = 0.0d0

            if(.not. allocated(coef)) then
                allocate(coef(0:nder, nlagr))
            end if

            call binsrc(r_prof, 1, iprof_length, val_rg, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, val_rg, r_prof(ibeg:iend), coef)

            ks_interp = sum(coef(0,:) * ks(ibeg:iend))
            kp_interp = sum(coef(0,:) * kp(ibeg:iend))
            om_E_interp = sum(coef(0,:) * om_E(ibeg:iend))

            do sigma = 0, number_of_ion_species

                if (sigma == 0) then ! electrons
                    vT_interp = sum(coef(0,:) * vTe(ibeg:iend))
                    omc_interp = sum(coef(0,:) * omce(ibeg:iend))
                    A1_interp = sum(coef(0,:) * A1e(ibeg:iend))
                    A2_interp = sum(coef(0,:) * A2e(ibeg:iend))
                    lambda_D_interp = sum(coef(0,:) * lambda_De(ibeg:iend))
                    z0_interp = sum(coef(0,:) * z0e(ibeg:iend))
                    nu_interp = sum(coef(0,:) * nue(ibeg:iend))
                else
                    vT_interp = sum(coef(0,:) * vTi(sigma, ibeg:iend))
                    omc_interp = sum(coef(0,:) * omci(sigma, ibeg:iend))
                    A1_interp = sum(coef(0,:) * A1i(sigma,ibeg:iend))
                    A2_interp = sum(coef(0,:) * A2i(sigma,ibeg:iend))
                    lambda_D_interp = sum(coef(0,:) * lambda_Di(sigma,ibeg:iend))
                    z0_interp = sum(coef(0,:) * z0i(sigma, ibeg:iend))
                    nu_interp = sum(coef(0,:) * nui(sigma, ibeg:iend))
                end if

                eval_bp = vT_interp**2.0d0 / (2.0d0 * omc_interp**2.0d0) * (2.0d0 * ks_interp**2.0d0 &
                        + val_kr**2.0d0 + val_krp**2.0d0)
                eval_bt = vT_interp**2.0d0 /(omc_interp**2.0d0) * sqrt(ks_interp**2.0d0 + val_kr**2.0d0)&
                        * sqrt(ks_interp**2.0d0 + val_krp**2.0d0)

                if (kernel_debye_case .eqv. .true.)then
                    eval_besselI0 = 0.0d0
                    eval_besselIm1 = 0.0d0
                    A1_interp = 0.0d0
                    A2_interp = 0.0d0
                    z0_interp = 0.0d0
                    a0 = 0.0d0
                    a1 = 0.0d0
                    a2 = 0.0d0
                else
                    if (real(eval_bt) > bessel_large_arg_limit) then
                        ! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                        ! use asymptotics for Bessel I functions
                        eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                        eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1/eval_bt**2.0d0)) &
                                / (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                    else
                        eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                        eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)
                    end if

                    a0 = eval_besselI0 * (- om_E_interp / omc_interp + ks_interp * vT_interp**2d0 &
                    / (omc_interp**2d0) * (A1_interp + (1.0d0 + eval_bp) * A2_interp)) &
                    + ks_interp * vT_interp**2d0 / (omc_interp**2d0) &
                    * A2_interp * eval_bt * eval_besselIm1 ! *exp(-eval_bp)

                    a1 = - kp_interp/omc_interp * eval_besselI0 ! * exp(-eval_bp)
                    a2 = ks_interp / (2d0 * omc_interp**2d0) * A2_interp * eval_besselI0 ! * exp(-eval_bp)

                end if

                if (abs(z0_interp) > large_z0_limit) then
                    ! limit close to resonant surface, z -> infinity, k_parallel -> 0
                    kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg - vT_interp**4.0d0 * kp_interp &
                        / (lambda_D_interp**2.0d0 * omc_interp)  &
                        / (om_E_interp - omega - com_unit * nu_interp)**2.0d0 * ((A1_interp + A2_interp * (1.0d0 + eval_bp)) &
                        * eval_besselI0 + A2_interp * eval_bt * eval_besselIm1)
                    !write(*,*) 'spec = ', sigma, ', r = ', r, ', lambda_D = ', lambda_D_interp, &!', omc_interp = ', omc_interp!eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1!',
                else
                    ! ideal region
                    kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg + vT_interp**2.0d0 &
                        / (lambda_D_interp**2.0d0 * omc_interp * kp_interp) &
                          * ((z0_interp * plasma_Z(z0_interp) + 1.0d0) &
                          * ((A1_interp + A2_interp * (1 + eval_bp + z0_interp**2.0d0)) * eval_besselI0 &
                          + A2_interp * eval_bp * eval_besselIm1) + 0.5d0 * A2_interp * eval_besselI0)
                end if

            end do

            deallocate(coef)

            kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg * com_unit / (8d0 * pi**2 * sol)

        end function kernel_rho_B_of_kr_krp_rg


        complex(dp) function kernel_j_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

            implicit none
            real(dp), intent(in) :: val_kr, val_krp, val_rg

            kernel_j_phi_of_kr_krp_rg = 0.0d0

        end function kernel_j_phi_of_kr_krp_rg

        complex(dp) function kernel_j_B_of_kr_krp_rg(val_kr, val_krp, val_rg)

            implicit none
            real(dp), intent(in) :: val_kr, val_krp, val_rg

            kernel_j_B_of_kr_krp_rg = 0.0d0

        end function kernel_j_B_of_kr_krp_rg


end module
