module kernel_functions

    use plasma_parameter
    use constants
    use config
    use grid
    use back_quants
    use omp_lib

    implicit none

    logical :: write_out

    double complex :: besselI ! complex bessel function from bessel.f90
    double complex :: plasma_Z ! plasma dispersion function

    !double complex, dimension(:,:), allocatable :: K_rho_phi_limit

    integer :: nlagr = 4
    integer :: nder = 0
    integer :: max_threads
    double complex :: a0, a1, a2
    !double precision :: eval_bp, eval_bt ! b_+ and b_\times
    double complex :: eval_bp, eval_bt ! b_+ and b_\times

    double precision :: b_times_limit = 3d0!3d0
    double precision :: z0_limit = 4.5d0


    contains 

        ! This is without the exp(i k_r(r_g - x_l)) factor
        double complex function kernel_rho_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

            use setup, only: omega
            use constants, only: pi
            implicit none

            double precision, intent(in) :: val_kr, val_krp, val_rg
            double complex :: a0, a1, a2
            double complex :: eval_bp, eval_bt ! b_+ and b_\times
            double complex :: besselI ! complex bessel function from bessel.f90
            double complex :: plasma_Z ! plasma dispersion function

            integer :: sigma ! for loop over species
            integer :: ibeg, iend
            integer :: ir

            double precision, dimension(:,:), allocatable :: coef

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res, nu_res
            double complex :: z0_res
            double complex :: eval_besselI0, eval_besselIm1

            kernel_rho_phi_of_kr_krp_rg = 0.0d0

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            call binsrc(r_prof, 1, iprof_length, val_rg, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, val_rg, r_prof(ibeg:iend), coef)

            ks_res = sum(coef(0,:) * ks(ibeg:iend))
            kp_res = sum(coef(0,:) * kp(ibeg:iend))
            om_E_res = sum(coef(0,:) * om_E(ibeg:iend))

            do sigma = 0, ispecies
                if (sigma == 0) then ! electrons
                    vT_res = sum(coef(0,:) * vTe(ibeg:iend))
                    omc_res = sum(coef(0,:) * omce(ibeg:iend))
                    A1_res = sum(coef(0,:) * A1e(ibeg:iend))
                    A2_res = sum(coef(0,:) * A2e(ibeg:iend))
                    lambda_D_res = sum(coef(0,:) * lambda_De(ibeg:iend))
                    z0_res = sum(coef(0,:) * z0e(ibeg:iend))
                    nu_res = sum(coef(0,:) * nue(ibeg:iend))
                else
                    vT_res = sum(coef(0,:) * vTi(sigma, ibeg:iend))
                    omc_res = sum(coef(0,:) * omci(sigma, ibeg:iend))
                    A1_res = sum(coef(0,:) * A1i(sigma,ibeg:iend))
                    A2_res = sum(coef(0,:) * A2i(sigma,ibeg:iend))
                    lambda_D_res = sum(coef(0,:) * lambda_Di(sigma,ibeg:iend))
                    z0_res = sum(coef(0,:) * z0i(sigma, ibeg:iend))
                    nu_res = sum(coef(0,:) * nui(sigma, ibeg:iend))
                end if

                eval_bp = vT_res**2.0d0 / (2.0d0 * omc_res**2.0d0) * (2.0d0 * ks_res**2.0d0 &
                        + val_kr**2.0d0 + val_krp**2.0d0)
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + val_kr**2.0d0)&
                        * sqrt(ks_res**2.0d0 + val_krp**2.0d0)

                if (real(eval_bt) > b_times_limit) then
                    ! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                    ! use asymptotics for Bessel I functions
                    eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                    eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1/eval_bt**2.0d0)) &
                            / (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                else
                    eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                    eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)
                end if 
                a0 = eval_besselI0 * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                    / (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) &
                    * A2_res * eval_bt * eval_besselIm1 ! *exp(-eval_bp) 
                a1 = - kp_res/omc_res * eval_besselI0 ! * exp(-eval_bp) 
                a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * eval_besselI0 ! * exp(-eval_bp) 

                ! large z limit applies for large k_parallel, i.e. close to the resonant surface
                if (abs(z0_res) > z0_limit) then
                !if (abs(kp_res) < 0.0002d0) then
                    ! with 1/z_0^2 and 1/z_0 terms
                    kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg &
                        + 1.0d0 / lambda_D_res**2.0d0  &
                        * ((a0 * (1.0d0 + (2.0d0 * vT_res**2.0d0 *kp_res**2.0d0)/ (om_E_res - omega - com_unit * nu_res)**2.0d0) &
                        - (vT_res**2.0d0 * kp_res * a1)/ (om_E_res - omega - com_unit * nu_res) &
                        + vT_res**2.0d0 * a2)&
                        * omc_res / (om_E_res - omega - com_unit * nu_res) &
                        - (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (val_krp - val_kr)**2d0) &
                        - eval_besselI0)) 
                else
                    ! ideal region, k_parallel not near zero
                    kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg + 1.0d0 /(lambda_D_res**2d0) &
                            * (omc_res/(sqrt(2.0d0) * kp_res * vT_res) * (plasma_Z(z0_res) &
                            * (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                            * sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) &
                            - (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (val_krp - val_kr)**2d0) - eval_besselI0))
                end if

            end do

            deallocate(coef)

            kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg / (2d0**3 * pi**2)

        end function kernel_rho_phi_of_kr_krp_rg


        ! TODO: implement the following functions
        double complex function kernel_rho_B_of_kr_krp_rg(val_kr, val_krp, val_rg)

            use setup, only: omega
            use constants, only: sol, com_unit, pi
            implicit none
            double precision, intent(in) :: val_kr, val_krp, val_rg
            double complex :: besselI
            double complex :: plasma_Z

            double complex :: eval_besselI0 = 0.0d0
            double complex :: eval_besselIm1 = 0.0d0

            integer :: sigma ! for loop over species
            integer :: ibeg, iend
            integer :: ir

            double precision, dimension(:,:), allocatable :: coef

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res, nu_res
            double complex :: z0_res
            double complex :: eval_bp, eval_bt


            kernel_rho_B_of_kr_krp_rg = 0.0d0

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            call binsrc(r_prof, 1, iprof_length, val_rg, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, val_rg, r_prof(ibeg:iend), coef)

            ks_res = sum(coef(0,:) * ks(ibeg:iend))
            kp_res = sum(coef(0,:) * kp(ibeg:iend))
            om_E_res = sum(coef(0,:) * om_E(ibeg:iend))

            do sigma = 0, ispecies
                if (sigma == 0) then ! electrons
                    vT_res = sum(coef(0,:) * vTe(ibeg:iend))
                    omc_res = sum(coef(0,:) * omce(ibeg:iend))
                    A1_res = sum(coef(0,:) * A1e(ibeg:iend))
                    A2_res = sum(coef(0,:) * A2e(ibeg:iend))
                    lambda_D_res = sum(coef(0,:) * lambda_De(ibeg:iend))
                    z0_res = sum(coef(0,:) * z0e(ibeg:iend))
                    nu_res = sum(coef(0,:) * nue(ibeg:iend))
                else
                    vT_res = sum(coef(0,:) * vTi(sigma, ibeg:iend))
                    omc_res = sum(coef(0,:) * omci(sigma, ibeg:iend))
                    A1_res = sum(coef(0,:) * A1i(sigma,ibeg:iend))
                    A2_res = sum(coef(0,:) * A2i(sigma,ibeg:iend))
                    lambda_D_res = sum(coef(0,:) * lambda_Di(sigma,ibeg:iend))
                    z0_res = sum(coef(0,:) * z0i(sigma, ibeg:iend))
                    nu_res = sum(coef(0,:) * nui(sigma, ibeg:iend))
                end if

                eval_bp = vT_res**2.0d0 / (2.0d0 * omc_res**2.0d0) * (2.0d0 * ks_res**2.0d0 &
                        + val_kr**2.0d0 + val_krp**2.0d0)
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + val_kr**2.0d0)&
                        * sqrt(ks_res**2.0d0 + val_krp**2.0d0)
                
                if (real(eval_bt) > b_times_limit) then
                    ! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                    !write(*,*) 'limit close to magnetic axis, val_rg = ', val_rg 
                    eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                    eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)) &
                            / (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                    !write(*,*) 'eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1 
                else
                    eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                    eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)
                end if

                if (abs(z0_res) > z0_limit) then
                    ! limit close to resonant surface, z -> infinity, k_parallel -> 0
                    kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg - vT_res**4.0d0 * kp_res &
                        / (lambda_D_res**2.0d0 * omc_res)  &
                        / (om_E_res - omega - com_unit * nu_res)**2.0d0 * ((A1_res + A2_res * (1.0d0 + eval_bp)) &
                        * eval_besselI0 + A2_res * eval_bt * eval_besselIm1)
                    !write(*,*) 'spec = ', sigma, ', r = ', r, ', lambda_D = ', lambda_D_res, &!', omc_res = ', omc_res!eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1!', 
                else
                    ! ideal region
                    kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg + vT_res**2.0d0 & 
                        / (lambda_D_res**2.0d0 * omc_res * kp_res) &
                          * ((z0_res * plasma_Z(z0_res) + 1.0d0) &
                          * ((A1_res + A2_res * (1 + eval_bp + z0_res**2.0d0)) * eval_besselI0 &
                          + A2_res * eval_bp * eval_besselIm1) + 0.5d0 * A2_res * eval_besselI0)
                end if

            end do

            deallocate(coef)

            kernel_rho_B_of_kr_krp_rg = kernel_rho_B_of_kr_krp_rg * com_unit / (8d0 * pi**2)! * sol)

        end function kernel_rho_B_of_kr_krp_rg


        double complex function kernel_j_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

            implicit none
            double precision, intent(in) :: val_kr, val_krp, val_rg

            kernel_j_phi_of_kr_krp_rg = 0.0d0

        end function kernel_j_phi_of_kr_krp_rg

        double complex function kernel_j_B_of_kr_krp_rg(val_kr, val_krp, val_rg)

            implicit none
            double precision, intent(in) :: val_kr, val_krp, val_rg

            kernel_j_B_of_kr_krp_rg = 0.0d0

        end function kernel_j_B_of_kr_krp_rg

end module