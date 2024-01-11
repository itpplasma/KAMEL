module kernel_functions

    use plas_parameter
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


    contains 

        ! This is without the exp(i k_r(r_g - x_l)) factor
        double complex function kernel_rho_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

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
                                A1_res, A2_res, lambda_D_res
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
                else
                    vT_res = sum(coef(0,:) * vTi(sigma, ibeg:iend))
                    omc_res = sum(coef(0,:) * omci(sigma, ibeg:iend))
                    A1_res = sum(coef(0,:) * A1i(sigma,ibeg:iend))
                    A2_res = sum(coef(0,:) * A2i(sigma,ibeg:iend))
                    lambda_D_res = sum(coef(0,:) * lambda_Di(sigma,ibeg:iend))
                    z0_res = sum(coef(0,:) * z0i(sigma, ibeg:iend))
                end if

                eval_bp = vT_res**2.0d0 / (2.0d0 * omc_res**2.0d0) * (2.0d0 * ks_res**2.0d0 &
                        + val_kr**2.0d0 + val_krp**2.0d0)
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + val_kr**2.0d0)&
                        * sqrt(ks_res**2.0d0 + val_krp**2.0d0)
                    
                !write(*,*) 'species = ', sigma
                !write(*,*) 'vT_res = ', vT_res, '; omega_c = ', omc_res
                !write(*,*) 'eval_bt = ', eval_bt, '; val_kr = ', val_kr, '; val_krp = ', val_krp
                !write(*,*) 'besselI0 = ', besselI(0,eval_bt,0), '; besselIm1 = ', besselI(-1, eval_bt, 0)

                if (real(eval_bt) > 10d0) then
                    ! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                    ! use asymptotics for Bessel I functions
                    !write(*,*) 'special case'
                    eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                    eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1/eval_bt**2.0d0)) &
                            / (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                    
                    !write(*,*) 'eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1 
                else
                    eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                    eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)
                end if 
                a0 = eval_besselI0 * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                    / (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) &
                    * A2_res * eval_bt * eval_besselIm1 ! *exp(-eval_bp) 
                a1 = - kp_res/omc_res * eval_besselI0 ! * exp(-eval_bp) 
                a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * eval_besselI0 ! * exp(-eval_bp) 

                !a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                !    / (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) * A2_res &
                !    * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                !a1 = - kp_res/omc_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                !a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                
                !write(*,*) 'after as'
                
                kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg  &
                                + omc_res /(lambda_D_res**2d0 * vT_res) &
                                * (sqrt(pi)/kp_res * plasma_Z(z0_res) &
                                * (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                                * sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) - sqrt(2d0 * pi) * vT_res / omc_res &
                                * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (val_krp - val_kr)**2d0) - eval_besselI0)
                !write(*,*) 'after kernel'
            end do

            !write(*,*) 'after species loop'

            deallocate(coef)

        end function kernel_rho_phi_of_kr_krp_rg



        ! This is without the exp(i k_r(r_g - x_l)) factor
        !double complex function kernel_rho_phi_of_kr_krp_rg_same_grid(val_kr, val_krp)

            !implicit none

            !double precision, intent(in) :: val_kr, val_krp
            !double complex :: a0, a1, a2
            !double complex :: eval_bp, eval_bt ! b_+ and b_\times
            !double complex :: besselI ! complex bessel function from bessel.f90
            !double complex :: plasma_Z ! plasma dispersion function

            !integer :: sigma ! for loop over species

            !double complex :: z0_res
            !double complex :: eval_besselI0, eval_besselIm1

            !double precision, dimension(:) :: vT_prof, omc_prof, A1_prof, A2_prof, lambda_D_prof
            !double complex, dimension(:) :: z0_prof, eval_bp, eval_bt

            !allocate(vT_prof(iprof_length), omc_prof(iprof_length), A1_prof(iprof_length), A2_prof(iprof_length), &
                !lambda_D_prof(iprof_length), z0_prof(iprof_length))

            !kernel_rho_phi_of_kr_krp_rg = 0.0d0

            !do sigma = 0, ispecies
                !if (sigma == 0) then ! electrons
                    !vT_prof = vTe
                    !omc_prof = omce
                    !A1_prof = A1e(ibeg:iend)
                    !A2_prof = A2e(ibeg:iend)
                    !lambda_D_prof = lambda_De(ibeg:iend)
                    !z0_prof = z0e(ibeg:iend)
                !else
                    !vT_prof = vTi(sigma, ibeg:iend)
                    !omc_prof = omci(sigma, ibeg:iend)
                    !A1_prof = A1i(sigma,ibeg:iend)
                    !A2_prof = A2i(sigma,ibeg:iend)
                    !lambda_D_prof = lambda_Di(sigma,ibeg:iend)
                    !z0_prof = z0i(sigma, ibeg:iend)
                !end if

                !eval_bp = vT_prof**2.0d0 / (2.0d0 * omc_prof**2.0d0) * (2.0d0 * ks**2.0d0 &
                        !+ val_kr**2.0d0 + val_krp**2.0d0)
                !eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + val_kr**2.0d0)&
                        !* sqrt(ks_res**2.0d0 + val_krp**2.0d0)
                    
                !!write(*,*) 'species = ', sigma
                !!write(*,*) 'vT_res = ', vT_res, '; omega_c = ', omc_res
                !!write(*,*) 'eval_bt = ', eval_bt, '; val_kr = ', val_kr, '; val_krp = ', val_krp
                !!write(*,*) 'besselI0 = ', besselI(0,eval_bt,0), '; besselIm1 = ', besselI(-1, eval_bt, 0)

                !if (real(eval_bt) > 10d0) then
                    !! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                    !! use asymptotics for Bessel I functions
                    !!write(*,*) 'special case'
                    !eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                    !eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1/eval_bt**2.0d0)) &
                            !/ (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                    
                    !!write(*,*) 'eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1 
                !else
                    !eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                    !eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)
                !end if 
                !a0 = eval_besselI0 * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                    !/ (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) &
                    !* A2_res * eval_bt * eval_besselIm1 ! *exp(-eval_bp) 
                !a1 = - kp_res/omc_res * eval_besselI0 ! * exp(-eval_bp) 
                !a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * eval_besselI0 ! * exp(-eval_bp) 

                !!a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                !!    / (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) * A2_res &
                !!    * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                !!a1 = - kp_res/omc_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                !!a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                
                !!write(*,*) 'after as'
                
                !kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg  &
                                !+ omc_res /(lambda_D_res**2d0 * vT_res) &
                                !* (sqrt(pi)/kp_res * plasma_Z(z0_res) &
                                !* (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                                !* sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) - sqrt(2d0 * pi) * vT_res / omc_res &
                                !* (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (val_krp - val_kr)**2d0) - eval_besselI0)
                !!write(*,*) 'after kernel'
            !end do

            !!write(*,*) 'after species loop'

            !deallocate(coef)

        !end function kernel_rho_phi_of_kr_krp_rg

end module