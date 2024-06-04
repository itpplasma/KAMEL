module kernels

    use plasma_parameter
    use constants
    use config
    use grid
    use back_quants
    use omp_lib

    implicit none

    double complex, dimension(:,:,:), allocatable :: K_rho_phi_of_rg
    double complex, dimension(:,:,:), allocatable :: K_rho_B_of_rg

    double complex, dimension(:,:), allocatable :: K_rho_phi_llp
    double complex, dimension(:,:), allocatable :: K_rho_B_llp

    double complex, dimension(:,:), allocatable :: K_j_phi_llp
    double complex, dimension(:,:), allocatable :: K_j_B_llp

    logical :: write_out

    integer :: nlagr = 4
    integer :: nder = 0
    integer :: max_threads

    double precision :: bessel_large_arg_limit = 3d0
    double precision :: large_z0_limit = 4.5d0


    contains
        subroutine fill_rho_kernels

            use config, only: fstatus
            use loading_bar
            use grid, only: l_space_dim, r_space_dim, varphi_lkr, npoib, rb
            use kr_grid, only: k_space_dim, kr, krp

            implicit none
            integer :: i_kr, i_krp, i_rg
            integer :: count_loading = 0

            if (fstatus == 1) write(*,*) 'Status: Fill rho kernels'
            if (.not. allocated(K_rho_phi_of_rg)) allocate(K_rho_phi_of_rg(k_space_dim, k_space_dim, npoib))
            if (.not. allocated(K_rho_B_of_rg)) allocate(K_rho_B_of_rg(k_space_dim, k_space_dim, npoib))

            K_rho_phi_of_rg = 0.0d0
            K_rho_B_of_rg = 0.0d0

            !$OMP PARALLEL DO collapse(3) default(none) schedule(guided) &
            !$OMP PRIVATE(i_krp, i_kr, i_rg) &
            !$OMP SHARED(K_rho_phi_of_rg, kr, K_rho_B_of_rg, &
            !$OMP krp, rb, k_space_dim, npoib, count_loading)
            do i_krp = 1, k_space_dim
                do i_kr = 1, k_space_dim
                    do i_rg = 1, npoib
                        K_rho_phi_of_rg(i_krp, i_kr, i_rg) = kernel_rho_phi_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                        !K_rho_B_of_rg(i_krp, i_kr, i_rg) = kernel_rho_B_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 

                        if (isnan(real(K_rho_phi_of_rg(i_krp, i_kr, i_rg)))) then
                            write(*,*) "K_rho_phi_of_rg: NaN detected, kr = ", kr(i_kr), ", krp = ", krp(i_krp), ", rb = ", rb(i_rg)
                        end if
                        !write(*,*) 'K_rho_phi_of_rg: ', i_krp, ', ', i_kr, ', ', i_rg, ', ', K_rho_phi_of_rg(i_krp, i_kr, i_rg)
                    end do
                end do
            end do
            !$OMP END PARALLEL DO

            if (fstatus == 1) write(*,*) 'Status: Finished filling rho kernels'

        end subroutine fill_rho_kernels

        
        ! This is without the exp(i k_r(r_g - x_l)) factor
        double complex function kernel_rho_phi_of_kr_krp_rg(val_kr, val_krp, val_rg)

            use setup, only: omega
            use constants, only: pi
            implicit none

            double precision, dimension(:,:), allocatable :: coef
            integer :: ibeg, iend            
            double precision, intent(in) :: val_kr, val_krp, val_rg
            integer :: ir

            ! sub functions appearing in the kernels
            double complex :: a0, a1, a2
            double complex :: eval_bp, eval_bt ! b_+ and b_\times
            double complex :: z0_interp
            double complex :: eval_besselI0, eval_besselIm1

            double complex :: besselI ! complex bessel function from bessel.f90
            double complex :: plasma_Z ! plasma dispersion function

            integer :: sigma ! for loop over species

            ! interpolated values of the parameters
            double precision :: vT_interp, omc_interp, ks_interp, om_E_interp, kp_interp, &
                                        A1_interp, A2_interp, lambda_D_interp, nu_interp

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
                    eval_besselI0 = cmplx(0.0d0, 0.0d0)
                    eval_besselIm1 = cmplx(0.0d0, 0.0d0)
                    A1_interp = 0.0d0
                    A2_interp = 0.0d0
                    z0_interp = 0.0d0
                    a0 = 0.0d0
                    a1 = 0.0d0
                    a2 = 0.0d0

                    kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg + 1.0d0 / lambda_D_interp**2.0d0 &
                                            * exp(-vT_interp**2.0d0 / (2.0d0 * omc_interp ** 2.0d0) &
                                            * (val_krp - val_kr)**2.0d0)
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

                    ! large z limit applies for large k_parallel, i.e. close to the resonant surface
                    if (abs(z0_interp) > large_z0_limit) then
                        ! with 1/z_0^2 and 1/z_0 terms
                        kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg &
                            + 1.0d0 / lambda_D_interp**2.0d0  &
                            * ((a0 * (1.0d0 + (2.0d0 * vT_interp**2.0d0 *kp_interp**2.0d0) &
                            / (om_E_interp - omega - com_unit * nu_interp)**2.0d0) &
                            - (vT_interp**2.0d0 * kp_interp * a1)/ (om_E_interp - omega - com_unit * nu_interp) &
                            + vT_interp**2.0d0 * a2)&
                            * omc_interp / (om_E_interp - omega - com_unit * nu_interp) &
                            - (exp(-vT_interp**2d0/(2d0 * omc_interp**2d0) * (val_krp - val_kr)**2d0) &
                            - eval_besselI0)) 
                    else
                        ! ideal region, k_parallel not near zero
                        kernel_rho_phi_of_kr_krp_rg = kernel_rho_phi_of_kr_krp_rg + 1.0d0 /(lambda_D_interp**2d0) &
                                * (omc_interp/(sqrt(2.0d0) * kp_interp * vT_interp) * (plasma_Z(z0_interp) &
                                * (a0 + sqrt(2.0d0) * vT_interp * z0_interp * a1 + vT_interp**2d0 * a2 * 2d0 * z0_interp**2d0) &
                                + vT_interp * sqrt(2d0) * (a1 + vT_interp * sqrt(2d0) * z0_interp * a2)) &
                                - (exp(-vT_interp**2d0/(2d0 * omc_interp**2d0) * (val_krp - val_kr)**2d0) - eval_besselI0))
                    end if

                end if
            end do

            deallocate(coef)

            kernel_rho_phi_of_kr_krp_rg = -kernel_rho_phi_of_kr_krp_rg / (2d0**3 * pi**2)

            if (isnan(real(kernel_rho_phi_of_kr_krp_rg))) then
                write(*,*) "kernel_rho_phi_of_kr_krp_rg: NaN detected"
            end if

        end function kernel_rho_phi_of_kr_krp_rg


        ! TODO: implement the following functions
        double complex function kernel_rho_B_of_kr_krp_rg(val_kr, val_krp, val_rg)

            use setup, only: omega
            use constants, only: sol, com_unit, pi
            implicit none
            double precision, intent(in) :: val_kr, val_krp, val_rg

            double precision, dimension(:,:), allocatable :: coef
            integer :: ibeg, iend            
            integer :: ir
            
            ! sub functions appearing in the kernels
            double complex :: a0, a1, a2
            double complex :: eval_bp, eval_bt ! b_+ and b_\times
            double complex :: z0_interp
            double complex :: eval_besselI0, eval_besselIm1

            double complex :: besselI ! complex bessel function from bessel.f90
            double complex :: plasma_Z ! plasma dispersion function

            integer :: sigma ! for loop over species

            ! interpolated values of the parameters
            double precision :: vT_interp, omc_interp, ks_interp, om_E_interp, kp_interp, &
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
                    eval_besselI0 = cmplx(0.0d0, 0.0d0)
                    eval_besselIm1 = cmplx(0.0d0, 0.0d0)
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