module integrands 

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants
    use setup, only: omega

    implicit none

    integer :: nlagr = 4
    integer :: nder = 0

    contains

        subroutine integrand_K_rho_phi_rg(r, y, dydr, int_rg)

            use grid, only: npoib
            implicit none
            double precision, intent(in) :: r
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr
            double complex, dimension(npoib), intent(in) :: int_rg
            double complex :: int_intp

            integer :: ibeg, iend , ir
            double precision, dimension(:,:), allocatable :: coef
            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            dydr = 0.0d0
            int_intp = 0.0d0

            call binsrc(r, 1, npoib, r_prof, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. npoib) then
                iend = npoib
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, r, r_prof(ibeg:iend), coef)

            int_intp = sum(cmplx(coef(0,:)) * int_rg(ibeg:iend))

            dydr = dydr + int_intp
 

        end subroutine

        ! integrand K_{lp}(kr; r_g) used for integration over kr
        ! second integral after krp integration
        subroutine integrand_K_rho_phi_kr(c_kr, y, dydr, int_kr_r, j, l)

            implicit none

            integer, intent(in) :: j
            integer, intent(in) :: l
            double precision, intent(in) :: c_kr
            double complex, dimension(k_space_dim, npoib), intent(in) :: int_kr_r
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr
            double complex :: int_intp
            double complex :: varphi

            integer :: ibeg, iend , ikr
            double precision, dimension(:,:), allocatable :: coef
            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            dydr = 0.0d0
            int_intp = 0.0d0
            varphi = 0.0d0

            call binsrc(c_kr, 1, k_space_dim, kr, ikr)
            ibeg = max(1, ikr - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. k_space_dim) then
                iend = k_space_dim
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, c_kr, kr(ibeg:iend), coef)

            int_intp = sum(cmplx(coef(0,:)) * int_kr_r(ibeg:iend, j))
            varphi = sum(cmplx(coef(0,:)) * varphi_lkr(l, ibeg:iend))

            dydr = dydr + int_intp * varphi
            !write(*,*) 'dydr = ', dydr
            !write(*,*) 'int_intp = ', int_intp
            !write(*,*) 'varphi = ', varphi

        end subroutine

        ! integrand K(kr', kr; r_g) used for integration over kr'
        subroutine integrand_K_rho_phi_krp(c_krp, y, dydr, j, k1, lp)

            implicit none

            double precision, intent(in) :: c_krp
            integer, intent(in) :: j     ! index for r profile
            integer, intent(in) :: k1    ! index for kr profile
            integer, intent(in) :: lp ! indices for spline basis functions
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr
            double complex :: a0, a1, a2
            double complex :: eval_bp, eval_bt ! b_+ and b_\times
            double complex :: besselI
            double complex :: plasma_Z
            double complex :: eval_besselI0, eval_besselIm1

            integer :: sigma ! for loop over species
            ! for interpolation
            integer :: ibeg, iend , ikrp
            double precision, dimension(:,:), allocatable :: coef
            double complex :: varphi

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res, nu_res
            double complex :: z0_res

            double precision :: r

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            r = r_prof(j)
            dydr = 0.0d0

            ks_res = ks(j)
            kp_res = kp(j)

            do sigma = 0, ispecies
            !sigma = 0
                if (sigma == 0) then ! electrons
                    vT_res = vTe(j)
                    omc_res = omce(j)
                    A1_res = A1e(j)
                    A2_res = A2e(j)
                    lambda_D_res = lambda_De(j)
                    z0_res = z0e(j)
                    nu_res = nue(j)
                else
                    vT_res = vTi(sigma, j)
                    omc_res = omci(sigma, j)
                    A1_res = A1i(sigma, j)
                    A2_res = A2i(sigma, j)
                    lambda_D_res = lambda_Di(sigma, j)
                    z0_res = z0i(sigma, j)
                    nu_res = nui(sigma, j)
                end if

                eval_bp = vT_res**2.0d0 / (2.0d0 * omc_res**2.0d0) * (2.0d0 * ks_res**2.0d0 &
                        + kr(k1)**2.0d0 + c_krp**2.0d0)
                    
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + kr(k1)**2.0d0)&
                        * sqrt(ks_res**2.0d0 + c_krp**2.0d0)

                if (real(eval_bt) > 10d0) then
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

                !write(*,*) 'a0 = ', a0
                !write(*,*) 'a1 = ', a1
                !write(*,*) 'a2 = ', a2
                !write(*,*) 'eval_besselI0 = ', eval_besselI0
                !write(*,*) 'ks_res = ', ks_res

                ! large z limit
                if (abs(z0_res) > 4.5d0) then
                !if (abs(kp_res) < 0.0002d0) then
                    ! with 1/z_0^2 and 1/z_0 terms
                    !write(*,*) 'large z limit'
                    dydr = dydr + omc_res / lambda_D_res**2.0d0 * exp(com_unit * (kr(k1) - c_krp) * r) &
                        * ((a0 * (1.0d0 + (2.0d0 * vT_res**2.0d0 *kp_res**2.0d0)/ (om_E_res - omega - com_unit * nu_res)**2.0d0) &
                        - (vT_res**2.0d0 * kp_res * a1)/ (om_E_res - omega - com_unit * nu_res) &
                        + vT_res**2.0d0 * a2)&
                        / (om_E_res - omega - com_unit * nu_res) &
                        - 1.0d0 / omc_res * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - kr(k1))**2d0) &
                        - eval_besselI0)) !exp(-eval_bp) * eval_besselI0))
                else
                    ! ideal region, k_parallel not near zero
                    !write(*,*) 'NOT large z limit'
                    dydr = dydr + omc_res /(lambda_D_res**2d0 * vT_res) &
                            * exp(com_unit * (kr(k1) - c_krp) * r) * (sqrt(pi)/kp_res * plasma_Z(z0_res) &
                            * (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                            * sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) - sqrt(2d0 * pi) * vT_res / omc_res &
                            * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - kr(k1))**2d0) - eval_besselI0)!exp(-eval_bp) &
                            !* eval_besselI0)
                end if

                dydr = dydr * lambda_D_res**2d0

            end do

            call binsrc(c_krp, 1, k_space_dim, krp, ikrp)
            ibeg = max(1, ikrp - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. k_space_dim) then
                iend = k_space_dim
                ibeg = iend - nlagr + 1
            end if
            call plag_coeff(nlagr, nder, c_krp, krp(ibeg:iend), coef)

            varphi = sum(cmplx(coef(0,:)) * varphi_lkr(lp, ibeg:iend))
            dydr = dydr * conjg(varphi)
            !write(*,*) 'dydr = ', dydr
            !write(*,*) 'varphi = ', varphi_lkr(lp, ibeg:iend)
            !write(*,*) 'coef = ', cmplx(coef(0,:))
        end subroutine

        ! Kernel without varphis s
        subroutine integrand_K_rho_phi_krp_kr_rg(c_krp, y, dydr, j, k1, lp)

            implicit none

            double precision, intent(in) :: c_krp
            integer, intent(in) :: j     ! index for r profile
            integer, intent(in) :: k1    ! index for kr profile
            integer, intent(in) :: lp ! indices for spline basis functions
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr
            double complex :: a0, a1, a2
            double complex :: eval_bp, eval_bt ! b_+ and b_\times
            double complex :: besselI
            double complex :: plasma_Z
            double complex :: eval_besselI0, eval_besselIm1

            integer :: sigma ! for loop over species
            ! for interpolation
            integer :: ibeg, iend , ikrp
            double precision, dimension(:,:), allocatable :: coef
            double complex :: varphi

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res, nu_res
            double complex :: z0_res

            double precision :: r

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            r = r_prof(j)
            dydr = 0.0d0

            ks_res = ks(j)
            kp_res = kp(j)

            do sigma = 0, ispecies
            !sigma = 0
                if (sigma == 0) then ! electrons
                    vT_res = vTe(j)
                    omc_res = omce(j)
                    A1_res = A1e(j)
                    A2_res = A2e(j)
                    lambda_D_res = lambda_De(j)
                    z0_res = z0e(j)
                    nu_res = nue(j)
                else
                    vT_res = vTi(sigma, j)
                    omc_res = omci(sigma, j)
                    A1_res = A1i(sigma, j)
                    A2_res = A2i(sigma, j)
                    lambda_D_res = lambda_Di(sigma, j)
                    z0_res = z0i(sigma, j)
                    nu_res = nui(sigma, j)
                end if

                eval_bp = vT_res**2.0d0 / (2.0d0 * omc_res**2.0d0) * (2.0d0 * ks_res**2.0d0 &
                        + kr(k1)**2.0d0 + c_krp**2.0d0)
                    
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + kr(k1)**2.0d0)&
                        * sqrt(ks_res**2.0d0 + c_krp**2.0d0)

                if (real(eval_bt) > 10d0) then
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

                !write(*,*) 'a0 = ', a0
                !write(*,*) 'a1 = ', a1
                !write(*,*) 'a2 = ', a2
                !write(*,*) 'eval_besselI0 = ', eval_besselI0
                !write(*,*) 'ks_res = ', ks_res

                ! large z limit
                if (abs(z0_res) > 4.5d0) then
                !if (abs(kp_res) < 0.0002d0) then
                    ! with 1/z_0^2 and 1/z_0 terms
                    !write(*,*) 'large z limit'
                    dydr = dydr + omc_res / lambda_D_res**2.0d0 * exp(com_unit * (kr(k1) - c_krp) * r) &
                        * ((a0 * (1.0d0 + (2.0d0 * vT_res**2.0d0 *kp_res**2.0d0)/ (om_E_res - omega - com_unit * nu_res)**2.0d0) &
                        - (vT_res**2.0d0 * kp_res * a1)/ (om_E_res - omega - com_unit * nu_res) &
                        + vT_res**2.0d0 * a2)&
                        / (om_E_res - omega - com_unit * nu_res) &
                        - 1.0d0 / omc_res * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - kr(k1))**2d0) &
                        - eval_besselI0)) !exp(-eval_bp) * eval_besselI0))

                    !dydr = dydr / lambda_D_res**2d0
                else
                    ! ideal region, k_parallel not near zero
                    !write(*,*) 'NOT large z limit'
                    dydr = dydr + omc_res /(lambda_D_res**2d0 * vT_res) &
                            * exp(com_unit * (kr(k1) - c_krp) * r) * (sqrt(pi)/kp_res * plasma_Z(z0_res) &
                            * (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                            * sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) - sqrt(2d0 * pi) * vT_res / omc_res &
                            * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - kr(k1))**2d0) - eval_besselI0)!exp(-eval_bp) &
                            !* eval_besselI0)
                    !dydr = dydr / lambda_D_res**2d0
                end if

            end do

            call binsrc(c_krp, 1, k_space_dim, krp, ikrp)
            ibeg = max(1, ikrp - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. k_space_dim) then
                iend = k_space_dim
                ibeg = iend - nlagr + 1
            end if
            call plag_coeff(nlagr, nder, c_krp, krp(ibeg:iend), coef)

            varphi = sum(cmplx(coef(0,:)) * varphi_lkr(lp, ibeg:iend))
            dydr = dydr * conjg(varphi)
            write(*,*) 'dydr = ', dydr
            !write(*,*) 'varphi = ', varphi_lkr(lp, ibeg:iend)
            !write(*,*) 'coef = ', cmplx(coef(0,:))
        end subroutine



end module