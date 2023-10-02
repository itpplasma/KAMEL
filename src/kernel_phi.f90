subroutine kernel_phi(write_out)
! fill kernels correlating phi with rho and j
! both are filled simultaneously, since they share the same coefficients

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants
    use kernel, only: K_rho_phi, K_j_phi
    use adaptive_int, only: odeint_c

    implicit none

    logical, intent(in) :: write_out

    double complex :: besselI ! complex bessel function from bessel.f90
    double complex :: plasma_Z ! plasma dispersion function

    double complex, dimension(:,:), allocatable :: K_rho_phi_limit

    integer :: nlagr = 4
    integer :: nder = 0
    double complex :: a0, a1, a2
    !double precision :: eval_bp, eval_bt ! b_+ and b_\times
    double complex :: eval_bp, eval_bt ! b_+ and b_\times

    double complex, dimension(1) :: res

    integer :: sigma ! for loop over species
    integer :: i,j,n
    double precision :: int_fac
    double complex :: integr
    double precision :: c_kr, c_krp
    double precision :: eps = 1d-6
    double precision :: h1 = 0.1d0
    double precision :: hmin = 0.0d0
    integer :: nok, nbad
    double complex, dimension(1) :: temp

    integer :: choose_mode = 2


    if (fstatus == 1) write(*,*) 'Status: Generating kernels rho phi and j phi, write_out = ', write_out
    c_kr = kr(10)
    c_krp = krp(30)
    open(unit=77, file=trim(output_path)//'kernel/integrand_re.dat')
    open(unit=78, file=trim(output_path)//'kernel/integrand_im.dat')
    do i = 1, iprof_length
        call integrand_K_rho_phi_limit(r_prof(i), temp, res)
        write(77,*) r_prof(i), real(res(1))
        write(78,*) r_prof(i), dimag(res(1))
    end do
    close(77)
    close(78)
    !stop
    res = 0.0d0

    allocate(K_rho_phi(k_space_dim, k_space_dim), K_j_phi(k_space_dim, k_space_dim))
    allocate(K_rho_phi_limit(k_space_dim, k_space_dim))

    K_rho_phi = 0.0d0
    K_j_phi = 0.0d0

    if (choose_mode == 1) then ! trapezoidal integration

        ! fill kernel matrix and integrate over radius at the same time.
        ! For the integration, the trapezoidal method is used.
        do i=1, k_space_dim ! kr
            do j = 1, k_space_dim ! krp

                ! check integrand of rg integration
                if (i == 10 .and. j==10) then
                    write(*,*) 'Writing out integrand for testing, kr=', kr(i)
                    open(unit=192, file=trim(output_path)//'kernel/integrand_re.dat')
                    open(unit=193, file=trim(output_path)//'kernel/integrand_im.dat')
                    do n=1, iprof_length
                        eval_bp = vTe(n)**2d0 / (2d0*omce(n)**2d0) * &
                                (2d0*ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0)
                        eval_bt = vTe(n)**2d0 / (omce(n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                            sqrt(ks(n)**2d0 + krp(j)**2d0)

                        a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E(n) / omce(n) + ks(n) * vTe(n)**2d0 &
                            / (omce(n)**2d0) * (A1e(n) + (1.0d0 + eval_bp) * A2e(n))) + ks(n) * vTe(n)**2d0 / (omce(n)**2d0) &
                            * A2e(n)* eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                        a1 = - kp(n)/omce(n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                        a2 = ks(n) / (2d0 * omce(n)**2d0) * A2e(n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                        integr = omce(n) /(lambda_De(n)**2d0 * vTe(n)) &
                                * exp(com_unit * (kr(i) - krp(j)) * r_prof(n)) * (sqrt(pi)/kp(n) * (plasma_Z(z0e(n)) &
                                * (a0 + sqrt(2d0) * vTe(n) * z0e(n) * a1 + vTe(n)**2d0 * a2 * 2d0 * z0e(n)**2d0) + vTe(n) &
                                * sqrt(2d0) * (a1 + vTe(n) * sqrt(2d0) * z0e(n) * a2)) - sqrt(2d0 * pi) * vTe(n) / omce(n) &
                                * (exp(-vTe(n)**2d0/(2d0 * omce(n)**2d0) * (krp(j) - kr(i))**2d0) - exp(-eval_bp) &
                                * besselI(0, eval_bt, 0)))
                        write(192,*) r_prof(n), real(integr)
                        write(193,*) r_prof(n), dimag(integr)
                    end do
                    close(192)
                    close(193)
                end if

                do n = 1, iprof_length ! r_prof

                    ! electrons
                    eval_bp = vTe(n)**2d0 / (2d0*omce(n)**2d0) * &
                                (2d0*ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0)
                    eval_bt = vTe(n)**2d0 / (omce(n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                            sqrt(ks(n)**2d0 + krp(j)**2d0)

                    a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E(n) / omce(n) + ks(n) * vTe(n)**2d0 &
                        / (omce(n)**2d0) * (A1e(n) + (1.0d0 + eval_bp) * A2e(n))) + ks(n) * vTe(n)**2d0 / (omce(n)**2d0) * A2e(n)&
                        * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                    a1 = - kp(n)/omce(n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                    a2 = ks(n) / (2d0 * omce(n)**2d0) * A2e(n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                
                    ! for the trapezoidal integration; first and last summands need factor 1/2
                    if (n==1 .or. n==iprof_length) then
                        int_fac = 0.5d0
                    else
                        int_fac = 1.0d0
                    end if
                
                    K_rho_phi(i,j) = K_rho_phi(i,j) + int_fac * omce(n) /(lambda_De(n)**2d0 * vTe(n)) &
                                    * exp(com_unit * (kr(i) - krp(j)) * r_prof(n)) * (sqrt(pi)/kp(n) * (plasma_Z(z0e(n)) &
                                    * (a0 + sqrt(2d0) * vTe(n) * z0e(n) * a1 + vTe(n)**2d0 * a2 * 2d0 * z0e(n)**2d0) + vTe(n) &
                                    * sqrt(2d0) * (a1 + vTe(n) * sqrt(2d0) * z0e(n) * a2)) - sqrt(2d0 * pi) * vTe(n) / omce(n) &
                                    * (exp(-vTe(n)**2d0/(2d0 * omce(n)**2d0) * (krp(j) - kr(i))**2d0) - exp(-eval_bp) &
                                    * besselI(0, eval_bt, 0)))

                    K_j_phi(i,j) = K_j_phi(i,j) + int_fac * omce(n) / (lambda_De(n)**2d0 * kp(n)) &
                                    * exp(com_unit * (kr(i) - krp(j)) &
                                    * r_prof(n))* ((z0e(n) * plasma_Z(z0e(n)) + 1d0) * (a0 + sqrt(2d0) * vTe(n) * z0e(n) * a1 &
                                    + 2d0 * vTe(n) * z0e(n)**2d0 * a2) + 2d0 * vTe(n) * a2 + ks(n) * vTe(n)**2d0 &
                                    / (2d0 * omce(n)**2d0) * A2e(n) * (exp(- vTe(n)**2d0 / (2d0 * omce(n)**2d0) &
                                    * (krp(j)-kr(i))**2d0) - exp(- eval_bp) * besselI(0, eval_bt, 0))) 

                    ! ions
                    do sigma=1, ispecies

                        eval_bp = vTi(sigma, n)**2d0 / (2d0*omci(sigma, n)**2d0) * &
                                    (2d0*ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0)
                        eval_bt = vTi(sigma, n)**2d0 / (omci(sigma, n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                            sqrt(ks(n)**2d0 + krp(j)**2d0)

                        a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E(n) / omci(sigma, n) + ks(n) * vTi(sigma, n)**2d0&
                            / omci(sigma, n)**2d0 * (A1i(sigma, n) + (1d0 + eval_bp) * A2i(sigma, n))) + ks(n) &
                            * vTi(sigma, n)**2d0 &
                            / omci(sigma, n)**2d0 * A2i(sigma, n) * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                        a1 = - kp(n)/omci(sigma, n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                        a2 = ks(n) / (2d0 * omci(sigma, n)**2d0) * A2i(sigma, n) * exp(-eval_bp) * besselI(0,eval_bt,0)

                        ! for the trapezoidal integration; first and last summands need factor 1/2
                        if (n==1 .or. n==iprof_length) then
                            int_fac = 0.5d0
                        else
                            int_fac = 1.0d0
                        end if

                        K_rho_phi(i,j) = K_rho_phi(i,j) + int_fac * omci(sigma, n) /(lambda_Di(sigma, n)**2d0 * vTi(sigma, n)) &
                                        * exp(com_unit * (kr(i) -krp(j)) * r_prof(n)) * ((sqrt(pi)/kp(n) &
                                        * (plasma_Z(z0i(sigma, n))&
                                        * (a0 + sqrt(2d0) * vTi(sigma, n) * z0i(sigma, n) * a1 + vTi(sigma, n)**2d0 * a2 * 2d0 &
                                        * z0i(sigma, n)**2d0) + vTi(sigma, n) * sqrt(2d0) * (a1 + vTi(sigma, n) * sqrt(2d0) &
                                        * z0i(sigma, n) * a2))) - sqrt(2d0) * vTi(sigma, n) / omci(sigma, n) * &
                                        (exp(-vTi(sigma, n)**2d0 / (2d0 * omci(sigma, n)**2d0) * (krp(j) - kr(i))**2d0) &
                                        - exp(-eval_bp) * besselI(0, eval_bt, 0)))

                        K_j_phi(i,j) = K_j_phi(i,j) + int_fac * omci(sigma, n) / (lambda_Di(sigma, n)**2d0 * kp(n)) &
                                        * exp(com_unit * (kr(i) - krp(j)) * r_prof(n))* ((z0i(sigma, n) * plasma_Z(z0i(sigma, n))&
                                        + 1d0) * (a0 + sqrt(2d0) * vTi(sigma, n) * z0i(sigma, n) * a1 &
                                        + 2d0 * vTi(sigma, n) * z0i(sigma, n)**2d0 * a2) + 2d0 * vTi(sigma, n) * a2 &
                                        + ks(n) * vTi(sigma, n)**2d0 / (2d0 * omci(sigma, n)**2d0) * A2i(sigma, n) &
                                        * (exp(- vTi(sigma, n)**2d0 / (2d0 * omci(sigma, n)**2d0) * (krp(j)-kr(i))**2d0) & 
                                        - exp(- eval_bp) * besselI(0,eval_bt,0))) 

                    end do
                end do
            end do
        end do 

        K_rho_phi = K_rho_phi / (2d0**(7d0/2d0) * pi**(5d0/2d0)) * (r_prof(size(r_prof)) - r_prof(1))/iprof_length
    
    else if(choose_mode == 2) then ! odeint method from NR

        do i=1, k_space_dim ! kr
            do j = 1, k_space_dim ! krp
                res = 0.0d0
                c_kr = kr(i)
                c_krp = krp(j)
                !call odeint_c(TODO)
                !call odeint_c(res, r_prof(1), r_prof(10), eps, h1, hmin, nok, nbad, integrand_K_rho_phi)
                call odeint_c(res, r_prof(1), r_prof(iprof_length), eps, h1, hmin, nok, nbad, integrand_K_rho_phi_limit)
                !K_rho_phi(i,j) = res(1)
                K_rho_phi(i,j) = res(1)
            end do
            if (fstatus == 1) write(*,*) '    ', dble(i) * 100.0d0/dble(k_space_dim), '% of kernel phi filled'
        end do

        K_rho_phi = K_rho_phi / (2.0d0**3.0d0 * pi**2.0d0)

    end if

    !K_rho_phi = K_rho_phi / (2d0**(7d0/2d0) * pi**(5d0/2d0)) * (r_prof(size(r_prof)) - r_prof(1))/iprof_length
    !K_rho_phi = K_rho_phi / (2d0**(3d0) * pi**(2d0)) * (r_prof(size(r_prof)) - r_prof(1))/iprof_length
    
    if (write_out) call write_kernel_phi
    deallocate(K_rho_phi_limit)

    contains

        subroutine write_kernel_phi
            implicit none
            integer :: i,j
            logical :: ex

            inquire(file=trim(output_path)//'kernel', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'kernel')
            end if
            open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_kr_re.dat')
            open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_kr_im.dat')
            open(unit=79, file=trim(output_path)//'kernel/K_j_phi_kr_re.dat')
            open(unit=80, file=trim(output_path)//'kernel/K_j_phi_kr_im.dat')
            !open(unit=81, file=trim(output_path)//'kernel/K_rho_phi_kr_limit_re.dat')
            !open(unit=82, file=trim(output_path)//'kernel/K_rho_phi_kr_limit_im.dat')
            do i=1,k_space_dim
                do j=1,k_space_dim
                    write(77,*) real(K_rho_phi(i,j))
                    write(78,*) dimag(K_rho_phi(i,j))
                    write(79,*) real(K_j_phi(i,j))
                    write(80,*) dimag(K_j_phi(i,j))
                    write(81,*) real(K_rho_phi_limit(i,j))
                    write(82,*) dimag(K_rho_phi_limit(i,j))
                end do
            end do
            close(77)
            close(78)
            close(79)
            close(80)
            close(81)
            close(82)

        end subroutine

        ! used for the odeint RK integrator. Needs "derivs" subroutine that returns the right hand side of
        ! dy/dr = f(r,y)
        ! In this case, the right hand side is the integrand of the kernel
        ! y is just a dummy variable and is not used
        subroutine integrand_K_rho_phi(r, y, dydr)

            implicit none

            double precision, intent(in) :: r
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr

            integer :: sigma ! for loop over species
            integer :: ibeg, iend
            integer :: ir

            double precision, dimension(:,:), allocatable :: coef

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res
            double complex :: z0_res

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            call binsrc(r_prof, 1, iprof_length, r, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, r, r_prof(ibeg:iend), coef)


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
                        + c_kr**2.0d0 + c_krp**2.0d0)
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + c_kr**2.0d0)&
                        * sqrt(ks_res**2.0d0 + c_krp**2.0d0)


                a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                    / (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) * A2_res &
                    * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                a1 = - kp_res/omc_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                
                
                dydr = omc_res /(lambda_D_res**2d0 * vT_res) &
                                * exp(com_unit * (c_kr - c_krp) * r) * (sqrt(pi)/kp_res * plasma_Z(z0_res) &
                                * (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                                * sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) - sqrt(2d0 * pi) * vT_res / omc_res &
                                * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - c_kr)**2d0) - exp(-eval_bp) &
                                * besselI(0, eval_bt, 0))
            end do

            deallocate(coef)

        end subroutine

        subroutine integrand_K_rho_phi_limit(r, y, dydr)

            use setup, only: omega
            implicit none

            double precision, intent(in) :: r
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr

            integer :: sigma ! for loop over species
            integer :: ibeg, iend
            integer :: ir

            double precision, dimension(:,:), allocatable :: coef

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res, nu_res
            double complex :: z0_res

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            call binsrc(r_prof, 1, iprof_length, r, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, r, r_prof(ibeg:iend), coef)


            ks_res = sum(coef(0,:) * ks(ibeg:iend))
            kp_res = sum(coef(0,:) * kp(ibeg:iend))
            om_E_res = sum(coef(0,:) * om_E(ibeg:iend))

            dydr = 0.0d0

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
                        + c_kr**2.0d0 + c_krp**2.0d0)
                eval_bt = vT_res**2.0d0 /(omc_res**2.0d0) * sqrt(ks_res**2.0d0 + c_kr**2.0d0)&
                        * sqrt(ks_res**2.0d0 + c_krp**2.0d0)


                a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E_res / omc_res + ks_res * vT_res**2d0 &
                    / (omc_res**2d0) * (A1_res + (1.0d0 + eval_bp) * A2_res)) + ks_res * vT_res**2d0 / (omc_res**2d0) * A2_res &
                    * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                a1 = - kp_res/omc_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                a2 = ks_res / (2d0 * omc_res**2d0) * A2_res * exp(-eval_bp) * besselI(0,eval_bt,0)
                
                
                !dydr = omc_res /(lambda_D_res**2d0 * vT_res) &
                !                * exp(com_unit * (c_kr - c_krp) * r) * (sqrt(pi)/kp_res * plasma_Z(z0_res) &
                !                * (a0 + sqrt(2.0d0) * vT_res * z0_res * a1 + vT_res**2d0 * a2 * 2d0 * z0_res**2d0) + vT_res &
                !                * sqrt(2d0) * (a1 + vT_res * sqrt(2d0) * z0_res * a2)) - sqrt(2d0 * pi) * vT_res / omc_res &
                !                * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - c_kr)**2d0) - exp(-eval_bp) &
                !                * besselI(0, eval_bt, 0))

                dydr = dydr + omc_res / lambda_D_res**2.0d0 * exp(com_unit * (c_kr - c_krp) * r) &
                    * ((a0 * (1.0d0 + (2.0d0 * vT_res**2.0d0 *kp_res**2.0d0)/ (om_E_res - omega - com_unit * nu_res)**2.0d0) &
                    - (vT_res**2.0d0 * kp_res * a1)/ (om_E_res - omega - com_unit * nu_res) &
                    + vT_res**2.0d0 * a2)&
                       / (om_E_res - omega - com_unit * nu_res) &
                       - 1.0d0 / omc_res * (exp(-vT_res**2d0/(2d0 * omc_res**2d0) * (c_krp - c_kr)**2d0) &
                       - exp(-eval_bp) * besselI(0, eval_bt, 0)))
            end do

            deallocate(coef)

        end subroutine

end subroutine

