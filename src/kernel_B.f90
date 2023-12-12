subroutine kernel_B(write_out)

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants
    use kernel, only: K_rho_B, K_j_B
    !use adaptive_int, only: odeint_c
    use omp_lib

    implicit none

    double complex :: besselI
    double complex :: plasma_Z
    logical, intent(in) :: write_out

    double complex :: eval_besselI0 = 0.0d0
    double complex :: eval_besselIm1 = 0.0d0

    !double precision :: eval_bp, eval_bt
    double complex :: eval_bp, eval_bt

    integer :: max_threads
    integer :: i, j, n, sigma
    double precision :: int_fac
    integer :: choose_mode = 2
    double precision :: c_kr, c_krp
    double precision :: eps = 1d-4
    double precision :: h1 = 0.1d0
    double precision :: hmin = 0.0d0
    integer :: nok, nbad
    integer :: nlagr = 4
    integer :: nder = 0
    double complex, dimension(1) :: temp

    double complex, dimension(1) :: res

    max_threads = OMP_GET_MAX_THREADS()
    if (fdebug == 1) write(*,*) ' Debug: Number of threads = ', max_threads
    if (fstatus == 1) write(*,*) 'Status: Generating kernel rho B and j B, write_out=',write_out

    c_kr = 0.0d0 !kr(10)
    c_krp = 0.0d0 !krp(30)
    open(unit=77, file=trim(output_path)//'kernel/integrand_B_re.dat')
    open(unit=78, file=trim(output_path)//'kernel/integrand_B_im.dat')
    do i = 1, iprof_length
        !if (r_prof(i) < int_rmin) cycle
        call integrand_K_rho_B_stitching(r_prof(i), temp, res, c_kr, c_krp)
        write(77,*) r_prof(i), real(res(1))
        write(78,*) r_prof(i), dimag(res(1))
    end do
    close(77)
    close(78)
    stop
    allocate(K_rho_B(k_space_dim, k_space_dim))
    allocate(K_j_B(k_space_dim, k_space_dim))

    K_rho_B = 0.0d0
    K_j_B = 0.0d0

    if (choose_mode == 1) then
        do i=1, k_space_dim ! kr
            do j=1, k_space_dim ! kr prime
                do n=1, iprof_length ! r_prof

                    ! electrons
                    eval_bp = vTe(n)**2d0 / (2d0 * omce(n)**2d0) * (2d0 * ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0) 
                    eval_bt = vTe(n)**2d0 / (omce(n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                            sqrt(ks(n)**2d0 + krp(j)**2d0)

                    eval_besselI0 = besselI(0, eval_bt, 0)
                    eval_besselIm1 = besselI(-1, eval_bt, 0)

                    ! for trapezoidal integration
                    if (n==1 .or. n==iprof_length) then
                        int_fac = 0.5d0
                    else
                        int_fac = 1.0d0
                    end if

                    K_rho_B(i,j) = K_rho_B(i,j) + int_fac * exp(com_unit * (kr(i) - krp(j))*r_prof(n)) *vTe(n)**2d0 &
                                / (lambda_De(n)**2d0 * omce(n) * kp(n)) * exp(-eval_bp) * ((z0e(n) * plasma_Z(z0e(n)) + 1d0) &
                                * ((A1e(n) + A2e(n) * (1d0 + eval_bp + z0e(n)**2d0)) * eval_besselI0 + A2e(n) * eval_bt &
                                * eval_besselIm1) + 0.5d0 * A2e(n) * eval_besselI0)

                
                    K_j_B(i,j) = K_j_B(i,j) + int_fac * vTe(n)**3d0 / (lambda_De(n)**2d0 * omce(n) * kp(n)) &
                                * exp(com_unit * (kr(i) - krp(j)) * r_prof(n)) * exp(-eval_bp) * z0e(n) &
                                * ((z0e(n) * plasma_Z(z0e(n)) + 1d0) * ((A1e(n) + A2e(n) * (1d0 + eval_bp + z0e(n)**2d0)) &
                                * eval_besselI0 + A2e(n) * eval_bt * eval_besselIm1) + 0.5d0 * A2e(n) &
                                * eval_besselI0)

                    ! ions
                    do sigma=1, ispecies

                        eval_bp = vTi(sigma, n)**2d0 / (2d0 * omci(sigma, n)**2d0) * (2d0 * ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0) 
                        eval_bt = vTi(sigma, n)**2d0 / (omci(sigma, n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                                sqrt(ks(n)**2d0 + krp(j)**2d0)

                        eval_besselI0 = besselI(0, eval_bt, 0)
                        eval_besselIm1 = besselI(-1, eval_bt, 0)

                        K_rho_B(i,j) = K_rho_B(i,j) + int_fac * exp(com_unit * (kr(i) - krp(j))*r_prof(n))* vTi(sigma, n)**2d0 &
                                    / (lambda_Di(sigma, n)**2d0 * kp(n) * omci(sigma, n)) * exp(-eval_bp) * ((z0i(sigma, n) &
                                    * plasma_Z(z0i(sigma, n)) + 1d0) * ((A1i(sigma, n) + A2i(sigma, n) * (1d0 + eval_bp &
                                    + z0i(sigma, n)**2d0)) * eval_besselI0 + A2i(sigma, n) * eval_bt &
                                    * eval_besselIm1) + 0.5d0 * A2i(sigma, n) * eval_besselI0)

                        K_j_B(i,j) = K_j_B(i,j) + int_fac * vTi(sigma, n)**3d0 / (lambda_Di(sigma, n)**2d0 * omci(sigma, n) &
                                    * kp(n))&
                                    * exp(com_unit * (kr(i) - krp(j)) * r_prof(n)) * exp(-eval_bp) * z0i(sigma, n) &
                                    * ((z0i(sigma, n) * plasma_Z(z0i(sigma, n)) + 1d0) * ((A1i(sigma, n) + A2i(sigma, n) &
                                    * (1d0 + eval_bp + z0i(sigma, n)**2d0)) * eval_besselI0 + A2i(sigma, n) * eval_bt &
                                    * eval_besselIm1) + 0.5d0 * A2i(sigma, n) * eval_besselI0)
                    end do
                end do
            end do
        end do
        K_rho_B = K_rho_B * com_unit / (8d0 * pi**2d0 * sol)
        K_j_B = K_j_B * com_unit / (2d0**(5d0/2d0) * pi**2d0 * sol)

    else if(choose_mode == 2)then

        !$OMP PARALLEL DO collapse(2) default(none) schedule(guided) &
        !$OMP PRIVATE(res, nok, nbad, i,j) &
        !$OMP SHARED(c_kr, c_krp, K_rho_B, kr, krp, r_prof, iprof_length, k_space_dim, eps,&
        !$OMP h1, hmin, fstatus, max_threads) 
        do i=1, k_space_dim ! kr
            do j = 1, k_space_dim ! krp
                res = 0.0d0
                call integrate_rg_c(res, size(res), kr(i), krp(j), r_prof(1), r_prof(iprof_length), eps, h1,&
                                    hmin, nok, nbad, integrand_K_rho_B_stitching)
                K_rho_B(i,j) = res(1)
            end do
            !!$OMP critical
            !if (fstatus == 1 .and. OMP_GET_THREAD_NUM() == 0) then
            !    write(*,*) '    ', dble(i) * 100.0d0/dble(k_space_dim) * max_threads, '% of kernel B filled'
            !end if
            !!$OMP end critical
        end do
        !$OMP END PARALLEL DO

        K_rho_B = K_rho_B  * (-1.0d0) * com_unit /(8d0 * pi**2.0d0 * sol)
    end if

    
    if (write_out) call write_kernel
    !deallocate(K_rho_B_limit)

    contains

        subroutine write_kernel

            implicit none
            integer :: i,j
            logical :: ex

            inquire(file=trim(output_path)//'kernel', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'kernel')
            end if
            open(unit=77, file=trim(output_path)//'kernel/K_rho_B_kr_im.dat')
            open(unit=78, file=trim(output_path)//'kernel/K_rho_B_kr_re.dat')
            open(unit=79, file=trim(output_path)//'kernel/K_j_B_kr_im.dat')
            open(unit=80, file=trim(output_path)//'kernel/K_j_B_kr_re.dat')
            !open(unit=81, file=trim(output_path)//'kernel/K_rho_B_kr_limit_re.dat')
            !open(unit=82, file=trim(output_path)//'kernel/K_rho_B_kr_limit_im.dat')
            do i=1,k_space_dim
                do j=1,k_space_dim
                    write(77,*) real(K_rho_B(i,j))
                    write(78,*) dimag(K_rho_B(i,j))
                    write(79,*) real(K_j_B(i,j))
                    write(80,*) dimag(K_j_B(i,j))
                    !write(81,*) real(K_rho_B_limit(i,j))
                    !write(82,*) dimag(K_rho_B_limit(i,j))
                end do
            end do
            close(77)
            close(78)
            close(79)
            close(80)
            !close(81)
            !close(82)

        end subroutine

        subroutine integrand_K_rho_B_limit(r, y, dydr)

            use setup, only: omega
            implicit none

            double complex :: besselI
            double complex :: plasma_Z

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
            double complex :: eval_bp, eval_bt

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
                
                dydr = dydr + vT_res**4.0d0 * kp_res / (lambda_D_res**2.0d0 * omc_res) * exp(com_unit * (c_kr - c_krp) * r) &
                      * exp(-eval_bp) &
                      / (om_E_res - omega - com_unit * nu_res)**2.0d0 * ((A1_res + A2_res * (1 + eval_bp)) &
                      * besselI(0,eval_bt,0) + A2_res * eval_bt * besselI(-1, eval_bt, 0))

            end do

            deallocate(coef)

        end subroutine

        subroutine integrand_K_rho_B_limit_kr(r, y, dydr, c_kr, c_krp)

            use setup, only: omega
            implicit none

            double complex :: besselI
            double complex :: plasma_Z

            double precision, intent(in) :: r, c_kr, c_krp
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr

            integer :: sigma ! for loop over species
            integer :: ibeg, iend
            integer :: ir

            double precision, dimension(:,:), allocatable :: coef

            double precision :: vT_res, omc_res, ks_res, om_E_res, kp_res, &
                                A1_res, A2_res, lambda_D_res, nu_res
            double complex :: z0_res
            double complex :: eval_bp, eval_bt

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
                
                dydr = dydr + vT_res**4.0d0 * kp_res / (lambda_D_res**2.0d0 * omc_res) * exp(com_unit * (c_kr - c_krp) * r) &
                      * exp(-eval_bp) &
                      / (om_E_res - omega - com_unit * nu_res)**2.0d0 * ((A1_res + A2_res * (1 + eval_bp)) &
                      * besselI(0,eval_bt,0) + A2_res * eval_bt * besselI(-1, eval_bt, 0))

            end do

            deallocate(coef)

        end subroutine

        subroutine integrand_K_rho_B_stitching(r, y, dydr, c_kr, c_krp)

            use setup, only: omega
            implicit none

            double complex :: besselI
            double complex :: plasma_Z

            double precision, intent(in) :: r, c_kr, c_krp
            double complex, dimension(:), intent(in) :: y
            double complex, dimension(:), intent(out) :: dydr
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
                
                if (real(eval_bt) > 5d0) then
                    ! limit close to magnetic axis (k_s -> infinity) and large k_r and k_rp
                    eval_besselI0 = exp(eval_bt - eval_bp) /(sqrt(2.0d0 * pi * eval_bt))
                    eval_besselIm1 = exp(- eval_bp + asinh(-1.0d0/eval_bt) + eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)) &
                            / (sqrt(2.0d0*pi*eval_bt * sqrt(1.0d0 + 1.0d0/eval_bt**2.0d0)))
                    !write(*,*) 'eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1 
                else
                    eval_besselI0 = besselI(0, eval_bt, 0) * exp(-eval_bp)
                    eval_besselIm1 = besselI(-1, eval_bt, 0) * exp(-eval_bp)
                end if

                if (abs(z0_res) > 3d0) then
                    ! limit close to resonant surface, z -> infinity, k_parallel -> 0
                    dydr = dydr - vT_res**4.0d0 * kp_res / (lambda_D_res**2.0d0 * omc_res) * exp(com_unit * (c_kr - c_krp) * r) &
                        !* exp(-eval_bp) &
                        / (om_E_res - omega - com_unit * nu_res)**2.0d0 * ((A1_res + A2_res * (1.0d0 + eval_bp)) &
                        * eval_besselI0 + A2_res * eval_bt * eval_besselIm1)
                    write(*,*) 'spec = ', sigma, ', r = ', r, ', lambda_D = ', lambda_D_res, &!', omc_res = ', omc_res!eval_besselI0 = ', eval_besselI0, ', eval_besselIm1 = ', eval_besselIm1!', 
                    ', vT = ', vT_res
                    !'eval_bt = ', eval_bt, ', eval_bp = ', eval_bp
                    !, k_p = ', kp_res!, ', omega_E = ', om_E_res
                else
                    ! ideal region
                    dydr = dydr + vT_res**2.0d0 / (lambda_D_res**2.0d0 * omc_res * kp_res) &
                          * exp(com_unit * (c_kr - c_krp) * r) * ((z0_res * plasma_Z(z0_res) + 1.0d0) &
                          * ((A1_res + A2_res * (1 + eval_bp + z0_res**2.0d0)) * eval_besselI0 &
                          + A2_res * eval_bp * eval_besselIm1) + 0.5d0 * A2_res * eval_besselI0)
                end if

            end do

            deallocate(coef)

        end subroutine

end subroutine