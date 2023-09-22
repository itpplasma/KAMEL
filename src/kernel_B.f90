subroutine kernel_B(write_out)

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants
    use kernel, only: K_rho_B, K_j_B

    implicit none

    logical, intent(in) :: write_out

    double complex :: besselI
    double complex :: plasma_Z

    double complex :: eval_besselI0 = 0.0d0
    double complex :: eval_besselIm1 = 0.0d0

    !double precision :: eval_bp, eval_bt
    double complex :: eval_bp, eval_bt

    integer :: i, j, n, sigma
    double precision :: int_fac

    if (fstatus == 1) write(*,*) 'Status: Generating kernel rho B and j B, write_out=',write_out

    allocate(K_rho_B(k_space_dim, k_space_dim))
    allocate(K_j_B(k_space_dim, k_space_dim))

    K_rho_B = 0.0d0
    K_j_B = 0.0d0

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
                                    * plasma_Z(z0i(sigma, n)) + 1) * ((A1i(sigma, n) + A2i(sigma, n) * (1d0 + eval_bp &
                                    + z0i(sigma, n)**2)) * eval_besselI0 + A2i(sigma, n) * eval_bt &
                                    * eval_besselIm1) + 0.5d0 * A2i(sigma, n) * eval_besselI0)

                    K_j_B(i,j) = K_j_B(i,j) + int_fac * vTi(sigma, n)**3d0 / (lambda_Di(sigma, n)**2d0 * omci(sigma, n) * kp(n)) &
                                    * exp(com_unit * (kr(i) - krp(j)) * r_prof(n)) * exp(-eval_bp) * z0i(sigma, n) &
                                    * ((z0i(sigma, n) * plasma_Z(z0i(sigma, n)) + 1d0) * ((A1i(sigma, n) + A2i(sigma, n) &
                                    * (1d0 + eval_bp + z0i(sigma, n)**2d0)) * eval_besselI0 + A2i(sigma, n) * eval_bt &
                                    * eval_besselIm1) + 0.5d0 * A2i(sigma, n) * eval_besselI0)
                end do
            end do
        end do
    end do

    K_rho_B = K_rho_B / (8d0 * pi**2d0 * sol)
    K_j_B = K_j_B * com_unit / (2d0**(5d0/2d0) * pi**2d0 * sol)

    if (write_out) call write_kernel

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
            do i=1,k_space_dim
                do j=1,k_space_dim
                    write(77,*) real(K_rho_B(i,j))
                    write(78,*) dimag(K_rho_B(i,j))
                    write(79,*) real(K_j_B(i,j))
                    write(80,*) dimag(K_j_B(i,j))
                end do
            end do
            close(77)
            close(78)
            close(79)
            close(80)

        end subroutine


end subroutine