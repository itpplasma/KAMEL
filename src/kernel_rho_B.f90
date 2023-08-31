subroutine kernel_rho_B(write_out)

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants

    implicit none

    logical, intent(in) :: write_out

    double complex :: besselI
    double complex :: plasma_Z

    double precision, dimension(:,:,:), allocatable :: bp, bt
    double complex, dimension(:,:), allocatable :: K_rho_B

    integer :: i, j, n, sigma
    double precision :: int_fac


    if (fstatus == 1) write(*,*) 'Status: Generating kernel rho phi'

    allocate(bp(k_space_dim, k_space_dim, iprof_length), bt(k_space_dim, k_space_dim, iprof_length))
    allocate(K_rho_B(k_space_dim, k_space_dim))

    do i=1, k_space_dim ! kr
        do j=1, k_space_dim ! kr prime
            do n=1, iprof_length ! r_prof
                ! electrons
                bp(i,j,n) = vTe(n)**2 / (2*omce(n)**2) * &
                            (2*ks(n)**2 + kr(i)**2 + krp(j)**2) 
                bt(i,j,n) = vTe(n)**2 / (omce(n)**2) * sqrt(ks(n)**2 + kr(i)**2) *&
                           sqrt(ks(n)**2 + krp(j)**2)

                if (n==1 .or. n==iprof_length) then
                    int_fac = 0.5d0
                else
                    int_fac = 1.0d0
                end if

                K_rho_B(i,j) = K_rho_B(i,j) + int_fac * exp(com_unit * (kr(i) - krp(j))*r_prof(n)) *&
                                vTe(n)**2 /(lambda_De(n)**2 * &
                                omce(n)*kp(n)) * exp(-bp(i,j,n)) * ((z0e(n) * plasma_Z(z0e(n)) &
                                + 1) * ((A1e(n) + A2e(n) * (1+bp(i,j,n) + z0e(n)**2)) * &
                                besselI(0,bt(i,j,n),0) + A2e(n) * bt(i,j,n) * &
                                besselI(-1, bt(i,j,n), 0)) + 0.5d0 * A2e(n) * besselI(0,bt(i,j,n),0))


                do sigma=1, ispecies
                    bp(i,j,n) = vTi(sigma, n)**2 / (2*omci(sigma, n)**2) * &
                            (2*ks(n)**2 + kr(i)**2 + krp(j)**2) 
                    bt(i,j,n) = vTi(sigma, n)**2 / (omci(sigma, n)**2) * sqrt(ks(n)**2 + kr(i)**2) *&
                           sqrt(ks(n)**2 + krp(j)**2)


                    K_rho_B(i,j) = K_rho_B(i,j) + int_fac * exp(com_unit * (kr(i) - krp(j))*r_prof(n))*&
                                 vTi(sigma, n)**2 / &
                                (lambda_Di(sigma, n)**2 * kp(n) *&
                                omci(sigma, n)) * exp(-bp(i,j,n)) * ((z0i(sigma, n) * &
                                plasma_Z(z0i(sigma, n)) + 1) * ((A1i(sigma, n) + A2i(sigma, n) * &
                                (1+bp(i,j,n) + z0i(sigma, n)**2)) * besselI(0,bt(i,j,n),0) + &
                                A2i(sigma, n) * bt(i,j,n) * besselI(-1, bt(i,j,n), 0)) + 0.5d0 * A2i(sigma, n) * &
                                besselI(0,bt(i,j,n),0))
                end do

            end do
        end do
    end do

    K_rho_B = K_rho_B / (8d0 * pi**2 * sol)


    if (write_out) call write_kernel

    contains

        subroutine write_kernel

            implicit none
            integer :: i,j

            open(unit=77, file=trim(output_path)//'kernel/K_rho_B_im.dat')
            open(unit=78, file=trim(output_path)//'kernel/K_rho_B_re.dat')
            do i=1,k_space_dim
                do j=1,k_space_dim
                    write(77,*) real(K_rho_B(i,j))
                    write(78,*) dimag(K_rho_B(i,j))
                end do
            end do
            close(77)
            close(78)

        end subroutine


end subroutine