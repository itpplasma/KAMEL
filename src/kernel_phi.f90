subroutine kernel_phi(write_out)
! fill kernels correlating phi with rho and j
! both are filled simultaneously, since they share the same coefficients

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants
    use kernel, only: K_rho_phi, K_j_phi

    implicit none

    logical, intent(in) :: write_out

    double complex :: besselI ! complex bessel function from bessel.f90
    double complex :: plasma_Z ! plasma dispersion function

    integer :: nder
    double precision :: a0, a1, a2
    double precision :: eval_bp, eval_bt ! b_+ and b_\times

    integer :: sigma ! for loop over species
    integer :: i,j,n
    double precision :: int_fac

    if (fstatus == 1) write(*,*) 'Status: Generating kernels rho phi and j phi, write_out=', write_out

    npoi_der = 4 ! number of polynomials for derivative
    nder = 1 ! first derivative

    allocate(K_rho_phi(k_space_dim, k_space_dim), K_j_phi(k_space_dim, k_space_dim))

    K_rho_phi = 0.0d0
    K_j_phi = 0.0d0

    ! fill kernel matrix and integrate over radius at the same time.
    ! For the integration, the trapezoidal method is used.
    do i=1, k_space_dim ! kr
        do j = 1, k_space_dim ! krp
            do n = 1, iprof_length ! r_prof

                ! electrons
                eval_bp = vTe(n)**2d0 / (2d0*omce(n)**2d0) * &
                            (2d0*ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0)
                eval_bt = vTe(n)**2d0 / (omce(n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                           sqrt(ks(n)**2d0 + krp(j)**2d0)

                a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E(n) / omce(n) + ks(n) * vTe(n)**2d0 &
                    / omce(n)**2d0 * (A1e(n) + (1 + eval_bp) * A2e(n))) + ks(n) * vTe(n)**2d0 / omce(n)**2d0 * A2e(n) &
                    * eval_bt * exp(-eval_bp) * besselI(-1, eval_bt, 0)
                a1 = - kp(n)/omce(n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                a2 = ks(n) / (2d0 * omce(n)**2d0) * A2e(n) * exp(-eval_bp) * besselI(0,eval_bt,0)
                
                ! for the trapezoidal integration; first and last summands need factor 1/2
                if (n==1 .or. n==iprof_length) then
                    int_fac = 0.5d0
                else
                    int_fac = 1.0d0
                end if
                
                K_rho_phi(i,j) = K_rho_phi(i,j) + int_fac * omce(n) /(lambda_De(n)**2d0 * &
                                vTe(n)) * exp(com_unit * (kr(i) - krp(j)) * r_prof(n)) * &
                                ((sqrt(pi)/kp(n) * (plasma_Z(z0e(n)) * (a0 + sqrt(2d0) * &
                                 vTe(n) * z0e(n) * a1 + vTe(n)**2 * a2 * 2 * z0e(n)**2) +&
                                 vTe(n) * sqrt(2d0) * (a1 + vTe(n) * sqrt(2d0) * z0e(n) * a2))) &
                                - sqrt(2d0 * pi) * vTe(n) / omce(n) &
                                * (exp(-vTe(n)**2d0/(2d0 * omce(n)**2d0) * (krp(j) - kr(i))**2) &
                                - exp(-eval_bp) * besselI(0, eval_bt, 0)))

                K_j_phi(i,j) = K_j_phi(i,j) + int_fac * omce(n) / (lambda_De(n)**2d0 * kp(n)) * exp(com_unit * (kr(i) - krp(j)) &
                                * r_prof(n))* ((z0e(n) * plasma_Z(z0e(n)) + 1) * (a0 &
                                + sqrt(2d0) * vTe(n) * z0e(n) * a1 + 2d0 * vTe(n) * z0e(n)**2d0 * a2) &
                                + 2d0 * vTe(n) * a2 + ks(n) * vTe(n)**2d0 / (2d0 * omce(n)**2d0) * A2e(n) &
                                * (exp(- vTe(n)**2d0 / (2d0 * omce(n)**2d0) * (krp(j)-kr(i))**2d0) & 
                                - exp(- eval_bp) * besselI(0, eval_bt, 0))) 


                ! ions
                do sigma=1, ispecies

                    eval_bp = vTi(sigma, n)**2d0 / (2d0*omci(sigma, n)**2d0) * &
                                (2d0*ks(n)**2d0 + kr(i)**2d0 + krp(j)**2d0)
                    eval_bt = vTi(sigma, n)**2d0 / (omci(sigma, n)**2d0) * sqrt(ks(n)**2d0 + kr(i)**2d0) *&
                           sqrt(ks(n)**2d0 + krp(j)**2d0)

                    a0 = exp(-eval_bp) * besselI(0, eval_bt, 0) * (- om_E(n) / omci(sigma, n) + ks(n) * vTi(sigma, n)**2d0&
                        / omci(sigma, n)**2d0 * (A1i(sigma, n) + (1 + eval_bp) * A2i(sigma, n))) + ks(n) * vTi(sigma, n)**2d0&
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
                                    * exp(com_unit * (kr(i) -krp(j)) * r_prof(n)) * ((sqrt(pi)/kp(n) * (plasma_Z(z0i(sigma, n)) &
                                    * (a0 + sqrt(2d0) * vTi(sigma, n) * z0i(sigma, n) * a1 + vTi(sigma, n)**2 * a2 * 2 &
                                    * z0i(sigma, n)**2) + vTi(sigma, n) * sqrt(2d0) * (a1 + vTi(sigma, n) * sqrt(2d0) &
                                    * z0i(sigma, n) * a2))) - sqrt(2d0) * vTi(sigma, n) / omci(sigma, n) * &
                                    (exp(-vTi(sigma, n)**2/(2*omci(sigma, n)**2) * (krp(j) - kr(i))**2) &
                                    - exp(-eval_bp) * besselI(0, eval_bt, 0)))

                    K_j_phi(i,j) = K_j_phi(i,j) + int_fac * omci(sigma, n) / (lambda_Di(sigma, n)**2d0 * kp(n)) &
                                    * exp(com_unit * (kr(i) - krp(j)) * r_prof(n))* ((z0i(sigma, n) * plasma_Z(z0i(sigma, n)) &
                                    + 1) * (a0 + sqrt(2d0) * vTi(sigma, n) * z0i(sigma, n) * a1 &
                                    + 2d0 * vTi(sigma, n) * z0i(sigma, n)**2d0 * a2) + 2d0 * vTi(sigma, n) * a2 &
                                    + ks(n) * vTi(sigma, n)**2d0 / (2d0 * omci(sigma, n)**2d0) * A2i(sigma, n) &
                                    * (exp(- vTi(sigma, n)**2d0 / (2d0 * omci(sigma, n)**2d0) * (krp(j)-kr(i))**2d0) & 
                                    - exp(- eval_bp) * besselI(0,eval_bt,0))) 

                end do
            end do
        end do
    end do 

    K_rho_phi = K_rho_phi / (2**(7d0/2d0) * pi**(5d0/2d0)) * (r_prof(size(r_prof)) - r_prof(1))/iprof_length
    K_rho_phi = K_rho_phi / (2**(3d0) * pi**(2d0)) * (r_prof(size(r_prof)) - r_prof(1))/iprof_length

    if (write_out) call write_kernel_phi

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
            do i=1,k_space_dim
                do j=1,k_space_dim
                    write(77,*) real(K_rho_phi(i,j))
                    write(78,*) dimag(K_rho_phi(i,j))
                    write(79,*) real(K_j_phi(i,j))
                    write(80,*) dimag(K_j_phi(i,j))
                end do
            end do
            close(77)
            close(78)
            close(79)
            close(80)

        end subroutine

end subroutine