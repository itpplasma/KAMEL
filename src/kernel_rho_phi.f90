subroutine kernel_rho_phi(write_out)
! fill kernel correlating rho and phi

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants
    use kernel, only: K_rho_phi

    implicit none

    logical, intent(in) :: write_out

    double complex :: besselI ! complex bessel function from bessel.f90
    double complex :: plasma_Z ! plasma dispersion function

    integer :: nder
    double precision, dimension(:,:), allocatable :: coef
    double precision, dimension(:), allocatable :: x

    ! \mathcal{ I }^{m_\phi=0} ... integrand from the zeroth order cyclotron harmonic
    double complex, dimension(:,:,:), allocatable :: calI_m0

    double precision, dimension(:,:,:), allocatable :: a0, a1, a2
    double precision, dimension(:,:,:), allocatable :: bp, bt ! b_+ and b_\times

    integer :: sigma ! for loop over species
    integer :: i,j,n
    double precision :: int_fac

    if (fstatus == 1) write(*,*) 'Status: Generating kernel rho phi, write_out=', write_out

    npoi_der = 4 ! number of polynomials for derivative
    nder = 1 ! first derivative

    allocate(x(npoi_der), coef(0:nder, npoi_der))
    allocate(a0(k_space_dim, k_space_dim, iprof_length))
    allocate(a1(k_space_dim, k_space_dim, iprof_length))
    allocate(a2(k_space_dim, k_space_dim, iprof_length))
    allocate(bp(k_space_dim, k_space_dim, iprof_length))
    allocate(bt(k_space_dim, k_space_dim, iprof_length))
    allocate(calI_m0(k_space_dim, k_space_dim, iprof_length))
    allocate(K_rho_phi(k_space_dim, k_space_dim))

    K_rho_phi = 0.0d0

    ! fill kernel matrix and integrate over radius at the same time.
    ! For the integration, the trapezoidal method is used.
    do i=1, k_space_dim ! kr
        do j = 1, k_space_dim ! krp
            do n = 1, iprof_length ! r_prof
                ! electrons
                bp(i,j,n) = vTe(n)**2 / (2*omce(n)**2) * &
                            (2*ks(n)**2 + kr(i)**2 + krp(j)**2) 
                bt(i,j,n) = vTe(n)**2 / (omce(n)**2) * sqrt(ks(n)**2 + kr(i)**2) *&
                           sqrt(ks(n)**2 + krp(j)**2)

                a0(i,j,n) = e_mass * vTe(n)**2 / omce(n) * exp(-bp(i,j,n)) * &
                            besselI(0, bt(i,j,n), 0) * (- om_E(n) / Te_prof(n) + ks(n)/(e_mass*omce(n)) * &
                            (A1e(n) + (1+bp(i,j,n)) * A2e(n))) + ks(n) * vTe(n)**2/(omce(n)**2) * A2e(n) *&
                            bt(i,j,n) * exp(-bp(i,j,n)) * besselI(-1, bt(i,j,n), 0)
                
                a1(i,j,n) = - e_mass * vTe(n)**2/omce(n) * kp(n) / Te_prof(n) * exp(-bp(i,j,n)) * besselI(0,bt(i,j,n),0)

                a2(i,j,n) = ks(n) / (2 * omce(n)**2) * A2e(n) * exp(-bp(i,j,n)) * besselI(0,bt(i,j,n),0)

                calI_m0(i,j,n) = sqrt(pi)/kp(n) * (plasma_Z(z0e(n)) * (a0(i,j,n) + sqrt(2d0) * &
                                 vTe(n) * z0e(n) * a1(i,j,n) + vTe(n)**2 * a2(i,j,n) * 2 * z0e(n)**2) +&
                                 vTe(n) * sqrt(2d0) * (a1(i,j,n) + vTe(n) * sqrt(2d0) * z0e(n) * a2(i,j,n)))
                
                ! for the trapezoidal integration; first and last summands need factor 1/2
                if (n==1 .or. n==iprof_length) then
                    int_fac = 0.5d0
                else
                    int_fac = 1.0d0
                end if

                K_rho_phi(i,j) = K_rho_phi(i,j) + int_fac * omce(n) /(lambda_De(n)**2d0 * &
                                vTe(n)) * exp(com_unit * (kr(i) -krp(j)) * r_prof(n)) * &
                                (calI_m0(i,j,n) - sqrt(2d0) * vTe(n) / omce(n) * &
                                (exp(-vTe(n)**2/(2*omce(n)**2) * (kr(i) - krp(j))**2) - &
                                exp(-bp(i,j,n)) * besselI(0,bt(i,j,n),0)))

                ! ions
                do sigma=1, ispecies
                    bp(i,j,n) = vTi(sigma, n)**2 / (2*omci(sigma, n)**2) * &
                                (2*ks(n)**2 + kr(i)**2 + krp(j)**2) 
                    bt(i,j,n) = vTi(sigma, n)**2 / (omci(sigma, n)**2) * sqrt(ks(n)**2 + kr(i)**2) *&
                            sqrt(ks(n)**2 + krp(j)**2)

                    a0(i,j,n) = (p_mass * Ai(sigma)) * vTi(sigma, n)**2 / omci(sigma, n) * exp(-bp(i,j,n)) * &
                                besselI(0, bt(i,j,n), 0) * (- om_E(n) / Ti_prof(sigma, n) + &
                                ks(n)/((p_mass* Ai(sigma))*omci(sigma, n)) * (A1i(sigma, n) + (1+bp(i,j,n)) * &
                                A2i(sigma, n))) + ks(n) * vTi(sigma, n)**2/(omci(sigma, n)**2) * A2i(sigma, n) *&
                                bt(i,j,n) * exp(-bp(i,j,n)) * besselI(-1, bt(i,j,n), 0)
                
                    a1(i,j,n) = - (p_mass*Ai(sigma)) * vTi(sigma, n)**2/omci(sigma, n) * kp(n) / &
                                Ti_prof(sigma, n) * exp(-bp(i,j,n)) * besselI(0,bt(i,j,n),0)

                    a2(i,j,n) = ks(n) / (2 * omci(sigma, n)**2) * A2i(sigma, n) * exp(-bp(i,j,n)) * besselI(0,bt(i,j,n),0)

                    calI_m0(i,j,n) = sqrt(pi)/kp(n) * (plasma_Z(z0i(sigma, n)) * (a0(i,j,n) + sqrt(2d0) * &
                                 vTi(sigma, n) * z0i(sigma, n) * a1(i,j,n) + vTi(sigma, n)**2 * a2(i,j,n) * 2 * z0i(sigma, n)**2) +&
                                 vTi(sigma, n) * sqrt(2d0) * (a1(i,j,n) + vTi(sigma, n) * sqrt(2d0) * z0i(sigma, n) * a2(i,j,n)))
                    if (n==1 .or. n==iprof_length) then
                        int_fac = 0.5d0
                    else
                        int_fac = 1.0d0
                    end if

                    K_rho_phi(i,j) = K_rho_phi(i,j) + int_fac * omci(sigma, n) /(lambda_Di(sigma, n)**2d0 &
                                    * vTi(sigma, n)) * exp(com_unit * (kr(i) -krp(j)) * r_prof(n)) * &
                                    (calI_m0(i,j,n) - sqrt(2d0) * vTi(sigma, n) / omci(sigma, n) * &
                                    (exp(-vTi(sigma, n)**2/(2*omci(sigma, n)**2) * (kr(i) - krp(j))**2) &
                                    - exp(-bp(i,j,n)) * besselI(0,bt(i,j,n),0)))
                end do
            end do
        end do
    end do 

    K_rho_phi = K_rho_phi / (2**(7d0/2d0) * pi**(5d0/2d0)) * (r_prof(size(r_prof)) - r_prof(1))/iprof_length

    if (write_out) call write_kernel_rho_phi

    contains

        subroutine write_kernel_rho_phi
            implicit none
            integer :: i,j
            logical :: ex

            inquire(file=trim(output_path)//'kernel', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'kernel')
            end if
            open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_kr_re.dat')
            open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_kr_im.dat')
            do i=1,k_space_dim
                do j=1,k_space_dim
                    write(77,*) real(K_rho_phi(i,j))
                    write(78,*) dimag(K_rho_phi(i,j))
                end do
            end do
            close(77)
            close(78)

        end subroutine

end subroutine