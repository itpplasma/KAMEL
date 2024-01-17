! transform the integral kernels from Fourier basis to spline basis
subroutine basis_transform_kernel(write_out)

    use kernel
    use grid, only: npoib, kr, varphi_lkr, k_space_dim, krp
    use config
    use omp_lib

    implicit none
    logical, intent(in) :: write_out
    integer :: i,j
    integer :: k1, k2
    integer :: max_threads

    double precision :: int_fac1, int_fac2
    integer :: integration_mode = 2

    double complex, dimension(:), allocatable :: kr_int1, kr_int_B1
    double complex, dimension(1) :: res1, res2, res_B1, res_B2

    integer :: nlagr = 4
    integer :: nder = 0
    double precision :: eps = 1d-3
    double precision :: h1 = 0.1d0
    double precision :: hmin = 0.0d0
    integer :: nok, nbad

    max_threads = OMP_GET_MAX_THREADS()
    if (fdebug == 1) write(*,*) ' Debug: Number of threads = ', max_threads
    if (fstatus == 1) write(*,*) 'Status: Transforming basis of kernels, Fourier -> Spline, write_out=',write_out

    allocate(K_rho_phi_llp(npoib, npoib),&
             K_rho_B_llp(npoib, npoib))
    !allocate(kr_int1(npoib, k_space_dim))
    allocate(kr_int1(k_space_dim), kr_int_B1(k_space_dim))
    
    K_rho_phi_llp = 0.1d0
    K_rho_B_llp = 0.0d0

    if (.not. allocated(varphi_lkr)) then
        if (fstatus == 1) write(*,*) 'Status: Calculate Fourier transformed spline functions'
        call calculate_fourier_trans_spline_funcs(.true.)
    end if

    if (integration_mode == 1) then ! trapezoidal
        do i=1, npoib ! l'
            do j=1, npoib ! l
                do k1=1, k_space_dim ! kr
                    do k2 = 1, k_space_dim !kr'

                        if (k1==1 .or. k1==k_space_dim) then
                            int_fac1 = 0.5d0
                        else
                            int_fac1 = 1.0d0
                        end if
                        if (k2==1 .or. k2==k_space_dim) then
                            int_fac2 = 0.5d0
                        else
                            int_fac2 = 1.0d0
                        end if

                        K_rho_phi_llp(i,j) = K_rho_phi_llp(i,j) + int_fac1 * int_fac2 * varphi_lkr(i,k1) *&
                                         conjg(varphi_lkr(j,k2)) * K_rho_phi(k2,k1)
                        K_rho_B_llp(i,j) = K_rho_B_llp(i,j) + int_fac1 * int_fac2 * varphi_lkr(i,k1) *&
                                         conjg(varphi_lkr(j,k2)) * K_rho_B(k2,k1)
                        !if (isnan(real(K_rho_phi_llp(i,j)))) then
                        !   write(*,*) 'NAN: i=', i, ', j=', j, ', k1=', k1, ', k2=', k2
                        !end if
                    end do
                end do
            end do
        end do

        K_rho_phi_llp = K_rho_phi_llp * ((kr(k_space_dim) - kr(1)) / k_space_dim)**2d0
        K_rho_B_llp = K_rho_B_llp * ((kr(k_space_dim) - kr(1)) / k_space_dim)**2d0

    else if (integration_mode == 2) then ! adaptive RK
        if (fstatus==1) write(*,*) 'Status: going into adaptive RK for base transformation'

        !$OMP PARALLEL DO default(none) schedule(guided) &
        !$OMP PRIVATE(i,j,k1, res1, res_B1, res2, res_B2, nok, nbad, kr_int1, kr_int_B1) &
        !$OMP SHARED(k_space_dim, npoib, varphi_lkr, K_rho_phi_llp, K_rho_B_llp, kr, krp, eps,&
        !$OMP h1,hmin, max_threads, fstatus)
        do i = 1, npoib ! l'
            do j = 1, npoib ! l
                ! integrate over kr'
                do k1 = 1, k_space_dim ! kr'
                    res1 = krp(1)
                    res_B1 = krp(1)
                    call integrate_kr_c(res1, size(res1), varphi_lkr(i,:), size(varphi_lkr(i,:)), k1, krp(1), &
                                        krp(k_space_dim), eps, h1, hmin, nok, nbad, integrand_K_rho_phi_krp)
                    !call integrate_kr_c(res_B1, size(res1), varphi_lkr(i,:), size(varphi_lkr(i,:)), k1, krp(1), &
                    !                    krp(k_space_dim), eps, h1, hmin, nok, nbad, integrand_K_rho_B_krp)
                    kr_int1(k1) = res1(1)
                    !kr_int_B1(k1) = res_B1(1)
                end do
                res2 = kr(1)
                !res_B2 = kr(1)
                
                ! integrate over kr
                call integrate_kr_c_final(res2, size(res2), conjg(varphi_lkr(j,:)), size(varphi_lkr(i,:)), &
                                          kr_int1, size(kr_int1), kr(1), kr(k_space_dim), eps, h1, hmin, nok,&
                                          nbad, integrand_K_rho_phi_kr)
                !call integrate_kr_c_final(res_B2, size(res_B2), conjg(varphi_lkr(j,:)), size(varphi_lkr(i,:)), &
                !                          kr_int_B1, size(kr_int_B1), kr(1), kr(k_space_dim), eps, h1, hmin, nok,&
                !                          nbad, integrand_K_rho_B_kr)
 
                !$OMP critical
                K_rho_phi_llp(i,j) = res2(1)
                !K_rho_B_llp(i,j) = res_B2(1)
                !$OMP end critical

            end do
            !!$OMP critical
            !if (fstatus == 1 .and. OMP_GET_THREAD_NUM() == 0) then
            !    write(*,*) ' integration status : ', dble(i) / dble(npoib) * 100.0d0 * max_threads ,'%'
            !end if
            !!$OMP end critical
        end do
        !$OMP END PARALLEL DO
        if (fstatus == 1) write(*,*) 'Status: finished adaptive RK integration of base transformation'

    end if

    if (fstatus == 1) write(*,*) ' Status: finished basis transformation'

    deallocate(K_rho_phi, K_rho_B, kr_int1, kr_int_B1)

    if (write_out) call write_basis_trans_kernel

    contains

    subroutine write_basis_trans_kernel

        use config

        implicit none
        integer :: i,j
        logical :: ex

        inquire(file=trim(output_path)//'kernel', exist=ex)
        if(.not. ex) then
            call system('mkdir -p '//trim(output_path)//'kernel')
        end if

        open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_llp_re.dat')
        open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_llp_im.dat')

        !open(unit=79, file=trim(output_path)//'kernel/K_rho_B_llp_re.dat')
        !open(unit=80, file=trim(output_path)//'kernel/K_rho_B_llp_im.dat')

        do i=1, npoib
            do j=1, npoib
                write(77,*) real(K_rho_phi_llp(i,j))
                write(78,*) dimag(K_rho_phi_llp(i,j))
                !write(79,*) real(K_rho_B_llp(i,j))
                !write(80,*) dimag(K_rho_B_llp(i,j))
            end do
        end do

        close(77)
        close(78)
        close(79)
        close(80)

    end subroutine

    subroutine integrand_K_rho_phi_krp(kr_val, y, dydkr, varphi_l, k_ind)

        implicit none
        integer, intent(in) :: k_ind
        double precision, intent(in) :: kr_val
        double complex, dimension(:), intent(in) :: y
        double complex, dimension(:), intent(in) :: varphi_l
        double complex, dimension(:), intent(out) :: dydkr

        integer :: ibeg, iend , ikrp

        double precision, dimension(:,:), allocatable :: coef

        double complex :: K_intp, vphi_intp

        if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        dydkr = 0.0d0
        K_intp = 0.0d0
        vphi_intp = 0.0d0

        call binsrc(kr_val, 1, k_space_dim, krp, ikrp)
        ibeg = max(1, ikrp - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend .gt. k_space_dim) then
            iend = k_space_dim
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, kr_val, krp(ibeg:iend), coef)

        K_intp = sum(cmplx(coef(0,:)) * K_rho_phi(ibeg:iend, k_ind))
        vphi_intp = sum(cmplx(coef(0,:)) * varphi_l(ibeg:iend))

        dydkr = dydkr + K_intp * vphi_intp
        !deallocate(coef)

    end subroutine


    subroutine integrand_K_rho_phi_kr(kr_val, y, dydkr, varphi_l, kr_int)

        implicit none
        double precision, intent(in) :: kr_val
        double complex, dimension(:), intent(in) :: y
        double complex, dimension(:), intent(in) :: varphi_l
        double complex, dimension(:), intent(in) :: kr_int
        double complex, dimension(:), intent(out) :: dydkr

        integer :: ibeg, iend , ikr

        double precision, dimension(:,:), allocatable :: coef

        double complex :: vphi_intp, K_intp

        dydkr = 0.0d0

        if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        call binsrc(kr_val, 1, k_space_dim, kr, ikr)
        ibeg = max(1, ikr - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend .gt. k_space_dim) then
            iend = k_space_dim
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, kr_val, kr(ibeg:iend), coef)
        vphi_intp = sum(coef(0,:) * varphi_l(ibeg:iend))
        K_intp = sum(coef(0,:) * kr_int(ibeg:iend))

        dydkr = dydkr + K_intp * vphi_intp
        !deallocate(coef)

    end subroutine

    subroutine integrand_K_rho_B_krp(kr_val, y, dydkr, varphi_l, k_ind)

        implicit none
        integer, intent(in) :: k_ind
        double precision, intent(in) :: kr_val
        double complex, dimension(:), intent(in) :: y
        double complex, dimension(:), intent(in) :: varphi_l
        double complex, dimension(:), intent(out) :: dydkr

        integer :: ibeg, iend , ikrp

        double precision, dimension(:,:), allocatable :: coef

        double complex :: K_intp, vphi_intp

        if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        dydkr = 0.0d0
        K_intp = 0.0d0
        vphi_intp = 0.0d0

        call binsrc(kr_val, 1, k_space_dim, krp, ikrp)
        ibeg = max(1, ikrp - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend .gt. k_space_dim) then
            iend = k_space_dim
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, kr_val, krp(ibeg:iend), coef)

        K_intp = sum(cmplx(coef(0,:)) * K_rho_B(ibeg:iend, k_ind))
        vphi_intp = sum(cmplx(coef(0,:)) * varphi_l(ibeg:iend))
        !write(*,*) ' K_intp = ', K_intp, '; vphi_intp = ', vphi_intp

        dydkr = dydkr + K_intp * vphi_intp
        !write(*,*) 'dydkr = ', dydkr
        !deallocate(coef)

    end subroutine


    subroutine integrand_K_rho_B_kr(kr_val, y, dydkr, varphi_l, kr_int)

        implicit none
        double precision, intent(in) :: kr_val
        double complex, dimension(:), intent(in) :: y
        double complex, dimension(:), intent(in) :: varphi_l
        double complex, dimension(:), intent(in) :: kr_int
        double complex, dimension(:), intent(out) :: dydkr

        integer :: ibeg, iend , ikr

        double precision, dimension(:,:), allocatable :: coef

        double complex :: vphi_intp, K_intp

        dydkr = 0.0d0

        if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        call binsrc(kr_val, 1, k_space_dim, kr, ikr)
        ibeg = max(1, ikr - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend .gt. k_space_dim) then
            iend = k_space_dim
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, kr_val, kr(ibeg:iend), coef)

        vphi_intp = sum(coef(0,:) * varphi_l(ibeg:iend))
        K_intp = sum(coef(0,:) * kr_int(ibeg:iend))

        dydkr = dydkr + K_intp * vphi_intp
        !deallocate(coef)

    end subroutine

end subroutine